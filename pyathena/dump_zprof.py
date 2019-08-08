import pandas as pd
import pyathena as pa
import pickle as p
import os
import numpy as np
import xarray as xr

from pyathena import ath_hst

import glob
import argparse

# astropy
import astropy.constants as c
import astropy.units as u

solar_param={'surf_s':42*u.M_sun/u.pc**2,
            'rho_dm':0.0064*u.Msun/u.pc**3,
            'z0': 0.245*u.kpc,
            'R0': 8*u.kpc,}

class tigress_extpot(object):
    def __init__(self,param=solar_param):
        self.surf_s=param['surf_s']
        self.rho_dm=param['rho_dm']
        self.z0=param['z0']
        self.R0=param['R0']
        
    def gext(self,z):

        a1=2*np.pi*c.G*self.surf_s
        a2=4*np.pi*c.G*self.rho_dm
        g1=-a1*z/np.sqrt(z**2+self.z0**2)
        g2=-a2*z/(z**2/self.R0**2+1)
        g_new=g1+g2
    
        return g_new

    def phiext(self,z):
        phi=2*np.pi*c.G*(self.surf_s*(np.sqrt(z**2+self.z0**2)-self.z0)
                         +self.rho_dm*self.R0**2*np.log(1+z**2/self.R0**2))
        return phi
    
    def vesc(self,z):
        return np.sqrt(2*(self.phiext(z.max()).to('km^2/s^2')-self.phiext(z).to('km^2/s^2')))

def parse_par(rstfile):

    fp=open(rstfile,'rb')
    par={}
    line=fp.readline().decode('utf-8')

    while 1:

        if line.startswith('<'):
            block=line[1:line.rfind('>')]
            if block == 'par_end': break
            par[block]={}
        line=fp.readline().decode('utf-8')

        if block in ['problem','domain1','time']:
            sp = line.strip().split()
            if len(sp) >= 3: par[block][sp[0]]=eval(sp[2])
        else:
            sp=line.split('=')
            if len(sp) == 2: par[block][sp[0].strip()]=sp[1].split('#')[0].strip()

    par[block]=fp.tell()

    fp.close()

    return par

def set_field_list(par):
    shearing_box,cooling,MHD,nscalar,Gamma = get_configure(par)
    zprof_flist=['A','d','v1','v2','v3','M1','M2','M3','Ek1','Ek2','Ek3','P','T']
    if shearing_box: zprof_flist+=['dM2','dEk2']
    zprof_flist+=['Phie','gext','dWext','Phisg','gsg','dWsg','Ber']
    if cooling: zprof_flist+=['cool','heat']
    if MHD: 
        zprof_flist+='B1,B2,B3,PB1,PB2,PB3,vA1,vA2,vA3'.split(',')
        zprof_flist+='dB1,dB2,dB3,dPB1,dPB2,dPB3,dvA1,dvA2,dvA3'.split(',')
        zprof_flist+='S1,S2,S3'.split(',')
    for i in range(nscalar): zprof_flist+=['s%d' % (i+1)]
    zprof_flist+='pA,pd,pvz,pFzd,pFzM1,pFzM2,pFzM3,pFzE1,pFzE2,pFzE3,pFzP'.split(',')
    zprof_flist+='pFzEge,pFzEgsg,pFzEtidal'.split(',')
    if MHD: zprof_flist+='pSzEm1,pSzEm2,pSzvB1,pSzvB2'.split(',')
    for i in range(nscalar): zprof_flist+=['pFzs%d' % (i+1)]
    zprof_flist+='mA,md,mvz,mFzd,mFzM1,mFzM2,mFzM3,mFzE1,mFzE2,mFzE3,mFzP'.split(',')
    zprof_flist+='mFzEge,mFzEgsg,mFzEtidal'.split(',')
    if MHD: zprof_flist+='mSzEm1,mSzEm2,mSzvB1,mSzvB2'.split(',')
    for i in range(nscalar): zprof_flist+=['mFzs%d' % (i+1)]
    if shearing_box:
        zprof_flist+=['RxyL','RxyR']
        if MHD: zprof_flist+=['MxyL','MxyR']
            
    return zprof_flist

def get_configure(par):
    shearing_box=False
    if par['configure']['ShearingBox'] == 'yes': shearing_box=True 
    cooling=False
    if par['configure']['cooling'] == 'ON': cooling=True
    MHD=False
    if par['configure']['gas'] == 'mhd': MHD=True
    nscalar = eval(par['configure']['nscalars'])
    Gamma=par['problem']['gamma']
    return shearing_box,cooling,MHD,nscalar,Gamma

def get_zprof(data,domain,par,hst):
    Omega=par['problem']['Omega']
    qshear=par['problem']['qshear']
    shearing_box,cooling,MHD,nscalar,Gamma = get_configure(par)
    import pyathena as pa
    x,y,z,=pa.cc_arr(domain)
    Nx=domain['Nx']
    x3d=np.tile(x.reshape(1,1,Nx[0]),(Nx[2],Nx[1],1))
    z3d=np.tile(z.reshape(Nx[2],1,1),(1,Nx[1],Nx[0]))
    if shearing_box: 
        vy0=-qshear*Omega*x3d
        Phit=-qshear*Omega**2*x3d**2
    else:
        vy0=np.zeros(x3d.shape)
    
    if cooling: 
        heat_ratio=np.interp(domain['time'],hst['time'],hst['heat_ratio'])
        unit=pa.set_units(muH=1.4271)
        unitC=unit['density']*unit['velocity']**3/unit['length']
        unitC=unitC.cgs.value
        coolftn=pa.coolftn()
    
    dA=domain['dx'][0]*domain['dx'][1]
    dz=domain['dx'][2]
    
    d=np.copy(data['density'].data)
    v1=np.copy(data['velocity1'].data)
    v2=np.copy(data['velocity2'].data)
    v3=np.copy(data['velocity3'].data)
    dv2=v2-vy0
    if MHD:
        B1=np.copy(data['magnetic_field1'].data)
        B2=np.copy(data['magnetic_field2'].data)
        B3=np.copy(data['magnetic_field3'].data)
    P=np.copy(data['pressure'].data)
    Phi=np.copy(data['gravitational_potential'].data)
    T1=np.copy(data['T1'].data)
    
    if cooling:
        lam=coolftn.get_cool(T1)
        gam=coolftn.get_heat(T1)*heat_ratio

    zprof={}

    zprof['d']=d
    zprof['v1']=v1
    zprof['v2']=v2
    zprof['v3']=v3
    zprof['P']=P
    if MHD:
        zprof['B1']=B1
        zprof['B2']=B2
        zprof['B3']=B3
        dsqrt=np.sqrt(d)
        vdotB=B1*v1+B2*v2+B3*v3
        Emag=0.5*(B1**2+B2**2+B3**2)
    meanB={'1':B1.mean(axis=2).mean(axis=1),'2':B2.mean(axis=2).mean(axis=1),'3':B3.mean(axis=2).mean(axis=1)}
    for vec in ['1','2','3']:
        zprof['M%s' % vec] = zprof['d']*zprof['v%s' % vec]
        zprof['Ek%s' % vec] = 0.5*zprof['d']*zprof['v%s' % vec]**2
        if MHD:
            zprof['PB%s' % vec] = 0.5*zprof['B%s' % vec]**2
            zprof['vA%s' % vec] = zprof['B%s' % vec]/dsqrt
            Bbar=np.tile(meanB[vec].reshape(Nx[2],1,1),(1,Nx[1],Nx[0])).astype('double')
            zprof['dB%s' % vec] = zprof['B%s' % vec] - Bbar
            zprof['dPB%s' % vec] = 0.5*zprof['dB%s' % vec]**2
            zprof['dvA%s' % vec] = zprof['dB%s' % vec]/dsqrt
            zprof['S%s' % vec] = 2.0*Emag*zprof['v%s' % vec] - zprof['B%s' % vec]*vdotB
    if shearing_box:
        zprof['v2']=dv2
        zprof['dM2'] = zprof['d']*zprof['v2']
        zprof['dEk2'] = 0.5*zprof['d']*zprof['v2']**2
#        zprof['Phit'] = Phit
    zprof['T']=coolftn.get_temp(T1)
    if cooling:
        zprof['cool']=d**2*lam/unitC
        zprof['heat']=d*gam/unitC
    for ns in range(nscalar):
        zprof['s%s' % (ns +1)]=data['specific_scalar%s' % ns].data*d
        
    param={}
    param['R0']=par['problem']['R0']*c.pc
    param['rho_dm']=par['problem']['rhodm']*c.M_sun/c.pc**3
    param['surf_s']=par['problem']['SurfS']*c.M_sun/c.pc**2
    param['z0']=par['problem']['zstar']*c.pc
    phiext=tigress_extpot(param).phiext
    Phie=phiext(z3d*c.pc).to('km**2/s**2').value
    gext=phiext((z3d+dz/2)*c.pc)-phiext((z3d-dz/2)*c.pc)
    zprof['Phie']=Phie
    zprof['gext']=gext.to('km**2/s**2').value/dz
    zprof['dWext']=d*zprof['gext']
    zprof['Phisg']=Phi
    dPhi=np.copy(Phi)
    dPhi[1:-1,:,:]=(Phi[2:,:,:]-np.roll(Phi,2,axis=0)[2:,:,:])/2.0/dz
    dPhi[0,:,:]=(Phi[1,:,:]-Phi[0,:,:])/dz
    dPhi[-1,:,:]=(Phi[-1,:,:]-Phi[-2,:,:])/dz
    zprof['gsg']=dPhi
    zprof['dWsg']=d*zprof['gsg']
    zprof['Ber'] = 0.5*(v1**2+v2**2+v3**2) + Gamma/(Gamma-1)*P/d+Phie+Phi+Phit

    for pm in ['p','m']:
        zprof[pm+'A'] = np.ones(d.shape)*dA
        zprof[pm+'d'] = np.copy(d)
        zprof[pm+'vz'] = np.copy(v3)
        for f in ['d','M1','M2','M3']:
            zprof['%sFz%s' % (pm,f)] = zprof[f]*v3
        for f in ['E1','E2','E3','Ege','Egsg','Etidal']:
            if f in ['E1','E2','E3']: 
                zf='%sk%s' %(f[0],f[1])
                tmp = zprof[zf]*v3
            elif f == 'Ege': tmp = d*Phie*v3
            elif f == 'Egsg': tmp = d*Phi*v3
            elif f == 'Etidal': tmp = d*Phit*v3
            zprof['%sFz%s' % (pm,f)] = tmp
        zprof['%sFzP' % pm] = Gamma/(Gamma-1)*zprof['P']*v3
        if MHD:
            zprof['%sSzEm1' % pm] = 2.0*zprof['PB1']*v3
            zprof['%sSzEm2' % pm] = 2.0*zprof['PB2']*v3
            zprof['%sSzvB1' % pm] = -B3*B1*v1
            zprof['%sSzvB2' % pm] = -B3*B2*v2
        for ns in range(nscalar):
            zprof['%sFzs%s' % (pm,ns+1)]=data['specific_scalar%s' % ns].data*d*v3
    if shearing_box:
        zprof['Rxy']=d*v1*dv2
        zprof['Mxy']=-B1*B2
        
    nv3=zprof['v3']<0
    pv3=~nv3
    for k in zprof:
        if k.startswith('p'): zprof[k][nv3]=0
        if k.startswith('m'): zprof[k][pv3]=0

    zprof['A']=np.ones(d.shape)*dA

    return z,zprof

def get_phase(data):
    temp=data['T']

    idx={}
    idx['phase1']=temp < 184
    idx['phase2']=(temp >= 184) * (temp <5050)
    idx['phase3']=(temp >= 5050) * (temp <2.e4)
    idx['phase4']=(temp >= 2.e4) * (temp <5.e5)
    idx['phase5']=temp >= 5.e5
    return idx

def get_mean(dset,dA):
    
    zp=dset.sum(dim=['x','y'])*dA

    zp['A']=dset['A'].sum(dim=['x','y'])

    Nx=len(dset.x)

    zp=zp.drop('Rxy')
    zp=zp.drop('Mxy')

    stress=dset['Rxy'].sum(axis=1)
    zp['RxyL']=stress.isel(x=0).drop('x')*Nx*dA
    zp['RxyR']=stress.isel(x=-1).drop('x')*Nx*dA
    if 'Mxy' in dset:
        stress=dset['Mxy'].sum(axis=1)
        zp['MxyL']=stress.isel(x=0).drop('x')*Nx*dA
        zp['MxyR']=stress.isel(x=-1).drop('x')*Nx*dA
    return zp

def get_full_domain(par):
    full_domain={}
    full_domain['left_edge']=np.array([par['domain1']['x1min'],par['domain1']['x2min'],par['domain1']['x3min']]).astype('float')
    full_domain['right_edge']=np.array([par['domain1']['x1max'],par['domain1']['x2max'],par['domain1']['x3max']]).astype('float')
    full_domain['Nx']=np.array([par['domain1']['Nx1'],par['domain1']['Nx2'],par['domain1']['Nx3']])
    full_domain['Lx']=full_domain['right_edge']-full_domain['left_edge']
    full_domain['dx']=full_domain['Lx']/full_domain['Nx']
    return full_domain

def dump_zprof_one(f,dc,icm_field=None,outdir='zprof_icm',
                   plist=['phase1','phase2','phase3','phase4','phase5']):
    par=dc.par
    hst=dc.hst

    flist=set_field_list(par)
    full_domain=get_full_domain(par)
    x,y,z,=pa.cc_arr(full_domain)

    if icm_field is not None: icm=True

    print('Reading: ',f)

    ds = pa.AthenaDataSet(f)
    time = ds.domain['time']

    zp_slab=[]
    if icm: zp_slab_icm=[]
    print('Calculating zprof...') 
    for islab in range(ds.NGrids[2]):
        print('{}/{}'.format(islab,ds.NGrids[2]),end=' ')
        grids=ds._get_slab_grid(islab+1,verbose=False)
        slab_domain=ds._setup_domain(grids)
        xs,ys,zs,=pa.cc_arr(slab_domain)
        data=xr.Dataset()
        for fi in ds.field_list+['T1']:
            slab_data=ds.read_all_data(fi,slab=islab+1)
            if slab_data.ndim == 4:
                for ivec in [1,2,3]:
                    data['{}{}'.format(fi,ivec)]=xr.DataArray(slab_data[...,ivec-1],
                                                              dims=['z','y','x'],
                                                              coords=[zs,ys,xs])
            else:
                data[fi]=xr.DataArray(slab_data,dims=['z','y','x'],coords=[zs,ys,xs])

        z_part,zpw_part=get_zprof(data,slab_domain,par,hst)
        dset=xr.Dataset()
        for field in zpw_part:
            dset[field]=xr.DataArray(zpw_part[field],
                                     dims=['z','y','x'],coords=[zs,ys,xs])

        idx=get_phase(dset)
        dA=float(dset['A'][0,0,0].data)
        zp_part=[]
        if icm: 
            ficm=np.clip(dset[icm_field]/dset['d'],0,1)
            zp_part_icm=[]

        for phase in plist:
            zidx_part=pd.Series(z_part,name='z')

            zp_part.append(get_mean(dset*idx[phase],dA).expand_dims('phase'))

            if icm: 
                zp_part_icm.append(get_mean(dset*idx[phase]*ficm,dA).expand_dims('phase'))
        zp_slab.append(xr.concat(zp_part,dim='phase'))
        if icm: zp_slab_icm.append(xr.concat(zp_part_icm,dim='phase'))

        for g in grids:
            g['data'].clear()

    zpds=xr.concat(zp_slab,dim='z')
    zpds.coords['phase']=plist
    if icm: 
        zpds_icm=xr.concat(zp_slab_icm,dim='z')
        zpds_icm.coords['phase']=plist

    print('\nWriting at {}: '.format(os.path.dirname(f).replace('id0',outdir)))
    for phase in plist:
        print(phase,end=' ')
        zprof_fname=f.replace('vtk','%s.zprof' % phase).replace('id0',outdir)
        with open(zprof_fname,'w') as fp:
            fp.write('# Athena vertical profile at t={}\n'.format(time))
            zpdf=zpds.sel(phase=phase).drop('phase').to_dataframe()
            zpdf.to_csv(fp)
        if icm: 
            zprof_fname_icm=f.replace('vtk','%s-icm.zprof' % phase).replace('id0',outdir)
            with open(zprof_fname_icm,'w') as fp:
                fp.write('# Athena vertical profile at t={}\n'.format(time))
                zpdf=zpds_icm.sel(phase=phase).drop('phase').to_dataframe()
                zpdf.to_csv(fp)

    if icm: 
        return zpds,zpds_icm
    else:
        return zpds

class data_container(object):
    def __init__(self,pid,base='/tigress/changgoo/',pdir=None):
        if pdir is None: pdir=pid+'/'
        self.vtkfiles= glob.glob(base+pdir+"id0/"+pid+".????.vtk")
        self.vtkfiles.sort() 
        self.parfname='%s%s/%s.par' % (base,pdir,pid)
        self.par=parse_par(self.parfname)
        self.hstfname='%s%s/hst/%s.hst' % (base,pdir,pid)
        self.hst=pa.hst_reader(self.hstfname)
