
import xarray as xr
import pandas as pd
import os,glob
import shutil


from pyathena.parse_par import write_par_from_rst
import pyathena as pa
from pyathena.yt_analysis import yt_analysis

import astropy.constants as c
import astropy.units as u
import numpy as np

def cleanup_directory(base,problem_id,problem_dir=None):
    """
       move files into separate directories
       create par files from rst or out
    """
    if problem_dir is None: problem_dir = problem_id

    # parse parameters from rstfiles
    parfile='{}{}/{}.par'.format(base,problem_dir,problem_id)
    success={'par':True}
    if os.path.isfile(parfile):
        success['par']=False
    else:
        rstfiles=glob.glob('{}{}/id0/{}.????.rst'.format(base,problem_dir,problem_id))
        rstfiles+=glob.glob('{}{}/rst/{}.????.rst'.format(base,problem_dir,problem_id))
        if len(rstfiles):
            write_par_from_rst(rstfiles[0],parfile)
        else:
            print('cannot find restart file...')
            success['par']=False

    # move files to their directory
    nprocs=len(glob.glob('{}{}/id*'.format(base,problem_dir)))
    for ext in ['zprof','hst','rst','starpar']:
        subdir='{}{}/{}'.format(base,problem_dir,ext)
        success[ext]=True
        if not os.path.isdir(subdir): os.mkdir(subdir)
        if ext is 'zprof':
            files=glob.glob('{}{}/id0/{}.*.zprof'.format(base,problem_dir,problem_id))
            files+=glob.glob('{}{}/id0/{}.*.zprof.p'.format(base,problem_dir,problem_id))
        elif ext is 'hst':
            files=glob.glob('{}{}/id0/{}.hst'.format(base,problem_dir,problem_id))
            files+=glob.glob('{}{}/id0/{}.sn'.format(base,problem_dir,problem_id))
        elif ext is 'rst':
            files=glob.glob('{}{}/id0/{}.????.rst'.format(base,problem_dir,problem_id))
        elif ext is 'starpar':
            files=glob.glob('{}{}/id0/{}.????.starpar.vtk'.format(base,problem_dir,problem_id))
        else:
            print('{} extention specifier is not known'.format(ext))
            success[ext]=False
        if len(files):
            print('moving {} files...'.format(ext))
            for f in files:
                if ext is 'hst':
                    shutil.copy2(f,f.replace('id0/','{}/'.format(ext)))
                else:
                    shutil.move(f,f.replace('id0/','{}/'.format(ext)))
        else:
            success[ext]=False

        if ext is 'rst':
            for i in range(1,nprocs):
                files=glob.glob('{}{}/id{}/{}-id{}.????.rst'.format(base,problem_dir,i,problem_id,i))
                if len(files): 
                    for f in files: shutil.move(f,f.replace('id{}/'.format(i),'rst/'))

    return success

def zprof_to_xarray(base,problem_id):
    """
        merge zprof dumps along t-axis and create netCDF files 
        files should have moved to zprof/
    """    
    zpmerge_dir='{}{}/zprof_merged/'.format(base,problem_id)
    if not os.path.dirname(zpmerge_dir): os.mkdir(zpmerge_dir)
    plist=['phase1','phase2','phase3','phase4','phase5','whole']
    for phase in plist:
        zprof_fnames=glob.glob('{}{}/zprof/{}.*.{}.zprof'.format(base,problem_id,problem_id,phase))
        zprof_fnames.sort()

        taxis=[]
        dfall=None
        for f in zprof_fnames:
            fp=open(f,'r')
            hd=fp.readline()
            fp.close()
            time=float(hd[hd.rfind('t=')+2:])
            df=pd.read_csv(f,skiprows=1)
            zaxis=np.array(df['z'])
            fields=np.array(df.columns)
            taxis.append(time)
            if dfall is None:
                dfall=np.array(df)[np.newaxis,:]
            else:
                dfall=np.concatenate([dfall,np.array(df)[np.newaxis,:]],axis=0)
    
        da=xr.DataArray(dfall.T,coords={'fields':fields,'zaxis':zaxis,'taxis':taxis},
                                dims=('fields','zaxis','taxis'))
        zpfile='{}{}/zprof_merged/{}.{}.zprof.nc'.format(base,problem_id,problem_id,phase)
        da.to_netcdf(zpfile)
        print('{} is created'.format(zpfile))

def merge_xarray(base,problem_id):
    """
        merge all phases into a single file
    """    
    plist=['phase1','phase2','phase3','phase4','phase5']
    datasets = xr.Dataset()
    for phase in plist:
        path='{}{}/zprof_merged/{}.{}.zprof.nc'.format(base,problem_id,problem_id,phase)
        with xr.open_dataarray(path) as da: da.load()
        datasets[phase]=da
    return datasets
    #zpfile='{}{}/zprof_merged/{}.zprof.nc'.format(base,problem_id,problem_id)
    #datasets.to_netcdf(zpfile)
    #print('{} is created'.format(zpfile))

def panel_to_xarray(base,problem_id):
    """
        convert exsiting merged pickles for pandas Panel to netCDF for xarray DataArray
    """
    zpfile='{}{}/zprof_merged/{}.whole.zprof.p'.format(base,problem_id,problem_id)
    if os.path.isfile(zpfile):
        print('no merged pickle files')
        return

    if not os.path.isfile(zpfile.replace('.zprof.p','.zprof.nc')):
        print('merged pickle files are found and begin converting to netCDF')
        pn=pd.read_pickle(open(zpfile,'rb'))

        plist=['phase1','phase2','phase3','phase4','phase5','whole']
        for ph in plist:
            zpdata=pd.read_pickle(open(zpfile.replace('whole',ph),'rb'))
            taxis=zpdata.minor_axis.values
            zaxis=zpdata.major_axis.values
            fields=zpdata.items.values
            
            da=xr.DataArray(zpdata.values,coords={'fields':fields,'zaxis':zaxis,'taxis':taxis},
                                dims=('fields','zaxis','taxis'))
            ncfile=zpfile.replace('whole.zprof.p','{}.zprof.nc'.format(ph))
            da.to_netcdf(ncfile)
            print('{} is created'.format(ncfile))

def recal_history(base,problem_id):
    parfile='{}{}/{}.par'.format(base,problem_id,problem_id)
    params=pa.get_params(parfile)

    hstfile='{}{}/hst/{}.hst'.format(base,problem_id,problem_id)
    hst=pa.hst_reader(hstfile)

    hstrecalfile='{}_cal.p'.format(hstfile)
    h=processing_history_dump(hst,params,hstrecalfile)

    hstzpfile='{}_zp.p'.format(hstfile)
    if not os.path.isfile(hstzpfile): 
        zprof_ds=merge_xarray(base,problem_id)
        print("zprof dataset is loaded")
        h_zp=processing_zprof_dump(h,params,zprof_ds,hstzpfile)

def processing_history_dump(hst,params,hstfile):
    unit=pa.set_units(muH=1.4271)
    toMsun=unit['mass'].to('Msun').value
    toMyr=unit['time'].to('Myr').value
    toflux=(unit['density']*unit['velocity']).to('Msun/pc^2/Myr').value
    toPok=(unit['pressure']/c.k_B).cgs.value

    Lz=params['x3max']-params['x3min']
    Ly=params['x2max']-params['x2min']
    Lx=params['x1max']-params['x1min']
    Nx=params['Nx1']
    Ny=params['Nx2']
    Nz=params['Nx3']
    Ntot=Nx*Ny*Nz
    vol=Lx*Ly*Lz
    dz=Lz/Nz
    Omega=params['Omega']
    torb=2*np.pi/Omega*toMyr

    if 'x1ME' in hst: mhd=True

    # process history dump
    h=pd.DataFrame()

    h['tMyr']=hst['time']*toMyr
    h['torb']=h['tMyr']/torb
    h['surf']=hst['mass']*toMsun*Lz
    h['surfsp']=hst['msp']*toMsun*Lz
    for ph in ['c','u','w','h1','h2']:
        h['mf_{}'.format(ph)]=hst['M{}'.format(ph)]/hst['mass']
        h['vf_{}'.format(ph)]=hst['V{}'.format(ph)]
        h['H_{}'.format(ph)]=np.sqrt(hst['H2{}'.format(ph)]/hst['M{}'.format(ph)])
    h['H']=np.sqrt(hst['H2']/hst['mass'])

    h['mf_2p']=h['mf_c']+h['mf_u']+h['mf_w']
    h['vf_2p']=h['vf_c']+h['vf_u']+h['vf_w']
    h['H_2p']=np.sqrt((hst['H2c']+hst['H2u']+hst['H2w'])/(hst['Mc']+hst['Mu']+hst['Mw']))

    h['KE']=hst['x1KE']+hst['x2KE']+hst['x3KE']
    if mhd: h['ME']=hst['x1ME']+hst['x2ME']+hst['x3ME']
    hst['x2KE']=hst['x2dke']
    for ax in ['1','2','3']:
        Ekf = 'x{}KE'.format(ax)
        if ax == '2': Ekf = 'x2dke'
        h['v{}'.format(ax)]=np.sqrt(2*hst[Ekf]/hst['mass'])
        if mhd: h['vA{}'.format(ax)]=np.sqrt(2*hst['x{}ME'.format(ax)]/hst['mass'])
        h['v{}_2p'.format(ax)]=np.sqrt(2*hst['x{}KE_2p'.format(ax)]/hst['mass']/h['mf_2p'])
    h['cs']=np.sqrt(hst['P']/hst['mass'])
    h['Pth_mid']=hst['Pth']*toPok
    h['Pth_mid_2p']=hst['Pth_2p']*toPok/hst['Vmid_2p']
    h['Pturb_mid']=hst['Pturb']*toPok
    h['Pturb_mid_2p']=hst['Pturb_2p']*toPok/hst['Vmid_2p']
    h['nmid']=hst['nmid']
    h['nmid_2p']=hst['nmid_2p']/hst['Vmid_2p']
    h['sfr10']=hst['sfr10']
    h['sfr40']=hst['sfr40']
    h['sfr100']=hst['sfr100']
    h.index=hst['time']
    h.to_pickle(hstfile)
    print('{} is created'.format(hstfile))

    return h

def processing_zprof_dump(h,params,zprof_ds,hstfile):
    unit=pa.set_units(muH=1.4271)
    toMsun=unit['mass'].to('Msun').value
    toMyr=unit['time'].to('Myr').value
    toflux=(unit['density']*unit['velocity']).to('Msun/pc^2/Myr').value
    toPok=(unit['pressure']/c.k_B).cgs.value

    Lz=params['x3max']-params['x3min']
    Ly=params['x2max']-params['x2min']
    Lx=params['x1max']-params['x1min']
    Nx=params['Nx1']
    Ny=params['Nx2']
    Nz=params['Nx3']
    Ntot=Nx*Ny*Nz
    vol=Lx*Ly*Lz
    dz=Lz/Nz
    Omega=params['Omega']
    torb=2*np.pi/Omega*toMyr

    if 'vA1' in h: mhd=True

    # process zprof dump
    h_zp=pd.DataFrame()
    zpw=zprof_ds.to_array().sum(dim='variable')
    zpw['z']=zpw.zaxis
    zp2p=zprof_ds.to_array().sel(variable=['phase1','phase2','phase3']).sum(dim='variable')
    zp2p['z']=zp2p.zaxis

    h_zp['tMyr']=zpw.taxis*toMyr
    h_zp['torb']=h_zp['tMyr']/torb

    for i,ph in enumerate(['_c','_u','_w','_h1','_h2','_2p','']):
        dtot=zpw.sel(fields='d').sum(dim='zaxis')
        Atot=zpw.sel(fields='A').sum(dim='zaxis')
        if ph is '_2p':
            zp=zp2p
        elif ph is '':
            zp=zpw
        else:
            zp=zprof_ds['phase{}'.format(i+1)]
        dphase=zp.sel(fields='d').sum(dim='zaxis')
        h_zp['surf{}'.format(ph)]=(zp.sel(fields='d')/zpw.sel(fields='A')).sum(dim='zaxis')*toMsun*dz
        if ph is not '':
            h_zp['mf{}'.format(ph)]=zp.sel(fields='d').sum(dim='zaxis')/dtot
            h_zp['vf{}'.format(ph)]=zp.sel(fields='A').sum(dim='zaxis')/Atot
        h_zp['H{}'.format(ph)]=np.sqrt((zp.sel(fields='d')*zp.zaxis**2).sum(dim='zaxis')/dphase)
        for ax in ['1','2','3']:
            Ekf = 'Ek{}'.format(ax)
            if ax == '2': Ekf = 'dEk2'
            h_zp['v{}{}'.format(ax,ph)] = np.sqrt(2.0*zp.sel(fields=Ekf).sum(dim='zaxis')/dphase)
            if mhd:
                h_zp['vA{}{}'.format(ax,ph)] = np.sqrt(2.0*zp.sel(fields='PB{}'.format(ax)).sum(dim='zaxis')/dphase)
                h_zp['dvA{}{}'.format(ax,ph)] = np.sqrt(2.0*zp.sel(fields='dPB{}'.format(ax)).sum(dim='zaxis')/dphase)
        h_zp['cs{}'.format(ph)] = np.sqrt(zp.sel(fields='P').sum(dim='zaxis')/dphase)
        h_zp['sigma_eff{}'.format(ph)] = h_zp['cs{}'.format(ph)]**2 + h_zp['v3{}'.format(ph)]**2
        if mhd:
            h_zp['sigma_eff{}'.format(ph)] += 0.5*(h_zp['vA1{}'.format(ph)]**2 + 
                                                   h_zp['vA2{}'.format(ph)]**2 - 
                                                   h_zp['vA3{}'.format(ph)]**2)
        h_zp['sigma_eff{}'.format(ph)] = np.sqrt(h_zp['sigma_eff{}'.format(ph)])
        mid_idx=np.abs(zp.zaxis) < dz
        zpmid=zp[:,mid_idx,:]
        h_zp['vf_mid{}'.format(ph)] = zpmid.loc['A'].sum(axis=0)/zpw[:,mid_idx,:].loc['A'].sum(axis=0)
        h_zp['nmid{}'.format(ph)] = zpmid.loc['d'].sum(axis=0)/zpmid.loc['A'].sum(axis=0)
        h_zp['Pth_mid{}'.format(ph)] = (zpmid.loc['P']/zpmid.loc['A']).mean(axis=0)*toPok
        h_zp['Pturb_mid{}'.format(ph)] = (2.0*zpmid.loc['Ek3']/zpmid.loc['A']).mean(axis=0)*toPok
        h_zp['Pmid{}'.format(ph)] = h_zp['Pth_mid{}'.format(ph)]+h_zp['Pturb_mid{}'.format(ph)]
        if mhd:
            oPmag = 0.5*(zpmid.loc['B1']/zpmid.loc['A'])**2
            oPmag += 0.5*(zpmid.loc['B2']/zpmid.loc['A'])**2
            oPmag -= 0.5*(zpmid.loc['B3']/zpmid.loc['A'])**2
            h_zp['oPmag_mid{}'.format(ph)] = oPmag.mean(axis=0)*toPok
            tPmag = (zpmid.loc['dPB1']/zpmid.loc['A'])
            tPmag += (zpmid.loc['dPB2']/zpmid.loc['A'])
            tPmag -= (zpmid.loc['dPB3']/zpmid.loc['A'])
            h_zp['tPmag_mid{}'.format(ph)] = tPmag.mean(axis=0)*toPok
            h_zp['Pmid{}'.format(ph)] += h_zp['oPmag_mid{}'.format(ph)]+h_zp['tPmag_mid{}'.format(ph)]


        dWext=zp.sel(fields='dWext')/zpw.sel(fields='A')
        dWsg=zp.sel(fields='dWsg')/zpw.sel(fields='A')
        Wext1=-dWext.cumsum(dim='zaxis')
        Wext2=dWext[::-1,:].cumsum(dim='zaxis')
        Wsg1=-dWsg.cumsum(dim='zaxis')
        Wsg2=dWsg[::-1,:].cumsum(dim='zaxis')
        h_zp['Wext_upper{}'.format(ph)]=Wext1.max(dim='zaxis')*toPok*dz
        h_zp['Wext_lower{}'.format(ph)]=Wext2.max(dim='zaxis')*toPok*dz
        h_zp['Wsg_upper{}'.format(ph)]=Wsg1.max(dim='zaxis')*toPok*dz
        h_zp['Wsg_lower{}'.format(ph)]=Wsg2.max(dim='zaxis')*toPok*dz
        h_zp['Wext{}'.format(ph)]=0.5*(h_zp['Wext_upper{}'.format(ph)]+h_zp['Wext_lower{}'.format(ph)])
        h_zp['Wsg{}'.format(ph)]=0.5*(h_zp['Wsg_upper{}'.format(ph)]+h_zp['Wsg_lower{}'.format(ph)])
        h_zp['W{}'.format(ph)]=h_zp['Wext{}'.format(ph)]+h_zp['Wsg{}'.format(ph)]

    # data intepolated from history dump
    for field in ['sfr10','sfr40','sfr100','surfsp']:
        h_zp[field] = np.interp(zpw.taxis,h.index,h[field])

    # calculate zeta
    H=h_zp['surf']/toMsun/h_zp['nmid']/2.0
    zH=np.abs(zpw.zaxis.data[:,np.newaxis]-H[np.newaxis,:])
    zidx=zH.argmin(axis=0)
    gext=zpw.sel(fields='gext')/zpw.sel(fields='A')
    gext=gext[:,0]
    gextH=gext.isel(zaxis=zidx).data
    dWext=zpw.sel(fields='dWext')/zpw.sel(fields='A')
    Wext=-dWext.cumsum(dim='zaxis')
    Wext=Wext.max(dim='zaxis').data*dz
    surf=zpw.sel(fields='d')/zpw.sel(fields='A')
    surf=surf.sum(dim='zaxis').data*dz
    h_zp['zetad']=Wext/gextH/surf

    # calculate PDE
    unit_surf=c.M_sun/c.pc**2
    unit_rho=c.M_sun/c.pc**3
    unit_v=unit['velocity']
    Sigma_gas=h_zp['surf']
    Sigma_sp=h_zp['surfsp']
    Sigma_star=params['SurfS']
    H_star=params['zstar']
    rho_dm=params['rhodm']
    rho_sd=rho_dm+Sigma_star/H_star/2.0
    sigma_z=h_zp['sigma_eff_2p']
    h_zp['PDE_gas']=np.pi*Sigma_gas**2/2.0
    h_zp['PDE_gas']*=(c.G*unit_surf**2/c.k_B).to('K/cm**3').value
    h_zp['chi']=4*h_zp['zetad']*rho_sd/h_zp['nmid']*(unit_rho/unit['density']).cgs
    h_zp['PDE_ext']=h_zp['PDE_gas']*h_zp['chi']
    h_zp['PDE_ext2']=Sigma_gas*sigma_z*np.sqrt(2*rho_sd)
    h_zp['PDE_ext2']*=(unit_surf*unit_v*np.sqrt(c.G*unit_rho)/c.k_B).to('K/cm**3').value
    h_zp['PDE_sp']=np.pi*Sigma_gas*Sigma_sp
    h_zp['PDE_sp']*=(c.G*unit_surf**2/c.k_B).to('K/cm**3').value
    h_zp['PDE_sg']=h_zp['PDE_gas']+h_zp['PDE_sp']
    h_zp['PDE']=h_zp['PDE_gas']+h_zp['PDE_ext']+h_zp['PDE_sp']
    h_zp['PDE2']=h_zp['PDE_gas']+h_zp['PDE_ext2']+h_zp['PDE_sp']

    h_zp.index=zpw.taxis
     
    h_zp.to_pickle(hstfile)
    print('{} is created'.format(hstfile))
    return h_zp
    
def doall(base,problem_id,problem_dir=None,do_yt=True):
    """
       do all preprocessing for a model
    """
    if problem_dir is None: problem_dir = problem_id

    done=cleanup_directory(base,problem_id,problem_dir=problem_dir)
    for k in done:
        if done[k]: print('{}.{} is moved'.format(problem_id,k)) 

    parfile='{}{}/{}.par'.format(base,problem_dir,problem_id)
    params=pa.get_params(parfile)    

    zpfile='{}{}/zprof_merged/{}.whole.zprof.nc'.format(base,problem_id,problem_id)
    if not os.path.isfile(zpfile): zprof_to_xarray(base,problem_id)
    #zpfile='{}{}/zprof_merged/{}.zprof.nc'.format(base,problem_id,problem_id)
    #if not os.path.isfile(zpfile): merge_xarray(base,problem_id)

    hstzpfile='{}{}/hst/{}.hst_zp.p'.format(base,problem_id,problem_id)
    if not os.path.isfile(hstzpfile): recal_history(base,problem_id)

    if do_yt:
        kwargs={'base_directory':base,
            'directory':'{}/'.format(problem_dir),
            'id':problem_id,
            'parallel':False,
            'phase':False,
            'rotation':params['Omega']*1.e3,
            'range':''
        }
        yt_analysis.main(**kwargs)
