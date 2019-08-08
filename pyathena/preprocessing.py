
import xarray as xr
import pandas as pd
import os,glob
import shutil


from pyathena.parse_par import write_par_from_rst
import pyathena as pa

import astropy.constants as c
import astropy.units as u
import numpy as np



def cleanup_directory(base,problem_id,problem_dir=None,newbase=None):
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
    fields_to_move=['zprof','hst','rst','starpar']
    if newbase is None: 
        newbase = base
    else:
        fields_to_move += ['vtk']

    for ext in fields_to_move:
        subdir='{}{}/{}'.format(newbase,problem_dir,ext)
        success[ext]=True
        if (not os.path.isdir(subdir)) and (not ext is 'vtk'): os.mkdir(subdir)
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
        elif ext is 'star':
            files=glob.glob('{}{}/id0/{}.?????.star'.format(base,problem_dir,problem_id))
        elif ext is 'vtk':
            print('moving {} files...'.format(ext))
            files=[]
        else:
            print('{} extention specifier is not known'.format(ext))
            success[ext]=False

        if len(files):
            print('moving {} files...'.format(ext))
            for f in files:
                f2=f.replace('id0/','{}/'.format(ext))
                f2=f2.replace(base,newbase)
                if ext is 'hst':
                    shutil.copy(f,f2)
                    if base != base: 
                        f2=f.replace(base,newbase)
                        shutil.copy(f,f2)
                else:
                    shutil.move(f,f2)
        else:
            success[ext]=False

        if ext is 'rst':
            for i in range(1,nprocs):
                files=glob.glob('{}{}/id{}/{}-id{}.????.rst'.format(base,problem_dir,i,problem_id,i))
                if len(files): 
                    for f in files: 
                        f2=f.replace('id{}/'.format(i),'rst/').replace(base,newbase)
                        shutil.move(f,f2)

        if ext is 'vtk':
            files=glob.glob('{}{}/id0/{}.????.vtk'.format(base,problem_dir,problem_id))
            if len(files): 
                for f in files: 
                    f2=f.replace(base,newbase)
                    shutil.move(f,f2)
                success[ext]=True
            else:
                success[ext]=False

            for i in range(1,nprocs):
                files=glob.glob('{}{}/id{}/{}-id{}.????.vtk'.format(base,problem_dir,i,problem_id,i))
                if len(files): 
                    for f in files: 
                        f2=f.replace(base,newbase)
                        shutil.move(f,f2)

    return success

def zprof_to_xarray(base,problem_dir,problem_id,concatenated=True):
    """
        merge zprof dumps along t-axis and create netCDF files 
        files should have moved to zprof/
    """    
    zpmerge_dir='{}{}/zprof_merged/'.format(base,problem_dir)
    if not os.path.isdir(zpmerge_dir): os.mkdir(zpmerge_dir)
    plist=['phase1','phase2','phase3','phase4','phase5']
    for phase in plist:
        zprof_fnames=glob.glob('{}{}/zprof/{}.*.{}.zprof'.format(base,problem_dir,problem_id,phase))
        zprof_fnames.sort()
        zpfile='{}{}/zprof_merged/{}.{}.zprof.nc'.format(base,problem_dir,problem_id,phase)
        if os.path.isfile(zpfile) and concatenated: 
            with xr.open_dataarray(zpfile) as da: da.load()
            nfiles_in_zpmerged=len(da.taxis)
            nappend=len(zprof_fnames)-nfiles_in_zpmerged
            if nappend>0: 
                da_part=read_zprof(zprof_fnames[nfiles_in_zpmerged:])
                da=xr.concat((da,da_part),dim='taxis')
                print('{} files have been merged. {} will be appended'.format(
                  nfiles_in_zpmerged,nappend))
        else:
            da=read_zprof(zprof_fnames)
        da.to_netcdf(zpfile)
        print('{} is created'.format(zpfile))

def read_zprof_one(zprof_fname):
    with open(zprof_fname,'r') as fp:
      hd=fp.readline()
      time=float(hd[hd.rfind('t=')+2:])
    df=pd.read_csv(zprof_fname,skiprows=1)
    return df,time

def read_zprof(zprof_fnames):
    taxis=[]
    dfall=None
    for f in zprof_fnames:
        df,time=read_zprof_one(f)
        zaxis=np.array(df['z'])
        fields=np.array(df.columns.get_values())
        taxis.append(time)
        if dfall is None:
            dfall=np.array(df)[np.newaxis,:]
        else:
            dfall=np.concatenate([dfall,np.array(df)[np.newaxis,:]],axis=0)
    da=xr.DataArray(dfall.T,coords={'fields':fields,'zaxis':zaxis,'taxis':taxis},
                        dims=('fields','zaxis','taxis'))
    return da

def merge_xarray(base,problem_dir,problem_id):
    """
        merge all phases into a single file
    """    
    plist=['phase1','phase2','phase3','phase4','phase5']
    datasets = xr.Dataset()
    for phase in plist:
        path='{}{}/zprof_merged/{}.{}.zprof.nc'.format(base,problem_dir,problem_id,phase)
        with xr.open_dataarray(path) as da: da.load()
        datasets[phase]=da
    return datasets
    #zpfile='{}{}/zprof_merged/{}.zprof.nc'.format(base,problem_id,problem_id)
    #datasets.to_netcdf(zpfile)
    #print('{} is created'.format(zpfile))

def panel_to_xarray(base,problem_dir,problem_id):
    """
        convert exsiting merged pickles for pandas Panel to netCDF for xarray DataArray
    """
    zpfile='{}{}/zprof_merged/{}.whole.zprof.p'.format(base,problem_dir,problem_id)
    if not os.path.isfile(zpfile):
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

def recal_rates(h,sn,base,problem_dir,problem_id):
    import glob
    from pyathena.set_plt import units
    from pyathena.parse_par import parse_par
    import pandas as pd
    
    par,blocks,fields,=parse_par('{}{}/{}.par'.format(base,problem_dir,problem_id))
    dt=float(par['output1']['dt'][0])
    Lx=float(par['domain1']['x1max'][0])-float(par['domain1']['x1min'][0])
    Ly=float(par['domain1']['x2max'][0])-float(par['domain1']['x2min'][0])
    area=Lx*Ly
    
    nhst=len(h.index)
    
    starvtk=glob.glob('{}{}/starpar/{}.*.starpar.vtk'.format(base,problem_dir,problem_id))
    starvtk.sort()
    sp=pa.read_starvtk(starvtk[-1])
    cl=sp[sp.mass!=0]
    cl_birth_time=(cl.time-cl.age)#*units['Myr']
    cl_mass=cl.mass*units['Msun']
    
    time=np.arange(nhst+1)*dt+h.index[0]*units['Myr']
    Msp,t=np.histogram(cl_birth_time,bins=time,weights=cl_mass)
    sfr_inst=pd.Series(Msp/area/dt)

    rates=pd.DataFrame()
    rates['sfr']=sfr_inst

    typeIa=sn[(sn.runaway==2) | (sn.runaway==-1)]
    runaway=sn[sn.runaway==1]
    cluster=sn[sn.runaway==0]
    for lab,sndata in zip(['','_run','_cl','_Ia'],[sn,runaway,cluster,typeIa]):
        sn_time=sndata.time*units['Myr']
        NSN,t=np.histogram(sn_time,bins=time)
        rates['snr{}'.format(lab)]=NSN/dt/area

    for rinst_key in ['sfr','snr','snr_run','snr_cl','snr_Ia']:
        rinst=rates[rinst_key] 
        for tbin in [1.,10.,40.,100.]:
            window=int(tbin/dt)
            rates['{}{}'.format(rinst_key,int(tbin))]=rinst.rolling(window,min_periods=1).mean()

    rates.index=time[:-1]
    return rates

def recal_history(base,problem_dir,problem_id):
    parfile='{}{}/{}.par'.format(base,problem_dir,problem_id)
    params=pa.get_params(parfile)

    hstfile='{}{}/hst/{}.hst'.format(base,problem_dir,problem_id)
    hst=pa.hst_reader(hstfile)

    hstrecalfile='{}_cal.p'.format(hstfile)
    h=processing_history_dump(hst,params,hstrecalfile)

    snfile='{}{}/hst/{}.sn'.format(base,problem_dir,problem_id)
    sn=pa.hst_reader(snfile)
    rates=recal_rates(h,sn,base,problem_dir,problem_id)

    hstzpfile='{}_zp.p'.format(hstfile)
    zprof_ds=merge_xarray(base,problem_dir,problem_id)
    #print("zprof dataset is loaded")
    h_zp=processing_zprof_dump(h,rates,params,zprof_ds,hstzpfile)

def processing_history_dump(hst,params,hstfile):
    from pyathena.set_plt import unit,units
    toMsun=units['Msun']
    toMyr=units['Myr']
    toflux=units['massflux']
    toPok=units['P']

    Lz=params['x3max']-params['x3min']
    Ly=params['x2max']-params['x2min']
    Lx=params['x1max']-params['x1min']
    Nx=params['Nx1']
    Ny=params['Nx2']
    Nz=params['Nx3']
    Ntot=Nx*Ny*Nz
    vol=Lx*Ly*Lz
    area=Lx*Ly
    dz=Lz/Nz
    Omega=params['Omega']
    torb=2*np.pi/Omega*toMyr

    mhd=False
    if 'x1ME' in hst: mhd=True
    #if mhd: print("this is MHD output")

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
    h['heat_ratio']=hst['heat_ratio']
    if 'ftau' in hst: h['ftau']=hst['ftau']

    h.index=hst['time']
    h.to_pickle(hstfile)
    print('{} is created'.format(hstfile))

    return h

def processing_zprof_dump(h,rates,params,zprof_ds,hstfile):
    from pyathena.set_plt import unit,units
    toMsun=units['Msun']
    toMyr=units['Myr']
    toflux=units['massflux']
    toPok=units['P']

    Lz=params['x3max']-params['x3min']
    Ly=params['x2max']-params['x2min']
    Lx=params['x1max']-params['x1min']
    Nx=params['Nx1']
    Ny=params['Nx2']
    Nz=params['Nx3']
    Ntot=Nx*Ny*Nz
    vol=Lx*Ly*Lz
    area=Lx*Ly
    dz=Lz/Nz
    Omega=params['Omega']
    torb=2*np.pi/Omega*toMyr

    mhd=False
    if 'B1' in zprof_ds.fields.to_index(): mhd=True
    nscal=0
    for f in zprof_ds.fields.data:
        if f.startswith('s'): nscal +=1
    #if mhd: print("this is MHD output")

    # process zprof dump
    h_zp=pd.DataFrame()
    zpw=zprof_ds.to_array().sum(dim='variable')
    zp2p=zprof_ds.to_array().sel(variable=['phase1','phase2','phase3']).sum(dim='variable')
    zph=zprof_ds.to_array().sel(variable=['phase4','phase5']).sum(dim='variable')

    h_zp['tMyr']=zpw.taxis*toMyr
    h_zp['torb']=h_zp['tMyr']/torb

    for i,ph in enumerate(['_2p','_h','','_c','_u','_w','_h1','_h2',]):
        dtot=zpw.sel(fields='d').sum(dim='zaxis')
        Atot=zpw.sel(fields='A').sum(dim='zaxis')
        if ph is '_2p':
            zp=zp2p
        elif ph is '_h':
            zp=zph
        elif ph is '':
            zp=zpw
        else:
            zp=zprof_ds['phase{}'.format(i-2)]
        dphase=zp.sel(fields='d').sum(dim='zaxis')
        h_zp['surf{}'.format(ph)]=(zp.sel(fields='d')/zpw.sel(fields='A')).sum(dim='zaxis')*toMsun*dz
        if ph is not '':
            h_zp['mf{}'.format(ph)]=zp.sel(fields='d').sum(dim='zaxis')/dtot
            h_zp['vf{}'.format(ph)]=zp.sel(fields='A').sum(dim='zaxis')/Atot
        h_zp['H{}'.format(ph)]=np.sqrt((zp.sel(fields='d')*zp.zaxis**2).sum(dim='zaxis')/dphase)
        for ax in ['1','2','3']:
            Ekf = 'Ek{}'.format(ax)
            if (ax == '2') and ('dEk2' in np.array(zp.fields)): Ekf = 'dEk2'
            h_zp['v{}{}'.format(ax,ph)] = np.sqrt(2.0*zp.sel(fields=Ekf).sum(dim='zaxis')/dphase)
            h_zp['Ek{}{}'.format(ax,ph)] = zp.sel(fields=Ekf).sum(dim='zaxis')/Atot
            if mhd:
                h_zp['vA{}{}'.format(ax,ph)] = np.sqrt(2.0*zp.sel(fields='PB{}'.format(ax)).sum(dim='zaxis')/dphase)
                h_zp['dvA{}{}'.format(ax,ph)] = np.sqrt(2.0*zp.sel(fields='dPB{}'.format(ax)).sum(dim='zaxis')/dphase)
                h_zp['EB{}{}'.format(ax,ph)] = zp.sel(fields='PB{}'.format(ax)).sum(dim='zaxis')/Atot
                h_zp['EB{}_turb{}'.format(ax,ph)] = zp.sel(fields='dPB{}'.format(ax)).sum(dim='zaxis')/Atot

        h_zp['v{}'.format(ph)] = np.sqrt(h_zp['v1{}'.format(ph)]**2 + \
                                         h_zp['v2{}'.format(ph)]**2 + \
                                         h_zp['v3{}'.format(ph)]**2)
        h_zp['Ek{}'.format(ph)] = h_zp['Ek1{}'.format(ph)] + \
                                  h_zp['Ek2{}'.format(ph)] + \
                                  h_zp['Ek3{}'.format(ph)]
        if mhd:
            h_zp['vA{}'.format(ph)] = np.sqrt(h_zp['vA1{}'.format(ph)]**2 + \
                                              h_zp['vA2{}'.format(ph)]**2 + \
                                              h_zp['vA3{}'.format(ph)]**2)
            h_zp['dvA{}'.format(ph)] = np.sqrt(h_zp['dvA1{}'.format(ph)]**2 + \
                                               h_zp['dvA2{}'.format(ph)]**2 + \
                                               h_zp['dvA3{}'.format(ph)]**2)
            h_zp['EB{}'.format(ph)] = h_zp['EB1{}'.format(ph)] + \
                                      h_zp['EB2{}'.format(ph)] + \
                                      h_zp['EB3{}'.format(ph)]
            h_zp['EB_turb{}'.format(ph)] = h_zp['EB1_turb{}'.format(ph)] + \
                                           h_zp['EB2_turb{}'.format(ph)] + \
                                           h_zp['EB3_turb{}'.format(ph)]


        h_zp['cs{}'.format(ph)] = np.sqrt(zp.sel(fields='P').sum(dim='zaxis')/dphase)
        h_zp['Eth{}'.format(ph)] = 1.5*zp.sel(fields='P').sum(dim='zaxis')/Atot
        h_zp['sigma_eff{}'.format(ph)] = h_zp['cs{}'.format(ph)]**2 + h_zp['v3{}'.format(ph)]**2
        if mhd:
            h_zp['sigma_eff{}'.format(ph)] = h_zp['sigma_eff{}'.format(ph)] + \
                                             0.5*(h_zp['vA1{}'.format(ph)]**2 + 
                                                  h_zp['vA2{}'.format(ph)]**2 - 
                                                  h_zp['vA3{}'.format(ph)]**2)
        h_zp['sigma_eff{}'.format(ph)] = np.sqrt(h_zp['sigma_eff{}'.format(ph)])

        zp_upper=zp[:,-1,:]
        zp_lower=zp[:,0,:]
        h_zp['v3_ubd{}'.format(ph)]=(zp_upper.sel(fields='pFzd')/zp_upper.loc['pd'])
        h_zp['v3_lbd{}'.format(ph)]=(-zp_lower.sel(fields='mFzd')/zp_lower.loc['md'])
        h_zp['v3_bd{}'.format(ph)]=0.5*(h_zp['v3_ubd{}'.format(ph)]+h_zp['v3_lbd{}'.format(ph)])

        for ns in range(nscal+1):
            if ns is 0: var='d'
            else: var='s{}'.format(ns)
            field_name_u='massflux_{}_{}{}'.format('ubd',var,ph)
            flx=zp_upper.loc['pFz{}'.format(var)]/area
            h_zp[field_name_u]=flx*toflux
            field_name_l='massflux_{}_{}{}'.format('lbd',var,ph)
            flx=-zp_lower.loc['mFz{}'.format(var)]/area
            h_zp[field_name_l]=flx*toflux
            field_name='massflux_{}_{}{}'.format('bd',var,ph)
            h_zp[field_name]=h_zp[field_name_u]+h_zp[field_name_l]

        for var in ['M1','M2','M3','E1','E2','E3','P']:
            field_name_u='flux_{}_{}{}'.format('ubd',var,ph)
            flx=zp_upper.loc['pFz{}'.format(var)]/area
            h_zp[field_name_u]=flx
            field_name_l='flux_{}_{}{}'.format('lbd',var,ph)
            flx=-zp_lower.loc['mFz{}'.format(var)]/area
            h_zp[field_name_l]=flx
            field_name='flux_{}_{}{}'.format('bd',var,ph)
            h_zp[field_name]=h_zp[field_name_u]+h_zp[field_name_l]

        for z0 in [500,1000]:
            zstr='{:02d}'.format(int(z0/100))
            zp_upper=zp.sel(zaxis=z0,method='nearest')
            zp_lower=zp.sel(zaxis=-z0,method='nearest')
            for ns in range(nscal+1):
                if ns is 0: var='d'
                else: var='s{}'.format(ns)
                field_name_u='massflux_out_u{}_{}{}'.format(zstr,var,ph)
                flx=zp_upper.loc['pFz{}'.format(var)]/area
                h_zp[field_name_u]=flx*toflux
                field_name_l='massflux_out_l{}_{}{}'.format(zstr,var,ph)
                flx=-zp_lower.loc['mFz{}'.format(var)]/area
                h_zp[field_name_l]=flx*toflux
                field_name='massflux_out_{}_{}{}'.format(zstr,var,ph)
                h_zp[field_name]=h_zp[field_name_u]+h_zp[field_name_l]
 
                field_name_u='massflux_in_u{}_{}{}'.format(zstr,var,ph)
                flx=zp_upper.loc['mFz{}'.format(var)]/area
                h_zp[field_name_u]=flx*toflux
                field_name_l='massflux_in_l{}_{}{}'.format(zstr,var,ph)
                flx=-zp_lower.loc['pFz{}'.format(var)]/area
                h_zp[field_name_l]=flx*toflux
                field_name='massflux_in_{}_{}{}'.format(zstr,var,ph)
                h_zp[field_name]=h_zp[field_name_u]+h_zp[field_name_l]
 
                field_name='massflux_{}_{}{}'.format(zstr,var,ph)
                field_name_out='massflux_out_{}_{}{}'.format(zstr,var,ph)
                field_name_in='massflux_in_{}_{}{}'.format(zstr,var,ph)
                h_zp[field_name]=h_zp[field_name_out]+h_zp[field_name_in]

                h_zp['v3_out_u{}{}'.format(zstr,ph)]=\
                  (zp_upper.sel(fields='pFzd')/zp_upper.loc['pd'])
                h_zp['v3_out_l{}{}'.format(zstr,ph)]=\
                  (-zp_lower.sel(fields='mFzd')/zp_lower.loc['md'])
                h_zp['v3_out_{}{}'.format(zstr,ph)]=\
                  0.5*(h_zp['v3_out_u{}{}'.format(zstr,ph)]+ \
                       h_zp['v3_out_l{}{}'.format(zstr,ph)])

        mid_idx=np.abs(zp.zaxis) < dz
        zpmid=zp[:,mid_idx,:]
        h_zp['vf_mid{}'.format(ph)] = zpmid.loc['A'].sum(axis=0)/zpw[:,mid_idx,:].loc['A'].sum(axis=0)
        h_zp['nmid{}'.format(ph)] = zpmid.loc['d'].sum(axis=0)/zpmid.loc['A'].sum(axis=0)
        h_zp['Pth_mid{}'.format(ph)] = (zpmid.loc['P']/zpmid.loc['A']).mean(axis=0)*toPok
        h_zp['Pturb_mid{}'.format(ph)] = (2.0*zpmid.loc['Ek3']/zpmid.loc['A']).mean(axis=0)*toPok
        h_zp['Pmid{}'.format(ph)] = h_zp['Pth_mid{}'.format(ph)]+h_zp['Pturb_mid{}'.format(ph)]
        if mhd:
            oPmag = 0.5*(zpmid.loc['B1']/zpmid.loc['A'])**2
            oPmag = oPmag + 0.5*(zpmid.loc['B2']/zpmid.loc['A'])**2
            oPmag = oPmag + 0.5*(zpmid.loc['B3']/zpmid.loc['A'])**2
            h_zp['Pmag_mean_mid{}'.format(ph)] = oPmag.mean(axis=0)*toPok
            oPmag = oPmag - (zpmid.loc['B3']/zpmid.loc['A'])**2
            h_zp['Pimag_mean_mid{}'.format(ph)] = oPmag.mean(axis=0)*toPok
            h_zp['Pmid{}'.format(ph)] = h_zp['Pmid{}'.format(ph)] + \
                                        h_zp['Pimag_mean_mid{}'.format(ph)]

            tPmag = (zpmid.loc['dPB1']/zpmid.loc['A'])
            tPmag = tPmag + (zpmid.loc['dPB2']/zpmid.loc['A'])
            tPmag = tPmag + (zpmid.loc['dPB3']/zpmid.loc['A'])
            h_zp['Pmag_turb_mid{}'.format(ph)] = tPmag.mean(axis=0)*toPok
            tPmag = tPmag - 2.0*(zpmid.loc['dPB3']/zpmid.loc['A'])
            h_zp['Pimag_turb_mid{}'.format(ph)] = tPmag.mean(axis=0)*toPok
            h_zp['Pmid{}'.format(ph)] = h_zp['Pmid{}'.format(ph)] + \
                                        h_zp['Pimag_turb_mid{}'.format(ph)]


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
    for field in h:
        h_zp['{}_hst'.format(field)] = np.interp(zpw.taxis,h.index,h[field])

    for field in rates.keys():
        h_zp[field] = np.interp(zpw.taxis,rates.index,rates[field])

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
    Sigma_sp=h_zp['surfsp_hst']
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
    print('new {} is created'.format(hstfile))
    return h_zp
    
def draw_history(h_zp,metadata,figfname=''):
    from pyathena.set_plt import labels,plt
    nplot=len(metadata)
    fig=plt.figure(figsize=(10,4*nplot))
    iplot=1
    tmax=h_zp['tMyr'].max()
    for f,ylabel,yscale,subfields,ybottom,ytop in metadata:
        plt.subplot(nplot,1,iplot)
        if f in h_zp:
            l,=plt.plot(h_zp['tMyr'],h_zp[f])
            if len(subfields)>0:
                if f in labels: l.set_label(labels[f])
                for sf in subfields:
                    if sf in h_zp:
                        l,=plt.plot(h_zp['tMyr'],h_zp[sf])
                        if sf in labels: l.set_label(labels[sf])
                plt.legend(bbox_to_anchor=(1.02, 1),loc=2,fontsize='small')
        plt.ylabel(ylabel)
        plt.yscale(yscale)
        plt.ylim(bottom=ybottom,top=ytop)
        iplot+=1


    plt.xlabel('t [Myr]')
    if len(figfname) > 0:
        fig.savefig(figfname,bbox_inches='tight',dpi=150)

def doall(base,problem_id,problem_dir=None,do_pickling=True,use_yt=True,
  force_recal=False, force_redraw=False,vtkdir=None):
    """
        This function will do following tasks: 
        (1) reoranizing files 
           * .hst, .sn --> hst/
           * .zprof --> zprof/
           * .starpar.vtk --> starpar/
           * .rst --> rst/
        (2) creating a parameter file from a restart header
        (3) creating merged (along the time axis) zprof files for each thermal phase (using xarray)
        (4) creating history data using zprof dumps, including new fields
        (5) [if do_pickling] creating pickle files for slice and projection (will use yt [if use_yt])
    """
    if problem_dir is None: problem_dir = problem_id
    print('preparing metadata for {}...'.format(problem_id))
    done=cleanup_directory(base,problem_id,problem_dir=problem_dir)
    for k in done:
        if done[k]: print('{}.{} is moved'.format(problem_id,k)) 

    parfile='{}{}/{}.par'.format(base,problem_dir,problem_id)
    params=pa.get_params(parfile)    

    zpfile='{}{}/zprof_merged/{}.phase5.zprof.nc'.format(base,problem_dir,problem_id)
    from pyathena.utils import compare_files
    zpfiles=glob.glob('{}{}/zprof/{}.*.whole.zprof'.format(base,problem_dir,problem_id))
    zpfiles.sort()
    if len(zpfiles) > 0:
        if not compare_files(zpfiles[-1],zpfile): zprof_to_xarray(base,problem_dir,problem_id)
    #zpfile='{}{}/zprof_merged/{}.zprof.nc'.format(base,problem_id,problem_id)
    #if not os.path.isfile(zpfile): merge_xarray(base,problem_id)

    hstzpfile='{}{}/hst/{}.hst_zp.p'.format(base,problem_dir,problem_id)
    if os.path.isfile(zpfile): recal_history(base,problem_dir,problem_id)

    if do_pickling:
        kwargs={'base_directory':base,
            'directory':'{}/'.format(problem_dir),
            'id':problem_id,
            'parallel':False,
            'phase':False,
            'rotation':params['Omega']*1.e3,
            'range':''
        }

        if use_yt:
            print('slicing and projecting with yt ...')
            from pyathena.yt_analysis import yt_analysis
            yt_analysis.main(force_recal=force_recal,force_redraw=force_redraw,verbose=50,**kwargs)
        else:
            if not (vtkdir is None): kwargs['vtk_directory']=vtkdir
            print('slicing and projecting with pyathena ...')
            from pyathena.create_pickle import create_all_pickles
            create_all_pickles(force_recal=force_recal,force_redraw=force_redraw,verbose=True,**kwargs)
