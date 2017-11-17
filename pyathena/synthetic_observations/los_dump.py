import os
from .get_los import *
from .tools import get_hat
from .reader import read_data
import healpy as hp
import pandas as pd

def los_dump(data,domain,ithread=0,nthread=1,Nside=4,center=[0.,0.,0.],force_write=False):
    deltas=domain['dx'][2]/2.
    smax=domain['Lx'][2]/2

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
 
    npix=hp.nside2npix(Nside)
    npix_per_thread=int(npix/nthread)
    npix_min=npix_per_thread*ithread
    npix_max=npix_per_thread*(ithread+1)
    
    if nthread>1: 
        import time
        print("Starting %d" % (ithread))
        stime=time.time()
    for ipix in range(npix_min,npix_max):
        outfile='%s/%d.p' % (outdir,ipix)
        if not os.path.isfile(outfile) or force_write:
            los=get_los_all(data,domain,Nside,ipix,smax=smax,deltas=deltas,center=center)
            df=pd.DataFrame.from_dict(los)
            df.index=df.pop('sarr')
            pd.DataFrame(df).to_pickle(outfile)
    if nthread>1: 
        etime=time.time()
        print("Exiting %d: Wall time %g" % (ithread,etime-stime))

def los_dump_one(data,field,domain,ithread=0,nthread=1,Nside=4,center=[0.,0.,0.],force_write=False):
    deltas=domain['dx'][2]/2.
    smax=domain['Lx'][2]/2

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
 
    npix=hp.nside2npix(Nside)
    npix_per_thread=int(npix/nthread)
    npix_min=npix_per_thread*ithread
    npix_max=npix_per_thread*(ithread+1)
    
    if nthread>1: 
        import time
        print("Starting %d" % (ithread))
        stime=time.time()
        outfile='%s/%s-%d-%d.p' % (outdir,field,ithread,nthread)
    else:
        outfile='%s/%s.p' % (outdir,field)
    los=[]
    if not os.path.isfile(outfile) or force_write:
        for ipix in range(npix_min,npix_max):
            sarr,los_data=get_los_one(data,domain,Nside,ipix,smax=smax,deltas=deltas,center=center)
            if field is 'velocity2' and 'Omega' in domain:
                hat=get_hat(Nside,ipix)
                Omega=domain['Omega']
                qshear=domain['qshear'] 
                xarr=sarr*hat['Z'][0]+center[0]
                los_data -= qshear*Omega*xarr
            los.append(pd.Series(los_data,index=sarr))
        df=pd.DataFrame(los,index=range(npix_min,npix_max))
        pd.DataFrame(df).to_pickle(outfile)

    if nthread>1: 
        etime=time.time()
        print("Exiting %d: Wall time %g" % (ithread,etime-stime))

def merge_los_dump(domain,nthread,Nside=4,center=[0.,0.,0.],force_write=False):

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
 
    npix=hp.nside2npix(Nside)

    los={}
    
    for f in domain['fields']:
        outfile='%s/%s.p' % (outdir,f)
        for ithread in range(nthread):
            outfile_part='%s/%s-%d-%d.p' % (outdir,f,ithread,nthread)
            if ithread==0: 
                df=pd.read_pickle(outfile_part)
            else:
                df=df.append(pd.read_pickle(outfile_part))
        pd.DataFrame(df).to_pickle(outfile)

def los_dump_proj(domain,vec_field,Nside=4,center=[0.,0.,0.],force_write=False,ext='.npy'):

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
 
    npix=hp.nside2npix(Nside)

    los={}
    
    for f in ['1','2','3']:
        if ext=='.p': 
            outfile='%s/%s.p' % (outdir,vec_field+f)
            los[vec_field+f]=pd.read_pickle(outfile)
        elif ext=='.npy': 
            outfile='%s/%s.npy' % (outdir,vec_field+f)
            los[vec_field+f]=np.load(outfile)

    ipix = np.arange(npix)
    hat=get_hat(Nside,ipix)
    for axis in ['Z','X','Y']:
        los_out =hat[axis][0]*los[vec_field+'1'].T
        los_out+=hat[axis][1]*los[vec_field+'2'].T
        los_out+=hat[axis][2]*los[vec_field+'3'].T
        outfile='%s/%s%s' % (outdir,vec_field+axis,ext)

        if not os.path.isfile(outfile) or force_write:
            if ext=='.p':
                pd.DataFrame(los_out.T).to_pickle(outfile)
            elif ext=='.npy':
                np.save(outfile,los_out.T)
