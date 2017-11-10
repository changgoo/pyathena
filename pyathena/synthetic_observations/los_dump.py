import os
from .get_los import get_los_cy as get_los
import healpy as hp
import pandas as pd

def make_directory(domain,Nside=4,center=[0.,0.,0.]):
    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
    
    if not os.path.isdir(losdir): os.mkdir(losdir)
    if not os.path.isdir(losdir+step): os.mkdir(losdir+step)
    if not os.path.isdir(outdir): os.mkdir(outdir)


def los_dump(data,domain,ithread=0,nthread=1,Nside=4,center=[0.,0.,0.],force_write=False):
    deltas=domain['dx'][2]/2.
    smax=domain['Lx'][0]*2

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
            los=get_los(data,domain,Nside,ipix,smax=smax,deltas=deltas,center=center)
            df=pd.DataFrame.from_dict(los)
            df.index=df.pop('sarr')
            pd.DataFrame(df).to_pickle(outfile)
    if nthread>1: 
        etime=time.time()
        print("Exiting %d: Wall time %g" % (ithread,etime-stime))


