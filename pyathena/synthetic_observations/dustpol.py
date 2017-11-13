from .los_to_dustpol import los_to_dustpol
import healpy as hp
import pandas as pd
import numpy as np
import os
    
def load_los(domain,srange=None,bmin=-1,ithread=0,nthread=1,Nside=4,center=[0.,0.,0.]):
    deltas=domain['dx'][2]/2.

    losdir=domain['losdir']
    step=domain['step']
    outdir='%s%s/Nside%d' % (losdir,step,Nside)
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
    
    outfile='%s/%d.p' % (outdir,0)
    if not os.path.isfile(outfile): print("There is no correpnding LOS file: %s" % outfile)
    los=pd.read_pickle(outfile)
    if srange != None: sidx=(los.index >= srange[0]) & (los.index <= srange[1])

    npix=hp.nside2npix(Nside)
    npix_per_thread=int(npix/nthread)
    npix_min=npix_per_thread*ithread
    npix_max=npix_per_thread*(ithread+1)
    
    los_all=[]
    pix_arr=[]
    for ipix in range(npix_min,npix_max):
        angle = np.rad2deg(hp.pix2ang(Nside,ipix))
        if np.abs(90-angle[0]) > bmin:
            outfile='%s/%d.p' % (outdir,ipix)
            los=pd.read_pickle(outfile)
            if srange != None: los=los[sidx]
            los_all.append(los)
            pix_arr.append(ipix)
    return los_all,pix_arr
                
def make_pol_map(los_all,pix_arr,domain,Imap,Umap,Qmap,srange=None,Trange=None):
    deltas=domain['dx'][2]/2.

    los=los_all[0]
    if srange != None: sidx=(los.index >= srange[0]) & (los.index <= srange[1])

    args={'Bnu':41495.876171482356, 'sigma':1.e-26, 'p0':0.2, 'attenuation': 0}

    for ipix,los in list(zip(pix_arr,los_all)):
        if srange != None: los=los[sidx]
        if Trange != None: 
            Tidx=(los['temperature'] >= Trange[0]) & (los['temperature'] <= Trange[1])
            los=los[Tidx]
        nH=los['density']
        Bx=los['magnetic_field_X']
        By=los['magnetic_field_Y']
        Bz=los['magnetic_field_Z']
        I,Q,U=los_to_dustpol(nH,Bx,By,Bz,deltas,args)
        Imap[ipix]=I
        Qmap[ipix]=Q
        Umap[ipix]=U

def make_map(domain,Nside=4,center=[0,0,0],srange=None,Trange=None):
    deltas=domain['dx'][2]/2.
    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
    los=[]
    for f in ['density','magnetic_fieldX','magnetic_fieldY','magnetic_fieldZ','temperature']:
        outfile='%s/%s.p' % (outdir,f)
        los.append(pd.read_pickle(outfile))
    nH,Bx,By,Bz,temp=los
    sarr=nH.columns
    if srange != None: 
        sidx=(sarr >= srange[0]) & (sarr <= srange[1])
        nH=nH.iloc[:,sidx]
        Bx=Bx.iloc[:,sidx]
        By=By.iloc[:,sidx]
        Bz=Bz.iloc[:,sidx]
        temp=temp.iloc[:,sidx]
    if Trange != None: 
        Tidx=(temp >= Trange[0]) & (temp <= Trange[1])
        nH=nH[Tidx]
        Bx=Bx[Tidx]
        By=By[Tidx]
        Bz=Bz[Tidx]

    args={'Bnu':41495.876171482356, 'sigma':1.e-26, 'p0':0.2, 'attenuation': 0}
    Bnu=args['Bnu']
    p0=args['p0']
    sigma=args['sigma']

    Bperp2=Bx*Bx+By*By
    B2=Bperp2+Bz*Bz
    cos2phi=(Bx*Bx-By*By)/Bperp2
    sin2phi=Bx*By/Bperp2
    cosgam2=Bperp2/B2

    ds=deltas*3.085677581467192e+18
    dtau=sigma*nH*ds

    I=Bnu*(1.0-p0*(cosgam2-2./3.0))*dtau
    Q=p0*Bnu*cos2phi*cosgam2*dtau
    U=p0*Bnu*sin2phi*cosgam2*dtau

    return I.sum(axis=1),Q.sum(axis=1),U.sum(axis=1)
