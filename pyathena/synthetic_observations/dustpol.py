import os
from .get_los import get_los_cy as get_los
import healpy as hp
import pandas as pd
import numpy as np

def los_dump(data,domain,ithread=0,nthread=1,Nside=4,center=[0.,0.,0.]):
    deltas=domain['dx'][2]/2.
    smax=domain['Lx'][0]*2

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
    
    if not os.path.isdir(losdir): os.mkdir(losdir)
    if not os.path.isdir(losdir+step): os.mkdir(losdir+step)
    if not os.path.isdir(outdir): os.mkdir(outdir)

    npix=hp.nside2npix(Nside)
    npix_per_thread=int(npix/nthread)
    npix_min=npix_per_thread*ithread
    npix_max=npix_per_thread*(ithread+1)
    
    for ipix in range(npix_min,npix_max):
        outfile='%s/%d.p' % (outdir,ipix)
        if not os.path.isfile(outfile):
            los=get_los(data,domain,Nside,ipix,smax=smax,deltas=deltas,center=center)
            df=pd.DataFrame.from_dict(los)
            df.index=df.pop('sarr')
            pd.DataFrame(df).to_pickle(outfile)

def make_pol_map(domain,Imap,Umap,Qmap,srange=None,
                 bmin=-1,ithread=0,nthread=1,Nside=4,center=[0.,0.,0.]):
    deltas=domain['dx'][2]/2.

    losdir=domain['losdir']
    step=domain['step']
    outdir='%s%s/Nside%d' % (losdir,step,Nside)
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
    
    outfile='%s/%d.p' % (outdir,0)
    if os.path.isfile(outfile): print("There is no correpnding LOS file: %s" % outfile)
    los=pd.read_pickle(outfile)
    if srange != None: sidx=(los.index >= srange[0]) & (los.index <= srange[1])

    npix=hp.nside2npix(Nside)
    npix_per_thread=int(npix/nthread)
    npix_min=npix_per_thread*ithread
    npix_max=npix_per_thread*(ithread+1)
    
    for ipix in range(npix_min,npix_max):
        angle = np.rad2deg(hp.pix2ang(Nside,ipix))
        if np.abs(90-angle[0]) > bmin:
            outfile='%s/%d.p' % (outdir,ipix)
            los=pd.read_pickle(outfile)
            if srange != None: los=los[sidx]
            nH=los['density']
            Bfield=[los['magnetic_field_X'],los['magnetic_field_Y'],los['magnetic_field_Z']]
            I,Q,U=los_to_dustpol(nH,Bfield,deltas=deltas)
            Imap[ipix]=I
            Qmap[ipix]=Q
            Umap[ipix]=U
        else:
            Imap[ipix]=hp.UNSEEN
            Qmap[ipix]=hp.UNSEEN
            Umap[ipix]=hp.UNSEEN

def los_to_dustpol(nH,Bfield,\
                   Bnu=41495.876171482356,sigma=1.2e-26,deltas=1.,p0=0.2,attenuation=False):
    """
    inputs: 
        nH in units of cm^-3
        Bfield (X,Y,Z) following the Healpix convention
    parameters:
        Bnu(T_dust,nu0): planck function (default T_dust=18K, nu0=353GHz; Planck XX 2015)
        sigma: absorption crossection in units of cm^-2 at 353GHz (default 1.2e-26 cm^-2; Planck XX 2015)
        deltas: length of lins segments in units of pc
        p0: intrinsic polarization fraction (default 20%; Planck XX 2015)
        attenuation: if true, self-attenuation is considered (default False)
    output:
        I, Q, U stokes' parameters in units of the input Bnu units. (default is MJy/sr)
    """

    ds=deltas*3.085677581467192e+18

    Bx=Bfield[0]  # tengential component
    By=Bfield[1]  # tengential component
    Bz=Bfield[2]  # LOS component

    Bperp2=Bx**2+By**2
    B2=Bperp2+Bz**2
    cos2phi=(Bx**2-By**2)/Bperp2
    sin2phi=Bx*By/Bperp2
    cosgam2=Bperp2/B2

    dtau=sigma*nH*ds
 #   print nH.sum()*ds.cgs,tau[-1], np.exp(-tau[-1])

    I=Bnu*(1.0-p0*(cosgam2-2./3.0))*dtau#*np.exp(-tau)
    Q=p0*Bnu*cos2phi*cosgam2*dtau#*np.exp(-tau)
    U=p0*Bnu*sin2phi*cosgam2*dtau#*np.exp(-tau)
    if attenuation:
        tau=dtau.cumsum()
        I=I*np.exp(-tau)
        Q=Q*np.exp(-tau)
        U=U*np.exp(-tau)
    
    return I.sum(),Q.sum(),U.sum()

