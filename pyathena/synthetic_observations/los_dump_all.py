import os
from .tools import get_hat,get_joffset
from .reader import read_data
import healpy as hp
import numpy as np
import pandas as pd

from ..utils import cc_idx

def los_idx_all(hat,domain,smin=0.,smax=3000.,ds=1.,center=[0,0,0],zmax_cut=True):
    zmax=domain['right_edge'][2]-0.5*domain['dx'][2]
    zmin=domain['left_edge'][2]+0.5*domain['dx'][2]
    xhat=hat[0][:,np.newaxis]
    yhat=hat[1][:,np.newaxis]
    zhat=hat[2][:,np.newaxis]

    sarr=np.arange(smin,smax,ds)
    xarr=xhat*sarr + center[0]
    yarr=yhat*sarr + center[1]
    zarr=zhat*sarr + center[2]

    iarr = cc_idx(domain,[xarr,yarr,zarr])

    return iarr,[xarr,yarr,zarr],sarr


def los_dump(ds,domain,Nside=4,center=[0,0,0],force_write=False):
    deltas=domain['dx'][2]/2.
    smax=domain['Lx'][2]/2

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s/Nside%d-%s' % (losdir,step,Nside,cstring)
 
    x0,y0,z0,dx,dy,dz,vy0 = get_index_all(domain,Nside,center,smax,deltas)

    for f in domain['fields']:
        outfile='%s/%s.npy' % (outdir,f)
        if not os.path.isfile(outfile) or force_write:
            data=read_data(ds,f,domain)
            newdata=extend_data(domain,data)
            d000=newdata[z0  ,y0  ,x0  ]*(1-dz)*(1-dy)*(1-dx)
            d100=newdata[z0+1,y0  ,x0  ]*dz*(1-dy)*(1-dx)
            d010=newdata[z0  ,y0+1,x0  ]*(1-dz)*dy*(1-dx)
            d001=newdata[z0  ,y0  ,x0+1]*(1-dz)*(1-dy)*dx
            d110=newdata[z0+1,y0+1,x0  ]*dz*dy*(1-dx)
            d101=newdata[z0+1,y0  ,x0+1]*dz*(1-dy)*dx
            d011=newdata[z0  ,y0+1,x0+1]*(1-dz)*dy*dx
            d111=newdata[z0+1,y0+1,x0+1]*dz*dy*dx
 
            dlos=d000+d100+d010+d001+d110+d101+d011+d111
            if f is 'velocity2' and 'Omega' in domain: dlos += vy0
            np.save(outfile,dlos)

def get_index_all(domain,Nside,center,smax,ds):

    npix=hp.nside2npix(Nside)
    ipix=np.arange(npix)
    hat=get_hat(Nside,ipix)
    iarr,xarr,sarr=los_idx_all(hat['Z'],domain,smin=0,smax=smax,ds=ds,center=center)

    Nx,Ny,Nz=domain['Nx']

    if 'Omega' in domain:
        joffset=get_joffset(domain)
        Omega=domain['Omega']
        qshear=domain['qshear']
        vy0 = -qshear*Omega*xarr[0]
    else:
        joffset=0
        vy0 = None

    xdiv,xidx=np.divmod(iarr[0],Nx)
    yidx=np.remainder(iarr[1]+xdiv*joffset,Ny)
    zidx=iarr[2]

    x0 = xidx.astype(np.intp)
    y0 = yidx.astype(np.intp)
    z0 = zidx.astype(np.intp)
    dx = xidx - x0
    dy = yidx - y0
    dz = zidx - z0

    return x0,y0,z0,dx,dy,dz,vy0

def extend_data(domain,data):
    joffset=get_joffset(domain)

    dslicey=data[:,0,:]
    newdata=np.concatenate((data,dslicey[:,np.newaxis,:]),axis=1)
    d1=np.roll(newdata[:,:,0],-joffset.astype(np.int),axis=1)
    d2=np.roll(newdata[:,:,0],-(joffset.astype(np.int)+1),axis=1)
    dj=joffset-joffset.astype(np.int)
    dslicex=d1*(1-dj)+d2*dj

    newdata=np.concatenate((newdata,dslicex[:,:,np.newaxis]),axis=2)
    return newdata

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
