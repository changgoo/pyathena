import os
from .tools import get_hat,get_joffset
from .reader import read_data
import healpy as hp
import numpy as np
import pandas as pd

def cc_idx(le,dx,pos):
    if np.array(pos).ndim == 2:
        le=le[:,np.newaxis]
        dx=dx[:,np.newaxis]
    elif np.array(pos).ndim == 3:
        le=le[:,np.newaxis,np.newaxis]
        dx=dx[:,np.newaxis,np.newaxis]

    idx=(pos-le-0.5*dx)/dx
    return idx


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

    le=domain['left_edge']
    dx=domain['dx']
    if np.abs(le[2]) < smax: le[2]=-smax
    iarr = cc_idx(le,dx,[xarr,yarr,zarr])

    return iarr,[xarr,yarr,zarr],sarr


def los_dump(ds,domain,deltas,smax,fields=['density'],
             Nside=4,center=[0,0,0],force_write=False):

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s-%d/Nside%d-%s' % (losdir,step,smax,Nside,cstring)
 
    x0,y0,z0,dx,dy,dz,vy0 = get_index_all(domain,Nside,center,smax,deltas)

    for f in fields:
        outfile='%s/%s.npy' % (outdir,f)
        if not os.path.isfile(outfile) or force_write:
            print('interpolating and writing: %s' % f)
            data=read_data(ds,f,domain)
            data=extend_data(domain,data,smax)
            dlos=data[z0  ,y0  ,x0  ]*(1-dz)*(1-dy)*(1-dx) +\
                 data[z0+1,y0  ,x0  ]*dz*(1-dy)*(1-dx) +\
                 data[z0  ,y0+1,x0  ]*(1-dz)*dy*(1-dx) +\
                 data[z0  ,y0  ,x0+1]*(1-dz)*(1-dy)*dx +\
                 data[z0+1,y0+1,x0  ]*dz*dy*(1-dx) +\
                 data[z0+1,y0  ,x0+1]*dz*(1-dy)*dx +\
                 data[z0  ,y0+1,x0+1]*(1-dz)*dy*dx +\
                 data[z0+1,y0+1,x0+1]*dz*dy*dx
 
            #dlos=d000+d100+d010+d001+d110+d101+d011+d111
            if f is 'velocity2' and 'Omega' in domain: dlos += vy0
            np.save(outfile,dlos)

def los_dump_from_data(data,domain,deltas,smax,f,
             Nside=4,center=[0,0,0],force_write=False):

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s-%d/Nside%d-%s' % (losdir,step,smax,Nside,cstring)
 
    x0,y0,z0,dx,dy,dz,vy0 = get_index_all(domain,Nside,center,smax,deltas)

    outfile='%s/%s.npy' % (outdir,f)
    if not os.path.isfile(outfile) or force_write:
        print('interpolating and writing: %s' % f)
        data=extend_data(domain,data,smax)
        dlos=data[z0  ,y0  ,x0  ]*(1-dz)*(1-dy)*(1-dx) +\
             data[z0+1,y0  ,x0  ]*dz*(1-dy)*(1-dx) +\
             data[z0  ,y0+1,x0  ]*(1-dz)*dy*(1-dx) +\
             data[z0  ,y0  ,x0+1]*(1-dz)*(1-dy)*dx +\
             data[z0+1,y0+1,x0  ]*dz*dy*(1-dx) +\
             data[z0+1,y0  ,x0+1]*dz*(1-dy)*dx +\
             data[z0  ,y0+1,x0+1]*(1-dz)*dy*dx +\
             data[z0+1,y0+1,x0+1]*dz*dy*dx
 
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
    if domain['shear']:
        yidx=np.remainder(iarr[1]+xdiv*joffset,Ny)
    else:
        yidx=np.remainder(iarr[1],Ny)
    zidx=iarr[2]

    x0 = xidx.astype(np.intp)
    y0 = yidx.astype(np.intp)
    z0 = zidx.astype(np.intp)
    dx = xidx - x0
    dy = yidx - y0
    dz = zidx - z0

    return x0,y0,z0,dx,dy,dz,vy0

def extend_data(domain,data,smax):
    joffset=get_joffset(domain)

    dslicey=data[:,0,:]
    newdata=np.concatenate((data,dslicey[:,np.newaxis,:]),axis=1)
    d1=np.roll(newdata[:,:,0],-joffset.astype(np.int),axis=1)
    d2=np.roll(newdata[:,:,0],-(joffset.astype(np.int)+1),axis=1)
    dj=joffset-joffset.astype(np.int)
    dslicex=d1*(1-dj)+d2*dj

    newdata=np.concatenate((newdata,dslicex[:,:,np.newaxis]),axis=2)

    Nz, Ny, Nx = newdata.shape
    New_Nz=smax/domain['dx'][2]

    if New_Nz>(Nz/2):
        zeros=np.zeros((int(New_Nz-Nz/2),Ny,Nx))
        newdata=np.concatenate((zeros,newdata,zeros),axis=0)
    return newdata

def los_dump_proj(domain,vec_field,smax,Nside=4,center=[0.,0.,0.],force_write=False,ext='.npy'):

    losdir=domain['losdir']
    step=domain['step']
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s-%d/Nside%d-%s' % (losdir,step,smax,Nside,cstring)
 
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
