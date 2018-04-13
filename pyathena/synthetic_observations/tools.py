import numpy as np
import healpy as hp
import astropy.constants as c
import astropy.units as u
import os

from ..utils import cc_idx

def make_directory(domain,smax,Nside=4,center=[0.,0.,0.]):
    losdir=domain['losdir']
    step=domain['step']
    stepdir='%s%s-%d' % (losdir,step,smax)
    cstring='x%dy%dz%d' % (center[0],center[1],center[2])
    outdir='%s%s-%d/Nside%d-%s' % (losdir,step,smax,Nside,cstring)
    
    if not os.path.isdir(losdir): os.mkdir(losdir)
    if not os.path.isdir(stepdir): os.mkdir(stepdir)
    if not os.path.isdir(outdir): os.mkdir(outdir)


def sheared_periodic(domain,idx,joffset=0.0):
    Nx=domain['Nx'][0]
    nidx = idx[0] < 0
    pidx = idx[0] >= Nx
    while(nidx.any()):
        idx[1][nidx] -= joffset
        idx[0][nidx] += Nx
        nidx=idx[0]<0
    while(pidx.any()):
        idx[1][pidx] += joffset
        idx[0][pidx] -= Nx
        pidx = idx[0] >= Nx
    periodic(domain,idx,iaxis=[1])

def periodic(domain,idx,iaxis=[0,1,2]):
    Nx=domain['Nx']
    for i in iaxis:
        idx[i]=np.fmod(idx[i],Nx[i])
        idx[i][idx[i] < 0] += Nx[i]

def los_idx(hat,domain,smin=0.,smax=3000.,ds=1.,
            center=[0,0,0],vectorize=True,zmax_cut=True):
    zmax=domain['right_edge'][2]-0.5*domain['dx'][2]
    zmin=domain['left_edge'][2]+0.5*domain['dx'][2]
    xhat=hat[0]
    yhat=hat[1]
    zhat=hat[2]
# vectorlized version
#  smax = np.sqrt(smax**2-zmax**2)
    if vectorize:
        sarr=np.arange(smin,smax,ds)
        if zmax_cut:
            zarr=zhat*sarr + center[2]
            ind = np.where((zarr < zmax) & (zarr > zmin))
            sarr = sarr[ind]
        xarr=xhat*sarr + center[0]
        yarr=yhat*sarr + center[1]
        zarr=zhat*sarr + center[2]

        iarr = cc_idx(domain,[xarr,yarr,zarr])
    else:
# while loop for s from 0 to smax
        idx=[]
        s=ds
        sarr=[]
        while s < smax:
# get x,y,z by multiplying the path length s 
            x=xhat*s + center[0]
            y=yhat*s + center[1]
            z=zhat*s + center[2]
            if zmax_cut and abs(z)>zmax: break
# get i,j,k by applying periodic BCs in x and y
            idx.append(cc_idx(domain,[x,y,z]))
            sarr.append(s)
            s = s+ds
# store indices and path lengthes
        iarr=np.array(idx).T
        sarr=np.array(sarr)

    return iarr,sarr

def get_joffset(domain):
    Omega=domain['Omega']
    qshear=domain['qshear'] 
    Lx=domain['right_edge'][0]-domain['left_edge'][0]
    Ly=domain['right_edge'][1]-domain['left_edge'][1]
    qomL=qshear*Omega*Lx
    yshear=qomL*domain['time']
    deltay=np.mod(yshear,Ly)
    joffset=deltay/domain['dx'][1]

    return joffset

def get_hat(nside,ipix):

    theta,phi = hp.pix2ang(nside,ipix)
    rhat=[np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]
    thhat=[np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)]
    phhat=[-np.sin(phi),np.cos(phi),0]
    hat={'Z':rhat,'X':thhat,'Y':phhat}

    return hat

def load_planck_cmap(cmap_fname="../misc/Planck_Parchment_RGB.txt"):
    '''
    https://zonca.github.io/2013/09/Planck-CMB-map-at-high-resolution.html
    '''
    from matplotlib.colors import ListedColormap
    import numpy as np
    colombi1_cmap = ListedColormap(np.loadtxt(cmap_fname)/255.)
    colombi1_cmap.set_bad("gray") # color of missing pixels
    colombi1_cmap.set_under("white") # color of background, necessary if you want to use

    cmap = colombi1_cmap
    return cmap

def get_center(dname):
    cstr=dname.split('-x')[1]
    yidx=cstr.rfind('y')
    zidx=cstr.rfind('z')
    x0=int(cstr[:yidx])
    y0=int(cstr[yidx+1:zidx])
    z0=int(cstr[zidx+1:])
    center=[x0,y0,z0]
    return center
