# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True

import numpy as np
import healpy as hp
import astropy.constants as c
import astropy.units as u

Tdust=18
nu0=353
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
Bnu=blackbody_nu(nu0*u.GHz,Tdust*u.K).to('MJy/sr').value
sigma_353=1.2e-26 # Planck 2013 results XI. by Planck Collaboration XI (2014)
p0=0.2
pc_to_cm=c.pc.cgs.value

def cc_idx(domain,pos):
    le=domain['left_edge']
    dx=domain['dx']
    Nx=domain['Nx']
    if np.array(pos).ndim == 2:
        le=le[:,np.newaxis]
        dx=dx[:,np.newaxis]
        Nx=Nx[:,np.newaxis]
    idx=(pos-le-0.5*dx)/dx
    return idx

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
        idx[idx < 0] += Nx[i]

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
            zarr=zhat*sarr
            ind = np.where((zarr < zmax)*(zarr > zmin))
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

def los_to_dustpol(los,deltas=1.):

    ds=deltas*pc_to_cm
    
    Bperp2=los['magnetic_field_X']**2+los['magnetic_field_Y']**2
    B2=Bperp2+los['magnetic_field_Z']**2
    cos2phi=(los['magnetic_field_X']**2-los['magnetic_field_Y']**2)/Bperp2
    sin2phi=los['magnetic_field_X']*los['magnetic_field_Y']/Bperp2
    cosgam2=Bperp2/B2

    nH=los['density']
    dtau=sigma_353*nH*ds
    tau=dtau.cumsum()
 #   print nH.sum()*ds.cgs,tau[-1], np.exp(-tau[-1])

    I=Bnu*(1.0-p0*(cosgam2-2./3.0))*dtau#*np.exp(-tau)
    Q=p0*Bnu*cos2phi*cosgam2*dtau#*np.exp(-tau)
    U=p0*Bnu*sin2phi*cosgam2*dtau#*np.exp(-tau)
    
    return I.sum(),Q.sum(),U.sum()

# Cython

cimport numpy as np # access to Numpy from Cython layer

cimport cython

@cython.boundscheck(False) # won't check that index is in bounds of array
@cython.wraparound(False) # array[-1] won't work
@cython.nonecheck(False) # variables are never set to None
@cython.cdivision(True) # don't protect against dividing by zero

cdef interp3D_cy(double[:,:,::1] input_array,
             double[::1] x_indices,
             double[::1] y_indices,
             double[::1] z_indices,
             double[::1] output,
             int nindices, int Nx, int Ny, int Nz):
    cdef:
        int x0,y0,z0,x1,y1,z1
        double x,y,z
        size_t i
        
    for i in range(nindices):
        x0 = int(x_indices[i])
        y0 = int(y_indices[i])
        z0 = int(z_indices[i])
        x1 = x0 + 1
        y1 = y0 + 1
        z1 = z0 + 1

#Check if xyz1 is beyond array boundary:
        if x1 == Nx: x1 = 0
        if y1 == Ny: y1 = 0
        if z1 == Nz: z1 = 0
        x = x_indices[i] - x0
        y = y_indices[i] - y0
        z = z_indices[i] - z0
        output[i] = input_array[x0,y0,z0]*(1-x)*(1-y)*(1-z) +\
                    input_array[x1,y0,z0]*x*(1-y)*(1-z) +\
                    input_array[x0,y1,z0]*(1-x)*y*(1-z) +\
                    input_array[x0,y0,z1]*(1-x)*(1-y)*z +\
                    input_array[x1,y0,z1]*x*(1-y)*z +\
                    input_array[x0,y1,z1]*(1-x)*y*z +\
                    input_array[x1,y1,z0]*x*y*(1-z) +\
                    input_array[x1,y1,z1]*x*y*z
                    
cpdef get_los_cy(data,domain,nside,ipix,
                 smin=0.,smax=3000.,deltas=1.,center=[0.,0.,0.]):
    hat=get_hat(nside,ipix)
    idx,sarr = los_idx(hat['Z'],domain,smin=smin,smax=smax,ds=deltas,center=center)
    
    fields=['density','temperature']
    fields += ['velocity1','velocity2','velocity3']
    fields += ['magnetic_field1','magnetic_field2','magnetic_field3']
    if 'Omega' in domain:
        Omega=domain['Omega']
        qshear=domain['qshear'] 
        xarr=sarr*hat['Z'][0]+center[0]
        vy0=-qshear*Omega*xarr
        joffset=get_joffset(domain)
        sheared_periodic(domain,idx,joffset=joffset)
    else:
        periodic(domain,idx,iaxis=[0,1])
# Cythonized part
    cdef:
        int Nx,Ny,Nz,nidx
    Nz,Ny,Nx=data['density'].shape
    nidx=len(idx[0])
    los_out={} 
    
    cdef:
        double[::1] _xidx= np.array(idx[2], np.float64)
        double[::1] _yidx= np.array(idx[1], np.float64)
        double[::1] _zidx= np.array(idx[0], np.float64)
        double[::1] _output=np.empty(nidx, np.float64)
        double[:,:,::1] _input_den=np.array(data['density'], np.float64)
        double[:,:,::1] _input_temp=np.array(data['temperature'], np.float64)

        double[:,:,::1] _input_v1=np.array(data['velocity1'], np.float64)
        double[:,:,::1] _input_v2=np.array(data['velocity2'], np.float64)
        double[:,:,::1] _input_v3=np.array(data['velocity3'], np.float64)

        double[:,:,::1] _input_B1=np.array(data['magnetic_field1'], np.float64)
        double[:,:,::1] _input_B2=np.array(data['magnetic_field2'], np.float64)
        double[:,:,::1] _input_B3=np.array(data['magnetic_field3'], np.float64)

    interp3D_cy(_input_den,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)
    los_out['density']=np.array(_output)
    interp3D_cy(_input_temp,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)
    los_out['temperature']=np.array(_output)
    interp3D_cy(_input_v1,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)
    los_out['velocity1']=np.array(_output)
    interp3D_cy(_input_v2,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)
    los_out['velocity2']=np.array(_output)
    interp3D_cy(_input_v3,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)
    los_out['velocity3']=np.array(_output)
    interp3D_cy(_input_B1,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)
    los_out['magnetic_field1']=np.array(_output)
    interp3D_cy(_input_B2,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)
    los_out['magnetic_field2']=np.array(_output)
    interp3D_cy(_input_B3,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)
    los_out['magnetic_field3']=np.array(_output)
    
# end of Cythoinzed part

    if 'Omega' in domain:
        los_out['velocity2'] += vy0
    
    for axis in ['Z','X','Y']:
        los_out['velocity_'+axis] =hat[axis][0]*los_out['velocity1'] 
        los_out['velocity_'+axis]+=hat[axis][1]*los_out['velocity2']
        los_out['velocity_'+axis]+=hat[axis][2]*los_out['velocity3']
        los_out['magnetic_field_'+axis] =hat[axis][0]*los_out['magnetic_field1'] 
        los_out['magnetic_field_'+axis]+=hat[axis][1]*los_out['magnetic_field2']
        los_out['magnetic_field_'+axis]+=hat[axis][2]*los_out['magnetic_field3']
    los_out['sarr']=sarr
    return los_out
