# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True

# Cython

from .tools import *
import numpy as np
cimport numpy as np # access to Numpy from Cython layer

cdef void interp3D_cy(double[:,:,::1] input_array,
             double[::1] x_indices,
             double[::1] y_indices,
             double[::1] z_indices,
             double[::1] output,
             int nindices, int Nx, int Ny, int Nz):
    cdef:
        int x0,y0,z0,x1,y1,z1
        double x,y,z
        int i
        
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
                    
def get_los_all(data,domain,nside,ipix,smin=0.,smax=3000.,deltas=1.,center=[0.,0.,0.]):
    """
    cpdef get_los_all(data,domain,nside,ipix,smin=0.,smax=3000.,deltas=1.,center=[0.,0.,0.]):
    """
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

def get_los_one(data,domain,nside,ipix,smin=0.,smax=3000.,deltas=1.,center=[0.,0.,0.]):
    """
    cpdef get_los_one(data,domain,nside,ipix,smin=0.,smax=3000.,deltas=1.,center=[0.,0.,0.]):
    """
    hat=get_hat(nside,ipix)
    idx,sarr = los_idx(hat['Z'],domain,smin=smin,smax=smax,ds=deltas,center=center)

    if 'Omega' in domain:
        joffset=get_joffset(domain)
        sheared_periodic(domain,idx,joffset=joffset)
    else:
        periodic(domain,idx,iaxis=[0,1])

    Nz,Ny,Nx=data.shape
    nidx=len(idx[0])

    cdef:
        double[::1] _xidx= np.array(idx[2], np.float64)
        double[::1] _yidx= np.array(idx[1], np.float64)
        double[::1] _zidx= np.array(idx[0], np.float64)
        double[::1] _output=np.empty(nidx, np.float64)
        double[:,:,::1] _input=np.array(data, np.float64)

    interp3D_cy(_input,_xidx,_yidx,_zidx,_output,nidx,Nx,Ny,Nz)

    return sarr,np.array(_output)
