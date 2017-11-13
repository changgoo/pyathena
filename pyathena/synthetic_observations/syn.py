import numpy as np
import argparse
import healpy as hp
import astropy.constants as c
import astropy.units as u
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

def set_domain():
  domain={}
  domain['left_edge']=np.array([-256,-256,-512])
  domain['right_edge']=np.array([256,256,512])
  domain['dx']=np.array([2,2,2])
  domain['Lx']=domain['right_edge']-domain['left_edge']
  domain['Nx']=domain['Lx']/domain['dx']
  return domain

def interp3D(input_array,indices):
  output = np.empty(indices[0].shape)
  x_indices = indices[2]
  y_indices = indices[1]
  z_indices = indices[0]

  x0 = x_indices.astype(np.integer)
  y0 = y_indices.astype(np.integer)
  z0 = z_indices.astype(np.integer)
  x1 = x0 + 1
  y1 = y0 + 1
  z1 = z0 + 1

#Check if xyz1 is beyond array boundary:
  x1[np.where(x1==input_array.shape[0])] = 0
  y1[np.where(y1==input_array.shape[1])] = 0
  z1[np.where(z1==input_array.shape[2])] = 0

  x = x_indices - x0
  y = y_indices - y0
  z = z_indices - z0
  output = (input_array[x0,y0,z0]*(1-x)*(1-y)*(1-z) +
             input_array[x1,y0,z0]*x*(1-y)*(1-z) +
             input_array[x0,y1,z0]*(1-x)*y*(1-z) +
             input_array[x0,y0,z1]*(1-x)*(1-y)*z +
             input_array[x1,y0,z1]*x*(1-y)*z +
             input_array[x0,y1,z1]*(1-x)*y*z +
             input_array[x1,y1,z0]*x*y*(1-z) +
             input_array[x1,y1,z1]*x*y*z)
  return output

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

def get_rslice(data,domain,nside,s0,fields=['density','temperature'],\
  vel=True,mhd=True,smax=3000.,deltas=1.,center=[0.,0.,0.]):
    npix=hp.nside2npix(Nside)
    for ipix in range(npix): 
        hat=get_hat(nside,ipix)
        idx,sarr = los_idx(hat['Z'],domain,smax=smax,ds=deltas,center=center)
    if vel:
        fields += ['velocity1','velocity2','velocity3']
    if mhd:
        fields += ['magnetic_field1','magnetic_field2','magnetic_field3']
    if domain.has_key('Omega'):
        Omega=domain['Omega']
        qshear=domain['qshear'] 
        xarr=sarr*hat['Z'][0]+center[0]
        if vel: vy0=-qshear*Omega*xarr
        joffset=get_joffset(domain)
        sheared_periodic(domain,idx,joffset=joffset)
    else:
        periodic(domain,idx,iaxis=[0,1])
    los_out={} 
    for f in fields:
        los_out[f]=interp3D(data[f],idx)
    if vel and domain.has_key('Omega'):
        los_out['velocity2'] += vy0
    
    if vel:
        for axis in ['Z','X','Y']:
            los_out['velocity_'+axis] =hat[axis][0]*los_out['velocity1'] 
            los_out['velocity_'+axis]+=hat[axis][1]*los_out['velocity2']
            los_out['velocity_'+axis]+=hat[axis][2]*los_out['velocity3']
    if mhd: 
        for axis in ['Z','X','Y']:
            los_out['magnetic_field_'+axis] =hat[axis][0]*los_out['magnetic_field1'] 
            los_out['magnetic_field_'+axis]+=hat[axis][1]*los_out['magnetic_field2']
            los_out['magnetic_field_'+axis]+=hat[axis][2]*los_out['magnetic_field3']
    los_out['sarr']=sarr
    return los_out


def get_los(data,domain,nside,ipix,fields=['density','temperature'],\
  vel=True,mhd=True,smin=0.,smax=3000.,deltas=1.,center=[0.,0.,0.]):
    hat=get_hat(nside,ipix)
    idx,sarr = los_idx(hat['Z'],domain,smin=smin,smax=smax,ds=deltas,center=center)
    if vel:
        fields += ['velocity1','velocity2','velocity3']
    if mhd:
        fields += ['magnetic_field1','magnetic_field2','magnetic_field3']
    if domain.has_key('Omega'):
        Omega=domain['Omega']
        qshear=domain['qshear'] 
        xarr=sarr*hat['Z'][0]+center[0]
        if vel: vy0=-qshear*Omega*xarr
        joffset=get_joffset(domain)
        sheared_periodic(domain,idx,joffset=joffset)
    else:
        periodic(domain,idx,iaxis=[0,1])
    los_out={} 
    for f in fields:
        los_out[f]=interp3D(data[f],idx)
    if vel and domain.has_key('Omega'):
        los_out['velocity2'] += vy0
    
    if vel:
        for axis in ['Z','X','Y']:
            los_out['velocity_'+axis] =hat[axis][0]*los_out['velocity1'] 
            los_out['velocity_'+axis]+=hat[axis][1]*los_out['velocity2']
            los_out['velocity_'+axis]+=hat[axis][2]*los_out['velocity3']
    if mhd: 
        for axis in ['Z','X','Y']:
            los_out['magnetic_field_'+axis] =hat[axis][0]*los_out['magnetic_field1'] 
            los_out['magnetic_field_'+axis]+=hat[axis][1]*los_out['magnetic_field2']
            los_out['magnetic_field_'+axis]+=hat[axis][2]*los_out['magnetic_field3']
    los_out['sarr']=sarr
    return los_out

def los_to_HI(los,vmax=100,dvch=1.0,deltas=1.):
    nu0=1.420406e9*u.Hz
    A10=2.8843e-15/u.s
    mass=c.m_p
    c2=2.0*np.sqrt(np.log(2.0))
    c3=c.c/nu0
    c1=c2/np.sqrt(np.pi)*(c3)

    vchmin=-vmax
    vchmax=vmax
    Nchannel=int((vchmax-vchmin)/dvch)+1
    vchannel=np.tile((np.arange(Nchannel)*dvch+vchmin)[:,np.newaxis],(1,len(los['sarr'])))*u.km/u.s
    
    Tlos=np.tile(los['temperature'][np.newaxis,:],(Nchannel,1))*u.K
    vlos=np.tile(los['velocity_Z'][np.newaxis,:],(Nchannel,1))*u.km/u.s
    nlos=np.tile(los['density'][np.newaxis,:],(Nchannel,1))/u.cm**3
    Tspin=Tspin_WF(Tlos,nlos)
    
    alpha=np.sqrt(c.c**2*mass/(2*c.k_B*Tlos))
    nu_D=nu0/alpha
    nu_L=c2*nu_D
    v_L=c3*nu_L

    c4=3./32./np.pi*A10*c.h*c.c**2/c.k_B/nu0
    ds=deltas

    phi_v=c1/v_L*np.exp(-(c2*(vchannel-vlos)/v_L)**2)
    kappa_v=c4*nlos/Tspin*phi_v
    tau_los=kappa_v*ds*c.pc
    tau_cumul=tau_los.cumsum(axis=1)
    TB=np.nansum(Tspin*(1-np.exp(-tau_los))*np.exp(-tau_cumul),axis=1).cgs
    TBthin=np.nansum(Tspin*tau_los,axis=1).cgs
    tau_v=np.nansum(kappa_v*ds*c.pc,axis=1).cgs
    return {'TB':TB,'TBthin':TBthin,'tau':tau_v,'vchannel':vchannel[:,0]}

def los_to_dustpol(los,Tdust=18,nu0=353,deltas=1.):
    from astropy.analytic_functions import blackbody_lambda, blackbody_nu
    Bnu=blackbody_nu(nu0*u.GHz,Tdust*u.K)
    sigma_353=1.2e-26*u.cm**2 # Planck 2013 results XI. by Planck Collaboration XI (2014)
    p0=0.2
    ds=deltas*c.pc
    
    Bperp2=los['magnetic_field_X']**2+los['magnetic_field_Y']**2
    B2=Bperp2+los['magnetic_field_Z']**2
    cos2phi=(los['magnetic_field_X']**2-los['magnetic_field_Y']**2)/Bperp2
    sin2phi=los['magnetic_field_X']*los['magnetic_field_Y']/Bperp2
    cosgam2=Bperp2/B2

    nH=los['density']/u.cm**3
    dtau=(sigma_353*nH*ds).cgs
    tau=dtau.cumsum()
 #   print nH.sum()*ds.cgs,tau[-1], np.exp(-tau[-1])

    I=Bnu*(1.0-p0*(cosgam2-2./3.0))*dtau#*np.exp(-tau)
    Q=p0*Bnu*cos2phi*cosgam2*dtau#*np.exp(-tau)
    U=p0*Bnu*sin2phi*cosgam2*dtau#*np.exp(-tau)
    
    return I.sum().to('MJy/sr'),Q.sum().to('MJy/sr'),U.sum().to('MJy/sr')

def k10h(T2):
    k10_1=1.19e-10*T2**(0.74-0.20*np.log(T2))
    k10_2=2.24e-10*T2**0.207*np.exp(-0.876/T2)
    k10=k10_1
    idx=k10_2 > k10_1
    k10[idx]=k10_2[idx]
    return k10*u.cm**3/u.s

def k10e(T2):
    Temp=T2*1.e2
    k10=-9.607+0.5*np.log10(Temp)*np.exp(-(np.log10(Temp))**4.5/1.8e3)
    return 10.**k10*u.cm**3/u.s

def Tspin_WF(Temp,nH,nalpha=1.e-6):

    T2=Temp/(1.e2*u.K)
    nu0=1.420406e9*u.Hz
    A10=2.8843e-15/u.s
    T21=c.h*nu0/c.k_B
    TA=3.77*u.K
    Tk=Temp
    TL=Temp

    k10=k10h(T2)

    yc=(T21*k10*nH/A10/Temp).cgs
    yalpha=(7.76e11*nalpha/(TL/u.K)/np.sqrt(Tk/u.K)).cgs
    Ts=(TA+yc*Temp+yalpha*TL)/(1.0+yc+yalpha)

    return Ts


