from mpi4py import MPI

import matplotlib as mpl
mpl.use('agg')
import numpy as np
import sys
sys.path.insert(0,'/mnt/home/ckim/Sources/pyathena-TIGRESS')

import pyathena.synthetic_observations as syn
from pyathena.synthetic_observations.los_to_HI import *

import healpy as hp
import os

import pandas as pd
center=[0,0,0]
base='/mnt/ceph/users/ckim/'
pid='MHD_4pc_new'
Nside=128

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

tstart=300
dt=100/size
istart=rank*dt+tstart
iend=(rank+1)*dt+tstart

vchannel=np.linspace(-100,100,201)
#for itime in range(istart,iend,1):
#    los_dir='{}{}/los/{:04d}-3500/Nside{}-x{}y{}z{}/'.format(base,pid,itime,Nside,center[0],center[1],center[2])
#
#    if os.path.isfile('{}/vlos.npy'.format(los_dir)):
#        print los_dir
#        den=np.load('{}/density.npy'.format(los_dir))[:,40:875]
#        temp=np.load('{}/temperature.npy'.format(los_dir))[:,40:875]
#        vlos=np.load('{}/vlos.npy'.format(los_dir))[:,40:875]
#        idx=(temp > 2.e4) 
#        den[idx] = 0.
#
#        TB,tau=los_to_HI_small_mem(den,temp,vlos,vchannel,deltas=4)
#        np.save('{}/TB.npy'.format(los_dir),TB)
#        np.save('{}/tau.npy'.format(los_dir),tau)

for itime in range(istart,iend,1):
    los_dir='{}{}/los/{:04d}-3500/Nside{}-x{}y{}z{}/'.format(base,pid,itime,Nside,center[0],center[1],center[2])
    if os.path.isfile('{}/TB.npy'.format(los_dir)):
        print los_dir
        den=np.load('{}/density.npy'.format(los_dir))[:,40:875]
        temp=np.load('{}/temperature.npy'.format(los_dir))[:,40:875]
        vlos=np.load('{}/vlos.npy'.format(los_dir))[:,40:875]
        idx=(temp > 2.e4) 
        den[idx] = 0.

        TB=np.load('{}/TB.npy'.format(los_dir))
        tau=np.load('{}/tau.npy'.format(los_dir))

        ds=4*3.085677581467192e+18
        NH=den.sum(axis=1)*ds
        Tsavg=den.sum(axis=1)/(den/temp).sum(axis=1)

        idx=(temp > 2.e2) 
        den[idx] = 0.
        NHcold=den.sum(axis=1)*ds

        dv=1
        tauint=tau.sum(axis=0)*dv
        Tsobs=TB/(1.0-np.exp(-tau)) 
        Tsobs=np.nansum(tau*Tsobs,axis=0)*dv
        NHobs=1.813e18*Tsobs
        Tsobs=Tsobs/tauint
        NHthin=1.813e18*TB.sum(axis=0)*dv

        maps={'NH':NH,'Ts':Tsavg,'NHcold':NHcold,'Tsobs':Tsobs,\
              'NHobs':NHobs,'NHthin':NHthin}
        maps=pd.DataFrame(maps)
        maps.to_pickle('{}/HI_maps.p'.format(los_dir))
