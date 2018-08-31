from mpi4py import MPI

import pyathena.synthetic_observations as syn
from pyathena.synthetic_observations.los_dump_all import los_dump

import numpy as np
import os,glob


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

tstart=320
dt=20
istart=rank*dt+tstart
iend=(rank+1)*dt+tstart

base='/mnt/ceph/users/ckim/'
id='MHD_4pc_new'

Nside=128
center=[0,0,0]
smax=3500

print(rank, istart, iend)
for itime in range(istart,iend,1):
    fname='%s%s/id0/%s.%4.4d.vtk' % (base,id,id,itime)
    if os.path.isfile(fname):
        print(rank,fname)
        ds,domain=syn.setup_domain(fname,vel=True)
        deltas=domain['dx'][2]
        syn.make_directory(domain,smax,Nside=Nside,center=center)
        for f in ['density','temperature','velocity1','velocity2','velocity3']:
            los_dump(ds,domain,deltas,smax,[f],Nside,center,False)
