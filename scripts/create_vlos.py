
import sys
sys.path.insert(0,'/mnt/home/ckim/Sources/pyathena-TIGRESS')

import pyathena.synthetic_observations as syn
from pyathena.synthetic_observations.los_to_HI import *

import healpy as hp
import os
center=[0,0,0]
base='/mnt/ceph/users/ckim/'
pid='MHD_4pc_new'
Nside=128
npix=hp.nside2npix(Nside)
ipix = np.arange(npix)
hat=syn.get_hat(Nside,ipix)
 
for itime in range(300,400):
    los_dir='{}{}/los/{:04d}-3500/Nside{}-x{}y{}z{}/'.format(base,pid,itime,Nside,center[0],center[1],center[2])

    if (os.path.isfile('{}/velocity1.npy'.format(los_dir)) and \
        (not os.path.isfile('{}/vlos.npy'.format(los_dir)))):
        print los_dir
        v1=np.load('{}/velocity1.npy'.format(los_dir))
        v2=np.load('{}/velocity2.npy'.format(los_dir))
        v3=np.load('{}/velocity3.npy'.format(los_dir))
        vlos=hat['Z'][0][:,np.newaxis]*v1 + \
             hat['Z'][1][:,np.newaxis]*v2 + \
             hat['Z'][2][:,np.newaxis]*v3

        np.save('{}/vlos.npy'.format(los_dir),vlos)
