import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
sys.path.pop(1)
sys.path.insert(0,'../')
from matplotlib.colors import LogNorm
import pyathena as pa
coolftn=pa.coolftn()

import pyathena.synthetic_observations as syn
import healpy as hp
import glob
import os

## global parameters
base='/tigress/changgoo/'
id='MHD_2pc_S'
Nside=512
npix=hp.nside2npix(Nside)
write=False
center=[[0.,0.,0.],[256.,256.,0.],[-256.,256.,0.],[-256.,-256.,0.],[256.,-256.,0.]]
ifile=12
mhd=True
##

fnames=glob.glob('%s%s/id0/%s.????.vtk' % (base,id,id))
rstfnames=glob.glob('%s%s/id0/%s.????.rst' % (base,id,id))+glob.glob('%s%s/rst/%s.????.rst' % (base,id,id))
fnames.sort()
par,blocks,fields=pa.parse_par(rstfnames[0])

ds=pa.AthenaDataSet(fnames[ifile],ds=None)
print ds.domain['time']

domain=ds.domain
domain['qshear']=eval(par['problem']['qshear'][0])
domain['Omega']=eval(par['problem']['Omega'][0])

fields=['density','pressure','velocity1','velocity2','velocity3']
if mhd: 
    fields.append('magnetic_field1')
    fields.append('magnetic_field2')
    fields.append('magnetic_field3')

data={}
for f in fields:
    data[f] = ds.read_all_data(f)
    
r3d,x3d,y3d,z3d=pa.pos3d(ds.domain)
vy0=-domain['qshear']*domain['Omega']*x3d
data['velocity2'] -= vy0
data['temperature']=coolftn.get_temp(data['T1'])

fields.append('temperature')
print fields


deltas=domain['dx'][2]/2.
for i,c in enumerate(center):
    losdir='%s/dustpol/%s_%d_c%d/' % (ds.dir,ds.step,Nside,i)
    if not os.path.isdir(losdir): os.mkdir(losdir)
    Imap=[]
    Qmap=[]
    Umap=[]
    for ipix in xrange(npix):
        angle = np.rad2deg(hp.pix2ang(Nside,ipix))
        if ipix % Nside == 0: print '%d of %d' % (ipix,npix)
        if np.abs(90-angle[0]) > 5:
            los=syn.get_los(data,domain,Nside,ipix,deltas=deltas,center=c)
            if write:
                pd.DataFrame(los).to_pickle('%s%s.i%d.p' % (losdir,ds.step,ipix))
                
            I,Q,U=syn.los_to_dustpol(los)
            Imap.append(I.value)
            Qmap.append(Q.value)
            Umap.append(U.value)
        else:
            Imap.append(hp.UNSEEN)
            Qmap.append(hp.UNSEEN)
            Umap.append(hp.UNSEEN)
    Imap=np.array(Imap)
    Qmap=np.array(Qmap)
    Umap=np.array(Umap)
    fitsname='%s/dustpol/%s.%s.n%d.c%d.fits' % (ds.dir,ds.id,ds.step,Nside,i)
    print fitsname
    hp.write_map(fitsname,[Imap,Qmap,Umap],column_units=['MJy/sr','MJy/sr','MJy/sr'])
