from mpi4py import MPI

import matplotlib as mpl
mpl.use('agg')
import numpy as np
import sys
sys.path.insert(0,'/mnt/home/ckim/Sources/pyathena-TIGRESS')

import pyathena.synthetic_observations as syn
import os,glob
import pandas as pd
from astropy.io import fits

def to_map(IQU,kmin=0,kmax=-1):
    Nz,Ny,Nx=IQU[0].shape
    if kmax == -1 : kmax=Nz
    maps=(IQU[0][kmin:kmax,:,:].sum(axis=0),
          IQU[1][kmin:kmax,:,:].sum(axis=0),
          IQU[2][kmin:kmax,:,:].sum(axis=0),)
    return maps

def create_fits(data,domain,fieldname):
    hdr = fits.Header()
    hdr['field']=fieldname
    hdr['time']=domain['time']
    hdr['xmin']=(domain['left_edge'][0],'pc')
    hdr['xmax']=(domain['right_edge'][0],'pc')
    hdr['ymin']=(domain['left_edge'][1],'pc')
    hdr['ymax']=(domain['right_edge'][1],'pc')
    hdr['zmin']=(domain['left_edge'][2],'pc')
    hdr['zmax']=(domain['right_edge'][2],'pc')
    hdr['dx']=(domain['dx'][0],'pc')
    hdr['dy']=(domain['dx'][1],'pc')
    hdr['dz']=(domain['dx'][2],'pc')
    hdu = fits.PrimaryHDU(data,header=hdr)
    return hdu

base='/mnt/ceph/users/ckim/'
pid='MHD_4pc_new'

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

tstart=300
dt=100/size
istart=rank*dt+tstart
iend=(rank+1)*dt+tstart

proj_dir='{}{}/maps-XYproj/'.format(base,pid)

if not os.path.isdir(proj_dir): os.mkdir(proj_dir)

for itime in range(istart,iend,1):
    fname='%s%s/id0/%s.%4.4d.vtk' % (base,pid,pid,itime)
    print fname
    ds,domain=syn.setup_domain(fname,vel=False)

    fields=['density','magnetic_field1','magnetic_field2','magnetic_field3']

    losdata=[]
    for f in fields:
        print 'reading {} ...'.format(f)
        data=syn.read_data(ds,f,domain)
        losdata.append(data)

    hNz=domain['Nx'][2]/2
    hLz=domain['Lx'][2]/2.
    dz=domain['dx'][2]

    IQU=syn.calc_IQU_XY(losdata,domain,dz)

    for kmin in range(hNz,hNz+1024/int(dz)+1,128/int(dz)):
        z0=kmin*dz-hLz
        I,Q,U=to_map(IQU,kmin=kmin)
        domain['left_edge'][2] = z0
        proj_dir_zmin='{}/zmin{:d}/'.format(proj_dir,int(z0))
        if not os.path.isdir(proj_dir_zmin): os.mkdir(proj_dir_zmin)
        for field, fdata in zip(['I','Q','U'],[I,Q,U]):
            hdu = create_fits(fdata,domain,field)
            fbase = os.path.basename(fname)
            fitsname=proj_dir_zmin+fbase.replace('vtk','{}.fits'.format(field))
            hdu.writeto(fitsname,overwrite=True)
