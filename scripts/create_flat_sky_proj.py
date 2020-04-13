# %load /tigress/changgoo/pyathena-TIGRESS/scripts/create_XYsyn.py
from __future__ import print_function
#from mpi4py import MPI

import numpy as np
import sys

import pyathena as pa
import pyathena.synthetic_observations as syn
import os,glob
import pandas as pd
from astropy.io import fits

import gc

def to_map_alongy(IQU,jmin=0,jmax=-1):
    Nz,Ny,Nx=IQU[0].shape
    if jmax == -1 : jmax=Ny
    maps=(IQU[0][:,jmin:jmax,:].sum(axis=1),
          IQU[1][:,jmin:jmax,:].sum(axis=1),
          IQU[2][:,jmin:jmax,:].sum(axis=1),)
    return maps

def to_map_alongx(IQU,imin=0,imax=-1):
    Nz,Ny,Nx=IQU[0].shape
    if imax == -1 : imax=Nx
    maps=(IQU[0][:,:,imin:imax].sum(axis=2),
          IQU[1][:,:,imin:imax].sum(axis=2),
          IQU[2][:,:,imin:imax].sum(axis=2),)
    return maps

def to_map_alongz(IQU,kmin=0,kmax=-1):
    Nz,Ny,Nx=IQU[0].shape
    if kmax == -1 : kmax=Nz
    maps=(IQU[0][kmin:kmax,:,:].sum(axis=0),
          IQU[1][kmin:kmax,:,:].sum(axis=0),
          IQU[2][kmin:kmax,:,:].sum(axis=0),)
    return maps

def create_fits(domain):
    hdr = fits.Header()
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
    hdu = fits.PrimaryHDU(header=hdr)
    return hdu

def add_header_for_glue(hdu,hdr,axis='xyz'):
    for i,ax in enumerate(axis):
        hdu.header['CDELT{}'.format(i+1)]=hdr['d{}'.format(ax)]
        hdu.header['CTYPE{}'.format(i+1)]=ax
        hdu.header['CUNIT{}'.format(i+1)]=hdr.comments['d{}'.format(ax)]
        hdu.header['CRVAL{}'.format(i+1)]=hdr['{}min'.format(ax)]
        hdu.header['CRPIX{}'.format(i+1)]=hdr['{}max'.format(ax)]+hdr['{}min'.format(ax)]
    return 

def flat_sky_proj_kmin(pid,base='/tigress/changgoo/',istart=300,iend=301):
    '''
        flat sky projection along the z axis down to various zmin
    '''
    proj_dir='{}{}/maps-XYproj/'.format(base,pid)

    if not os.path.isdir(proj_dir): os.mkdir(proj_dir)

    vchannel=np.linspace(-100,100,201)

    for itime in range(istart,iend,1):
        fname='%s%s/id0/%s.%4.4d.vtk' % (base,pid,pid,itime)
        print('*** beginning projection with {} ***'.format(fname))
        ds,domain=syn.setup_domain(fname,vel=False)
        x,y,z,=pa.cc_arr(domain)

        fields=['density','magnetic_field1','magnetic_field2','magnetic_field3']

        hNz=domain['Nx'][2]/2
        hNy=domain['Nx'][1]/2
        hNx=domain['Nx'][0]/2
        hLz=domain['Lx'][2]/2.
        dx=domain['dx'][0]
        dy=domain['dx'][1]
        dz=domain['dx'][2]

        losdata=[]
        for f in fields:
            print('*** reading {} ...'.format(f))
            data=syn.read_data(ds,f,domain)
            losdata.append(data)

        nH=losdata[0]
        temp=syn.read_data(ds,'temperature',domain)
        vlos=syn.read_data(ds,'velocity3',domain)
        nH[(temp > 2.e4) | (temp < 10)] = 0.

        hNz=domain['Nx'][2]/2
        hLz=domain['Lx'][2]/2.
        dz=domain['dx'][2]

        IQU=syn.calc_IQU_XY(losdata,domain,dz)

        for kmin in range(hNz,hNz+256/int(dz)+1,128/int(dz)):
            # zmin
            z0=kmin*dz-hLz
            print('*** projectiong from {} to {} '.format(z0,domain['right_edge'][2]))

            domain['left_edge'][2] = z0
            proj_dir_zmin='{}/zmin{:d}/'.format(proj_dir,int(z0))
            if not os.path.isdir(proj_dir_zmin): os.mkdir(proj_dir_zmin)

            # save HI data
            TB,tau=syn.los_to_HI_small_mem(nH[kmin:,:,:],temp[kmin:,:,:],vlos[kmin:,:,:],
                                           vchannel,deltas=dz,los_axis=0)

            hdul = fits.HDUList()
            hdu = create_fits(domain)

            hdu.header['vmin']=(vchannel.min(),'km/s')
            hdu.header['vmax']=(vchannel.max(),'km/s')
            hdu.header['dv']=(vchannel[1]-vchannel[0],'km/s')

            hdul.append(hdu)
            for fdata,label in zip([TB,tau],['TB','tau']):
                hdul.append(fits.ImageHDU(name=label,data=fdata))

            hdr=hdu.header
            for hdu in hdul:
                add_header_for_glue(hdu,hdr,axis='xyv')

            fbase = os.path.basename(fname)
            fitsname=proj_dir_zmin+fbase.replace('vtk','xy.HI.fits')
            hdul.writeto(fitsname,overwrite=True)

            # save IQU data

            I,Q,U=to_map_alongz(IQU,kmin=kmin)

            hdul = fits.HDUList()
            hdu = create_fits(domain)
            hdul.append(hdu)
            for label,fdata in zip(['I','Q','U'],[I,Q,U]):
                hdul.append(fits.ImageHDU(name=label,data=fdata))

            hdr=hdu.header
            for hdu in hdul:
                add_header_for_glue(hdu,hdr,axis='xy')

            fbase = os.path.basename(fname)
            fitsname=proj_dir_zmin+fbase.replace('vtk','xy.IQU.fits')

            hdul.writeto(fitsname,overwrite=True)
    print('*** DONE: synthesized XY projection for {} ***'.format(pid))

def flat_sky_proj_midplane(pid,base='/tigress/changgoo/',istart=300,iend=301):
    '''
        flat sky projection along all x,y,z axes using data cube near the midplane
    '''

    proj_dir='{}{}/maps-XZproj/'.format(base,pid)

    if not os.path.isdir(proj_dir): os.mkdir(proj_dir)

    vchannel=np.linspace(-100,100,201)

    for itime in range(istart,iend,1):
        fname='%s%s/id0/%s.%4.4d.vtk' % (base,pid,pid,itime)
        print('*** beginning projection with {} ***'.format(fname))
        ds,domain=syn.setup_domain(fname,vel=False)
        x,y,z,=pa.cc_arr(domain)

        fields=['density','magnetic_field1','magnetic_field2','magnetic_field3']

        hNz=domain['Nx'][2]/2
        hNy=domain['Nx'][1]/2
        hNx=domain['Nx'][0]/2
        hLz=domain['Lx'][2]/2.
        dx=domain['dx'][0]
        dy=domain['dx'][1]
        dz=domain['dx'][2]
        zcut=z[hNz-hNx:hNz+hNx]
        domain['left_edge'][2] = zcut[0]-dz/2
        domain['right_edge'][2] = zcut[-1]+dz/2

        losdata=[]
        for f in fields:
            data=syn.read_data(ds,f,domain)
            data=data[hNz-hNx:hNz+hNx,:,:]
            print('*** reading {} ...'.format(f))
            losdata.append(data)

        nH=losdata[0]
        temp=syn.read_data(ds,'temperature',domain)[hNz-hNx:hNz+hNx,:,:]
        nH[(temp > 2.e4) | (temp < 10)] = 0.

        # save sim data
        hdul_sim = fits.HDUList()
        hdu_sim = create_fits(domain)

        hdul_sim.append(hdu_sim)
        for fdata,label in zip([nH,temp],['nH','Temp']):
            hdul_sim.append(fits.ImageHDU(name=label,data=fdata))

        for proj in ['xz','yz','xy']:

            if proj == 'xz':
                vlos=syn.read_data(ds,'velocity2',domain,vy0_subtract=False)[hNz-hNx:hNz+hNx,:,:]
                IQU=syn.calc_IQU_XZ(losdata,domain,dx)
                TB,tau=syn.los_to_HI_small_mem(losdata[0],temp,vlos,vchannel,deltas=dy,los_axis=1)
                hdul_sim.append(fits.ImageHDU(name='vy',data=vlos))
            elif proj == 'yz':
                vlos=syn.read_data(ds,'velocity1',domain)[hNz-hNx:hNz+hNx,:,:]
                IQU=syn.calc_IQU_YZ(losdata,domain,dx)
                TB,tau=syn.los_to_HI_small_mem(losdata[0],temp,vlos,vchannel,deltas=dx,los_axis=2)
                hdul_sim.append(fits.ImageHDU(name='vx',data=vlos))
            elif proj == 'xy':
                vlos=syn.read_data(ds,'velocity3',domain)[hNz-hNx:hNz+hNx,:,:]
                IQU=syn.calc_IQU_XY(losdata,domain,dx)
                TB,tau=syn.los_to_HI_small_mem(losdata[0],temp,vlos,vchannel,deltas=dz,los_axis=0)
                hdul_sim.append(fits.ImageHDU(name='vz',data=vlos))

            # save HI data
            hdul = fits.HDUList()
            hdu = create_fits(domain)

            hdu.header['vmin']=(vchannel.min(),'km/s')
            hdu.header['vmax']=(vchannel.max(),'km/s')
            hdu.header['dv']=(vchannel[1]-vchannel[0],'km/s')

            hdul.append(hdu)
            for fdata,label in zip([TB,tau],['TB','tau']):
                hdul.append(fits.ImageHDU(name=label,data=fdata))

            hdr=hdu.header
            for hdu in hdul:
                add_header_for_glue(hdu,hdr,axis='{}v'.format(proj))

            fbase = os.path.basename(fname)

            fitsname=proj_dir+fbase.replace('vtk','{}.HI.fits'.format(proj))
            hdul.writeto(fitsname,overwrite=True)

            # save IQU data
            if proj == 'xz':
                I,Q,U=to_map_alongy(IQU)
            elif proj == 'yz':
                I,Q,U=to_map_alongx(IQU)
            elif proj == 'xy':
                I,Q,U=to_map_alongz(IQU)

            hdul = fits.HDUList()
            hdu = create_fits(domain)
            hdul.append(hdu)
            for label,fdata in zip(['I','Q','U'],[I,Q,U]):
                hdul.append(fits.ImageHDU(name=label,data=fdata))

            hdr=hdu.header
            for hdu in hdul:
                add_header_for_glue(hdu,hdr,axis=proj)

            fbase = os.path.basename(fname)
            fitsname=proj_dir+fbase.replace('vtk','{}.IQU.fits'.format(proj))
            hdul.writeto(fitsname,overwrite=True)

        hdr_sim=hdu_sim.header
        for hdu_sim in hdul_sim:
            add_header_for_glue(hdu_sim,hdr_sim,axis='xyz')

        fbase = os.path.basename(fname)

        fitsname=proj_dir+fbase.replace('vtk','sim.fits')
        hdul_sim.writeto(fitsname,overwrite=True)
        print('*** DONE: synthesized Z projections for {} ***'.format(pid))

def flat_sky_proj_with_pp(pid,mybase='/tigress/changgoo/',
                          ppbase='/tigress/jk11/radps_postproc/',
                          xymax='.xymax1024',
                          itime=0,dt=1,overwrite=False):
    '''
        create synthetic HI projected along the z axis using post-processed data
    '''

    proj_dir='{}{}/maps-XYproj-pp/'.format(mybase,pid)

    if not os.path.isdir(proj_dir): os.mkdir(proj_dir)

    vchannel=np.linspace(-100,100,201)

    fname='%s%s%s/vtk/%s.%4.4d.vtk' % (ppbase,pid,xymax,pid,itime)
    if not os.path.isfile(fname): return
    print('*** beginning projection with {} ***'.format(fname))
    ds,domain=syn.setup_domain(fname,vel=False,shear=False,mhd=False)
    print(domain['time'])
    fitsname=proj_dir+pid+'.{:04d}'.format(round(domain['time']/dt))+xymax+'.HI.fits'
    if (not overwrite) and os.path.isfile(fitsname): return
    x,y,z,=pa.cc_arr(domain)

    losdata=[]

    print('*** reading density ...')
    nden=syn.read_data(ds,'density',domain)
    xH=syn.read_data(ds,'specific_scalar{}'.format(ds.domain['nscal']-1),domain)
    nH=nden*xH
    print('*** reading temperature ...')
    temp=syn.read_data(ds,'temperature',domain)
    print('*** reading LOS velocity ...')
    vlos=syn.read_data(ds,'velocity3',domain)

    dz=domain['dx'][2]
    # create synthetic HI data
    TB,tau=syn.los_to_HI_axis_proj(nH,temp,vlos,vchannel,deltas=dz,los_axis=0,memlim=-1)

    # save to FITS
        
    hdul = fits.HDUList()
    hdu = create_fits(domain)

    hdu.header['vmin']=(vchannel.min(),'km/s')
    hdu.header['vmax']=(vchannel.max(),'km/s')
    hdu.header['dv']=(vchannel[1]-vchannel[0],'km/s')

    hdul.append(hdu)
    for fdata,label in zip([TB,tau],['TB','tau']):
        hdul.append(fits.ImageHDU(name=label,data=fdata))

    hdr=hdu.header
    for hdu in hdul:
        add_header_for_glue(hdu,hdr,axis='xyv')

    fbase = os.path.basename(fname)
    hdul.writeto(fitsname,overwrite=overwrite)

    print('*** DONE: synthesized XY projection for {} ***'.format(fitsname))

def to_fits(pid,base='/tigress/changgoo/',istart=300,iend=301,midplane_cut=True):
    '''
        save raw 3D data in the FITS format
    '''
    if midplane_cut:
        proj_dir='{}{}/fits-midplane/'.format(base,pid)
    else:
        proj_dir='{}{}/fits/'.format(base,pid)

    if not os.path.isdir(proj_dir): os.mkdir(proj_dir)

    for itime in range(istart,iend,1):
        fname='%s%s/id0/%s.%4.4d.vtk' % (base,pid,pid,itime)
        print('*** Reading data from {} ***'.format(fname))
        ds,domain=syn.setup_domain(fname,vel=False)
        x,y,z,=pa.cc_arr(domain)

        fields=['density','temperature',
                'magnetic_field1','magnetic_field2','magnetic_field3',
                'velocity1','velocity2','velocity3']

        hNz=domain['Nx'][2]/2
        hNy=domain['Nx'][1]/2
        hNx=domain['Nx'][0]/2
        hLz=domain['Lx'][2]/2.
        dx=domain['dx'][0]
        dy=domain['dx'][1]
        dz=domain['dx'][2]
        if midplane_cut:
            zcut=z[hNz-hNx:hNz+hNx]
            domain['left_edge'][2] = zcut[0]-dz/2
            domain['right_edge'][2] = zcut[-1]+dz/2

        losdata=[]
        for f in fields:
            data=syn.read_data(ds,f,domain,vy0_subtract=False)
            if midplane_cut:
                data=data[hNz-hNx:hNz+hNx,:,:]
            print('*** reading {} ...'.format(f))
            losdata.append(data)

        # save nT data
        print('*** saving nH and T for {} ...'.format(pid))
        hdul_sim = fits.HDUList()
        hdu_sim = create_fits(domain)

        hdul_sim.append(hdu_sim)
        for fdata,label in zip(losdata[:2],['nH','Temp']):
            hdul_sim.append(fits.ImageHDU(name=label,data=fdata))

        hdr_sim=hdu_sim.header
        for hdu_sim in hdul_sim:
            add_header_for_glue(hdu_sim,hdr_sim,axis='xyz')

        fbase = os.path.basename(fname)

        fitsname=proj_dir+fbase.replace('vtk','nT.fits')
        hdul_sim.writeto(fitsname,overwrite=True)

        # save B data
        print('*** saving B for {} ...'.format(pid))
        hdul_sim = fits.HDUList()
        hdu_sim = create_fits(domain)

        hdul_sim.append(hdu_sim)
        for fdata,label in zip(losdata[2:5],['Bx','By','Bz']):
            hdul_sim.append(fits.ImageHDU(name=label,data=fdata))

        hdr_sim=hdu_sim.header
        for hdu_sim in hdul_sim:
            add_header_for_glue(hdu_sim,hdr_sim,axis='xyz')

        fbase = os.path.basename(fname)

        fitsname=proj_dir+fbase.replace('vtk','B.fits')
        hdul_sim.writeto(fitsname,overwrite=True)

        # save v data
        print('*** saving V for {} ...'.format(pid))
        hdul_sim = fits.HDUList()
        hdu_sim = create_fits(domain)

        hdul_sim.append(hdu_sim)
        for fdata,label in zip(losdata[5:8],['vx','vy','vz']):
            hdul_sim.append(fits.ImageHDU(name=label,data=fdata))

        hdr_sim=hdu_sim.header
        for hdu_sim in hdul_sim:
            add_header_for_glue(hdu_sim,hdr_sim,axis='xyz')

        fbase = os.path.basename(fname)

        fitsname=proj_dir+fbase.replace('vtk','V.fits')
        hdul_sim.writeto(fitsname,overwrite=True)
        print('*** DONE ***')



if __name__ == '__main__':
    #base='/tigress/changgoo/'
    #pid=sys.argv[1]

    #comm = MPI.COMM_WORLD
    #size = comm.Get_size()
    #rank = comm.Get_rank()

    #tstart=300
    #dt=100/size
    #istart=rank*dt+tstart
    #iend=(rank+1)*dt+tstart
    #iend=301

#    flat_sky_proj_with_pp('R8_8pc_rst',xymax='.xymax2048',itime=26)
#    for itime in np.arange(6,400,10):
    for itime in np.arange(101,500,10):
        flat_sky_proj_with_pp('R8_4pc_newacc',xymax='.xymax1024',itime=itime)
        #flat_sky_proj_with_pp('R8_8pc_rst',xymax='.xymax1024',itime=itime)
        gc.collect()

#flat_sky_proj_kmin(pid,base=base,istart=istart,iend=iend)
#flat_sky_proj_midplane(pid,base=base,istart=istart,iend=iend)
#to_fits(pid,base=base,istart=istart,iend=iend,midplane_cut=False)
