from __future__ import print_function
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from matplotlib.colors import LogNorm

import pyathena as pa
import pyathena.synthetic_observations as syn

import glob,os

base='/tigress/changgoo/'
sourcedir='../'
cmap=syn.load_planck_cmap('{}/misc/Planck_Parchment_RGB.txt'.format(sourcedir))
cmap.set_bad('white',1.)
plt.register_cmap(cmap=cmap)
plt.rcParams['image.cmap'] = 'planck'
plt.rcParams['font.size']=20
plt.rcParams['figure.figsize']=(12,6)
plt.rcParams['figure.dpi']=150

# magnetar
import sys
sys.path.append('{}/Sources/magnetar/'.format(base))
from bvisual import lic

from astropy.io import fits
from astropy.wcs import WCS

from reproject import reproject_from_healpix, reproject_to_healpix, reproject_interp


car_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                 800
NAXIS2  =                  400
CTYPE1  = 'GLON-CAR'
CRPIX1  =                  400
CRVAL1  =                0.0
CDELT1  =                 -0.45
CUNIT1  = 'deg     '
CTYPE2  = 'GLAT-CAR'
CRPIX2  =                  200
CRVAL2  =                  0.0
CDELT2  =                  0.45
CUNIT2  = 'deg     '
COORDSYS= 'Galactic'
""", sep='\n')

mol_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                 800
NAXIS2  =                  400
CTYPE1  = 'GLON-MOL'
CRPIX1  =                  400
CRVAL1  =                0.0
CDELT1  =                 -0.45
CUNIT1  = 'deg     '
CTYPE2  = 'GLAT-MOL'
CRPIX2  =                  200
CRVAL2  =                  0.0
CDELT2  =                  0.45
CUNIT2  = 'deg     '
COORDSYS= 'Galactic'
""", sep='\n')

def healpix_to_carteian(I,psi):
    xsize = 800
    ysize = xsize/2

    Nside=hp.get_nside(I)
    theta = np.linspace(np.pi, 0, ysize)
    phi   = np.linspace(-np.pi, np.pi, xsize)

    # project the map to a rectangular matrix xsize x ysize
    PHI, THETA = np.meshgrid(phi, theta)
    grid_pix = hp.ang2pix(Nside, THETA, PHI)

    I_car=I[grid_pix]
    psi_car=psi[grid_pix]
    
    return I_car,psi_car


def calc_lic(psi):
    By=np.cos(psi+np.pi/2)
    Bx=np.sin(psi+np.pi/2)
    sz=np.shape(psi)
    length=int(0.1*sz[0])
    licmap=lic(Bx, By, length=length, niter=2)
    
    return licmap

def IQU_to_lic(fitsname):
    I,Q,U=hp.read_map(fitsname,field=[0,1,2])
    U = U*2
    psi=np.arctan2(-U,Q)*0.5

    I_car, psi_car=healpix_to_carteian(I,psi)

    I_mol, footprint = reproject_interp((I_car, WCS(car_header)), mol_header)
    psi_mol, footprint = reproject_interp((psi_car, WCS(car_header)), mol_header)

    licmap_mol=calc_lic(psi_mol)
    licmap_car=calc_lic(psi_car)
    
    hdul=fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    hdul.append(fits.ImageHDU(name='Imap',data=I_mol[:,::-1],header=mol_header))
    hdul.append(fits.ImageHDU(name='LIC',data=licmap_mol[:,::-1],header=mol_header))
    hdul.append(fits.ImageHDU(name='Imap',data=I_car[:,::-1],header=car_header))
    hdul.append(fits.ImageHDU(name='LIC',data=licmap_car[:,::-1],header=car_header))

    hdul.writeto(fitsname.replace('.fits','.lic.fits'),overwrite=True)

def draw_lic_maps(lic_fitsname):
    hdul=fits.open(lic_fitsname)
    fig = plt.figure(0,figsize=(8,8))
    ax1 = plt.subplot(211,projection=WCS(hdul[1].header))
    licmap=hdul[2].data
    limits = np.nanmean(licmap.flatten()) + np.array([-1,1])*2*np.nanstd(licmap.flatten())
    im=ax1.imshow(hdul[1].data, origin='lower', cmap='planck',norm=LogNorm())
    ax1.imshow(hdul[2].data, origin='lower', alpha=0.2, cmap='binary', clim=limits)
    ax1.coords.frame.set_color('none')
    
    ax1 = plt.subplot(212,projection=WCS(hdul[3].header))
    im=ax1.imshow(hdul[3].data, origin='lower', cmap='planck',norm=LogNorm())
    ax1.imshow(hdul[4].data, origin='lower', alpha=0.2, cmap='binary', clim=limits)
    ax1.coords.frame.set_color('none')

    fig.savefig(lic_fitsname.replace('fits','png'),dpi=100,bbox_inches='tight')

base='/tigress/changgoo/'
pid='MHD_4pc_new'
fitsfiles=glob.glob('{}/{}/maps/s40_875/{}.03??.Nside128-x0y0z0.fits'.format(base,pid,pid))
fitsfiles.sort()
#fitsname='{}/{}/maps/s40_875/{}.0300.Nside128-x0y0z0.fits'.format(base,pid,pid)

for fitsname in fitsfiles:
    if not os.path.isfile(fitsname.replace('.fits','.lic.fits')):
        print('*** Calculating LIC map for {} ***'.format(fitsname))
        IQU_to_lic(fitsname)
    else:
        print('*** Skipping LIC map for {} ***'.format(fitsname))

    print('*** Drwaing for {} ***'.format(fitsname))
    draw_lic_maps(fitsname.replace('.fits','.lic.fits'))
