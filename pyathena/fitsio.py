
from astropy.io import fits
import pyathena as pa
import os
import numpy as np
def save_to_fits(ds,mhd=True):
    units=pa.set_units(muH=1.4271)
    coolftn=pa.coolftn()
    
    if not os.path.isdir(ds.dir+'fits'): os.mkdir(ds.dir+'fits')
    fields=['density','temperature','velocity']
    if mhd: fields.append('magnetic_field')
    for field in fields:
        if field is 'temperature':
            data=coolftn.get_temp(ds.read_all_data('T1'))
        else:
            data=ds.read_all_data(field)
        fitsname='%s/fits/%s.%s.%s.fits' % (ds.dir,ds.id,ds.step,field[:4])
        hdr = fits.Header()
        hdr['field']=field
        hdr['time']=ds.domain['time']
        hdr['tMyr']=(ds.domain['time']*units['time'].to('Myr').value,'Myr')
        hdr['xmin']=(ds.domain['left_edge'][0],'pc')
        hdr['xmax']=(ds.domain['right_edge'][0],'pc')
        hdr['ymin']=(ds.domain['left_edge'][1],'pc')
        hdr['ymax']=(ds.domain['right_edge'][1],'pc')
        hdr['zmin']=(ds.domain['left_edge'][2],'pc')
        hdr['zmax']=(ds.domain['right_edge'][2],'pc')
        hdr['dx']=(ds.domain['dx'][0],'pc')
        hdr['dy']=(ds.domain['dx'][1],'pc')
        hdr['dz']=(ds.domain['dx'][2],'pc')
        hdr['unit']=(units[field].value,units[field].unit)
        if ds.domain.has_key('qshear'):
            hdr['qshear']=ds.domain['qshear']
        if ds.domain.has_key('Omega'):
            hdr['Omega']=(ds.domain['Omega'],'km/s/pc')
        hdu = fits.PrimaryHDU(data,header=hdr)
        hdu.writeto(fitsname,overwrite=True)

def get_domain(fitsname,vel=True,mhd=True):
    hdulist=fits.open(fitsname,memmap=True)
    hdr=hdulist[0].header
    domain={}
    domain['time']=hdr['time']
    domain['qshear']=hdr['qshear']
    domain['Omega']=hdr['Omega']
    domain['left_edge']=np.array([hdr['xmin'],hdr['ymin'],hdr['zmin']])
    domain['right_edge']=np.array([hdr['xmax'],hdr['ymax'],hdr['zmax']])
    domain['dx']=np.array([hdr['dx'],hdr['dy'],hdr['dz']])
    domain['Nx']=np.array([hdr['naxis1'],hdr['naxis2'],hdr['naxis3']])
    domain['Lx']=domain['right_edge']-domain['left_edge']
    fields=['density','temperature']
    if vel: fields.append('velocity')
    if mhd: fields.append('magnetic_field')
    hdulist.close()
    return domain,fields

def get_data(fitsname,data={}):
    hdulist=fits.open(fitsname)
    field=hdulist[0].header['field']
    hdulist.readall()
    if field.startswith('velo') or field.startswith('magn'):
        for iaxis in range(3):
            data['%s%d' % (field,iaxis+1)]=hdulist[0].data[:,:,:,iaxis]
    else:
        data[field]=hdulist[0].data
    hdulist.close()

    return data
