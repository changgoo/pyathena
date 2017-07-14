def save_to_fits(ds,mhd=True):
    from astropy.io import fits
    import pyathena as pa
    units=pa.set_units(muH=1.4271)
    coolftn=pa.coolftn()
    
    import os
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
        hdu.writeto(fitsname,clobber=True)
