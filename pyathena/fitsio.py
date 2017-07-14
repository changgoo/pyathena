def save_to_fits(ds):
    from astropy.io import fits
    import pyathena as pa
    units=pa.set_units(muH=1.4271)
    coolftn=pa.coolftn()
    
    import os
    if not os.path.isdir(ds.dir+'fits'): os.mkdir(ds.dir+'fits')
    fields=['density','temperature','velocity','magnetic_field']
    for field in fields:
        if field is 'temperature':
            data=coolftn.get_temp(ds.read_all_data('T1'))
        else:
            data=ds.read_all_data(field)
        fitsname='%s/fits/%s.%s.%s.fits' % (ds.dir,ds.id,ds.step,field[:4])
        hdr = fits.Header()
        hdr['field']=field
        hdr['time']=ds.domain['time']
        hdr['xmin']=(ds.domain['left_edge'][0],'pc')
        hdr['xmax']=(ds.domain['right_edge'][0],'pc')
        hdr['ymin']=(ds.domain['left_edge'][1],'pc')
        hdr['ymax']=(ds.domain['right_edge'][1],'pc')
        hdr['zmin']=(ds.domain['left_edge'][2],'pc')
        hdr['zmax']=(ds.domain['right_edge'][2],'pc')
        hdr['dx']=(ds.domain['dx'][0],'pc')
        hdr['dy']=(ds.domain['dx'][1],'pc')
        hdr['dz']=(ds.domain['dx'][2],'pc')
        hdr['unit']=(units[field].cgs.value,units[field].cgs.unit)
        hdu = fits.PrimaryHDU(data,header=hdr)
        hdu.writeto(fitsname,clobber=True)
