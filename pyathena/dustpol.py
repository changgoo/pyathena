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
import copy

def data_read(fname,vel=True,mhd=True,usefits=True):
    dir, id, step, ext, mpi = pa.parse_filename(fname)
    if usefits:
        from astropy.io import fits
        fitsname='%s/fits/%s.%s.%s.fits' % (dir,id,step,'dens')
        
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
    else:
        ds=pa.AthenaDataSet(fname)
        rstfnames=glob.glob('%s/id0/%s.????.rst' % (ds.dir,ds.id)) \
                 +glob.glob('%s/rst/%s.????.rst' % (ds.dir,ds.id))
        par,blocks,fields=pa.parse_par(rstfnames[0])
 
        domain=ds.domain
        domain['qshear']=eval(par['problem']['qshear'][0])
        domain['Omega']=eval(par['problem']['Omega'][0])
        fields=['density']
        if vel:
            fields.append('velocity1')
            fields.append('velocity2')
            fields.append('velocity3')
        if mhd: 
            fields.append('magnetic_field1')
            fields.append('magnetic_field2')
            fields.append('magnetic_field3')

    print domain['time']

    data={}
    if usefits:
        for f in fields:
            fitsname='%s/fits/%s.%s.%s.fits' % (dir,id,step,f[:4])
            if f.startswith('velo') or f.startswith('magn'):
                hdulist=fits.open(fitsname,memmap=False,lazy_load_hdus=False)
                hdulist.readall()
                #fitsdata=fits.getdata(fitsname)
                for iaxis in range(3):
                    data['%s%d' % (f,iaxis+1)]=hdulist[0].data[:,:,:,iaxis]
                    #data['%s%d' % (f,iaxis+1)]=fitsdata[:,:,:,iaxis]
                hdulist.close()
            else:
                hdulist=fits.open(fitsname,memmap=False,lazy_load_hdus=False)
                hdulist.readall()
                #fitsdata=fits.getdata(fitsname)
                data[f]=hdulist[0].data
                #data[f]=fitsdata
    else:
        for f in fields:
            data[f] = ds.read_all_data(f)
        data['temperature']=coolftn.get_temp(ds.read_all_data('T1'))
        fields.append('temperature')
    
    if vel:
        r3d,x3d,y3d,z3d=pa.pos3d(domain)
        vy0=-domain['qshear']*domain['Omega']*x3d
        data['velocity2'] -= vy0

    print data.keys()

    return data,domain

def make_pol_map(data,domain,mapfname='nomap',Nside=4,center=[0.,0.,0.],write=False):
    deltas=domain['dx'][2]/2.
    rmax=domain['Lx'][0]*2
    print 'rmax', rmax
    print 'ds  ', deltas

    npix=hp.nside2npix(Nside)

    Imap=[]
    Qmap=[]
    Umap=[]
    for ipix in xrange(npix):
        angle = np.rad2deg(hp.pix2ang(Nside,ipix))
        if ipix % Nside == 0: print '%d of %d' % (ipix,npix)
        if np.abs(90-angle[0]) > 0:
            los=syn.get_los(data,domain,Nside,ipix,vel=False,mhd=True,\
                            rmax=rmax,deltas=deltas,center=center)
            if write:
                pd.DataFrame(los).to_pickle('%s%s.i%d.p' % (losdir,step,ipix))
                
            I,Q,U=syn.los_to_dustpol(los,deltas=deltas)
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
    if mapfname != 'nomap': 
        hp.write_map(mapfname,[Imap,Qmap,Umap],column_units=['MJy/sr','MJy/sr','MJy/sr'])

def make_pol_map_all(data,domain,mapfname='nomap',dl=5.,\
  Nside=4,center=[0.,0.,0.],write=False):
    deltas=domain['dx'][2]/2.
    rmax=1024.#domain['Lx'][2]/2
    print 'rmax', rmax
    print 'ds  ', deltas

    npix=hp.nside2npix(Nside)

    newfile=True
    nslice=8
    dr=int(rmax/nslice)
    nr=int(rmax/deltas)
    print dr, nr
    l0=0.0
    for ipix in xrange(npix):
        if newfile:
            Imap=[]
            Qmap=[]
            Umap=[]
            for map in [Imap,Qmap,Umap]: 
                for r in range(nslice): map.append({'ipix':[],'cnm':[],'unm':[],'wnm':[]})
        angle = np.rad2deg(hp.pix2ang(Nside,ipix))
        los=syn.get_los(data,domain,Nside,ipix,vel=False,mhd=True,\
                        rmax=rmax,deltas=deltas,center=center)
        for f in los.keys():
            los[f]=los[f][:nr].reshape(nslice,dr)
        idx={}
        idx['cnm']=los['temperature'] < 184.
        idx['unm']=(los['temperature'] >= 184.)&(los['temperature'] < 5050.)
        idx['wnm']=(los['temperature'] >= 5050.)&(los['temperature'] < 2.e4)
        los_p = {'cnm':copy.deepcopy(los),'unm':copy.deepcopy(los),'wnm':copy.deepcopy(los)}
        for p in ['cnm','unm','wnm']:
            los_p[p]['density'][~idx[p]] = 0.0
        for r in range(nslice):
            Imap[r]['ipix'].append(ipix)
            Qmap[r]['ipix'].append(ipix)
            Umap[r]['ipix'].append(ipix)
            for p in ['cnm','unm','wnm']:
                newlos={}
                for f in los.keys(): newlos[f]=los_p[p][f][r,:]
                I,Q,U=syn.los_to_dustpol(newlos,deltas=deltas)
                Imap[r][p].append(I.value)
                Qmap[r][p].append(Q.value)
                Umap[r][p].append(U.value)
        if angle[0] > l0+dl or ipix == (npix-1):
            for r in range(nslice):
               df=pd.DataFrame({'I':Imap[r],'Q':Qmap[r],'U':Umap[r]})
               df.to_pickle(mapfname+'.l%03d.r%d.p' % (int(l0),r))
            l0=angle[0]
            newfile=True
            print l0,ipix,len(Imap[0]['cnm'])
        else:
            newfile=False
## global parameters
base='/tigress/changgoo/'
id='MHD_2pc_S'
Nside=128
center=[[0.,0.,0.],[256.,256.,0.],[-256.,256.,0.],[-256.,-256.,0.],[256.,-256.,0.]]
ifile=12
##

fnames=glob.glob('%s%s/id0/%s.????.vtk' % (base,id,id))
fnames.sort()

if __name__ == '__main__':
    #def runall():
    fname=fnames[ifile]
    dir, id, step, ext, mpi = pa.parse_filename(fname)
    data,domain=data_read(fname,mhd=True,vel=False,usefits=False)
    for i,c in enumerate(center):
#        losdir='%s/dustpol/%s_%d_c%d/' % (dir,step,Nside,i)
#        mapfname='%s/dustpol/%s.%s.n%d.c%d.fits' % (dir,id,step,Nside,i)
#        if not os.path.isdir(losdir): os.mkdir(losdir)
#        make_pol_map(data,domain,mapfname=mapfname,center=c,Nside=Nside,write=False)

        mapfname='%s/dustpol/%s.%s.n%d.c%d' % (dir,id,step,Nside,i)
        print mapfname,i,c
        make_pol_map_all(data,domain,mapfname=mapfname,center=c,Nside=Nside,write=False,dl=5)

