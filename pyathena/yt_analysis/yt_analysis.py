import ytathena as ya
import yt
import glob
import argparse
import os
import cPickle as pickle

import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,SymLogNorm,NoNorm,Normalize
from pyathena import read_starvtk,texteffect,set_units
from plot_slices import slice2 as plot_slice
from plot_projection import plot_projection
import numpy as np
import string


fields=['cell_volume','cell_mass']
unit=set_units(muH=1.4271)
Myr=unit['time'].to('Myr').value

def compare_files(source, output):
    smtime=os.path.getmtime(source)
    if os.path.isfile(output):
        omtime=os.path.getmtime(output)
        if omtime < smtime:
            return False
        else:
            return True
    else:
        return False

class my_pdf(object):
    def __init__(self,box):
        self.time=box.ds.current_time.in_units('Myr').v
        self.left_edge=np.array(box.left_edge)
        self.right_edge=np.array(box.right_edge)
        self.pdf={}
        return

    def add_pdf(self,pdf,bin_field):
        key=string.join(bin_field,'-')
        #print key
        self.pdf[key]={}
        self.pdf[key]['xbin']=np.array(pdf.x_bins)
        self.pdf[key]['ybin']=np.array(pdf.y_bins)
        self.pdf[key]['cell_mass']=np.array(pdf['cell_mass'])
        self.pdf[key]['cell_volume']=np.array(pdf['cell_volume'])
        return

def draw_pdf(ax,pdf,field='cell_mass',key='nH-pok'):
    global aux
    bins=string.split(key,'-')
    xbin=pdf[key]['xbin']
    ybin=pdf[key]['ybin']
    data=pdf[key][field]
    im=ax.pcolormesh(xbin,ybin,data.T,norm=LogNorm())
    if aux[bins[0]]['log']: ax.set_xscale('log')
    if aux[bins[1]]['log']: ax.set_yscale('log')
    ax.set_xlabel(aux[bins[0]]['label'])
    ax.set_ylabel(aux[bins[1]]['label'])
    ax.set_xlim(aux[bins[0]]['limits'])
    ax.set_ylim(aux[bins[1]]['limits'])
    return im


def projection(ds,surfname):
    time = ds.current_time.in_units('Myr').v
    ds.coordinates.x_axis[1]=0
    ds.coordinates.x_axis['y']=0
    ds.coordinates.y_axis[1]=2
    ds.coordinates.y_axis['y']=2

    c=ds.domain_center
    dx=ds.domain_width/ds.domain_dimensions

    surf_data={}
    surf_data['time']=time
    for i,axis in enumerate(['x','y','z']):
        proj = ds.proj('density',axis=axis)
        ix=ds.coordinates.x_axis[axis]
        iy=ds.coordinates.y_axis[axis]
        res=(ds.domain_dimensions[ix],ds.domain_dimensions[iy])
        frb = proj.to_frb(ds.domain_width[ix],res,c,ds.domain_width[iy])
        surf=np.array(frb['density'].in_units('Msun/pc**2'))
        bounds = np.array(frb.bounds)
        surf_data[axis]={'data':surf,'bounds':bounds}
    if yt.is_root():
        pickle.dump(surf_data,open(surfname,'wb'),pickle.HIGHEST_PROTOCOL)

    scal_fields=ya.get_scalars(ds)
    for nscal,sf in enumerate(scal_fields):
        scal_data={}
        scal_data['time']=time
        for i,axis in enumerate(['x','y','z']):
            proj = ds.proj(sf,axis=axis,weight_field='density')
            ix=ds.coordinates.x_axis[axis]
            iy=ds.coordinates.y_axis[axis]
            res=(ds.domain_dimensions[ix],ds.domain_dimensions[iy])
            frb = proj.to_frb(ds.domain_width[ix],res,c,ds.domain_width[iy])
            scal=np.array(frb[sf])
            bounds = np.array(frb.bounds)
            scal_data[axis]={'data':scal,'bounds':bounds}
        if yt.is_root():
            scalfname=surfname.replace('surf.p','scal{}.p'.format(nscal))
            pickle.dump(scal_data,open(scalfname,'wb'),pickle.HIGHEST_PROTOCOL)

def phase(sq,phfname,bin_fields):
    global aux
    pdfs=my_pdf(sq)

    for bf in bin_fields:
        n_bins=(aux[bf[0]]['n_bins'],aux[bf[1]]['n_bins'])
        logs={}
        unit={}
        extrema={}
        for b in bf: 
            logs[b]=aux[b]['log'] 
            if aux[b].has_key('unit'): unit[b]=aux[b]['unit'] 
            if aux[b].has_key('limits'): extrema[b]=aux[b]['limits']
        pdf=yt.create_profile(sq,bf,fields=fields,
              n_bins=n_bins,logs=logs,extrema=extrema,units=unit,
              weight_field=None,fractional=True)
        pdfs.add_pdf(pdf,bf)
    if yt.is_root():
        pickle.dump(pdfs,open(phfname,'wb'),pickle.HIGHEST_PROTOCOL)

def slices(ds,slcfname,slc_fields):
    global aux
    ds.coordinates.x_axis[1]=0
    ds.coordinates.x_axis['y']=0
    ds.coordinates.y_axis[1]=2
    ds.coordinates.y_axis['y']=2

    time = ds.current_time.in_units('Myr').v
    c=ds.domain_center
    dx=ds.domain_width/ds.domain_dimensions

    slc_data={}
    slc_data['time']=time
    ya.check_aux(slc_fields)

    for i,axis in enumerate(['x','y','z']):
        slc=yt.SlicePlot(ds,axis,slc_fields)
        ix=ds.coordinates.x_axis[axis]
        iy=ds.coordinates.y_axis[axis]
#        res=(int(slc.width[1]/dx[ix]),int(slc.width[0]/dx[iy]))
        res=(ds.domain_dimensions[ix],ds.domain_dimensions[iy])
        if(axis == 'y'): print res
        slc_frb = slc.data_source.to_frb(slc.width[0],
                  res,c,slc.width[1])
        extent=np.array(slc_frb.bounds)/1.e3
        slc_data[axis]={} 
        slc_data[axis+'extent']=extent
        for f in slc_fields:
            if aux[f].has_key('unit'):
                slc_data[axis][f] = np.array(slc_frb[f].in_units(aux[f]['unit']))
            else:
                slc_data[axis][f] = np.array(slc_frb[f])
            if aux[f].has_key('factor'): slc_data[axis][f] *= aux[f]['factor']

    if yt.is_root():
        pickle.dump(slc_data,open(slcfname,'wb'),pickle.HIGHEST_PROTOCOL)


def plot_phase(phfname,bin_fields):
    global aux
    plt.rc('font',size=10)
    plt.rc('xtick',labelsize=10)
    plt.rc('ytick',labelsize=10)


    nbf=len(bin_fields)
    nrow=4
    ncol=nbf/nrow+1
    fig=plt.figure(2,figsize=(6*ncol,10))
#    gs = gridspec.GridSpec(nrow,3*ncol,width_ratios=[1,1,0.1]*ncol,
#            wspace=0.5,hspace=0.4)
    
    pdfs=pickle.load(open(phfname,'rb'))
    for j,bf in enumerate(bin_fields):
      for i in [0,1]:
        ax=fig.add_subplot(nrow,ncol*2,j*2+i+1)
        #ax=plt.subplot(gs[j%nrow,(j/nrow)*3+i])
        im = draw_pdf(ax,pdfs.pdf,field=fields[i],
                      key=string.join(bf,'-'))
        im.set_cmap(plt.cm.cubehelix_r)
        im.set_clim(1.e-7,1.e-1)
        if j<2:
            if i==0: ax.set_title('volume-PDF')
            if i==1: ax.set_title('mass-PDF')

    fig.tight_layout()

    cax=fig.add_axes([1.0, 0.05, 0.02, 0.90])
    fig.colorbar(im,cax=cax,orientation='vertical')
    
    pngfname=phfname[:-1]+'png'
    #canvas = mpl.backends.backend_agg.FigureCanvasAgg(fig)
    #canvas.print_figure(pngfname,bbox_inches='tight',num=2,dpi=150)
    plt.savefig(pngfname,bbox_inches='tight',num=2,dpi=150)
    plt.close()

def main(**kwargs):
    global aux
    #yt.funcs.mylog.setLevel(50) 
    dir = kwargs['base_directory']+kwargs['directory']
    fname=glob.glob(dir+'id0/'+kwargs['id']+'.????.vtk')
    fname.sort()
    aux=ya.set_aux(kwargs['id'])

    if yt.is_root(): print fname[0]

    if kwargs['range'] != '':
        sp=kwargs['range'].split(',')
        start = eval(sp[0])
        end = eval(sp[1])
        fskip = eval(sp[2])
    else:
        start = 0
        end = len(fname)
        fskip = 1
    fname=fname[start:end:fskip]

    ngrids=len(glob.glob(dir+'id*/'+kwargs['id']+'*'+fname[0][-8:]))
    if kwargs['parallel']:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        nprocs = comm.size
        rank = comm.rank
    else:
        nprocs = 1
        rank = 0
    print ngrids,rank,nprocs,yt.is_root()

    do_phase=kwargs['phase']

    if ngrids > nprocs: ds = yt.load(fname[0],units_override=ya.unit_base)
    else: ds = yt.load(fname[0],units_override=ya.unit_base, nprocs=nprocs)

    mhd=('athena','cell_centered_B_x') in ds.field_list
    cooling=('athena','pressure') in ds.field_list
    rotation=kwargs['rotation'] != 0.
    if rotation: 
      ya.Omega=ya.YTQuantity(kwargs['rotation'],'km/s/kpc')
      if kwargs['rotation']== 280: 
        aux=ya.set_aux('starburst')
    if rank == 0:
      print "phase plot:", do_phase
      print "MHD:", mhd
      print "cooling:", cooling
      print "rotation:", rotation, ya.Omega

    bin_fields=[]
    bin_fields.append(['nH','pok'])
    bin_fields.append(['nH','temperature'])
    bin_fields.append(['velocity_z','temperature'])
    if rotation: 
        bin_fields.append(['dvelocity_magnitude','temperature'])
    else:
        bin_fields.append(['velocity_magnitude','temperature'])
    if mhd:
        bin_fields.append(['nH','mag_pok'])
        bin_fields.append(['nH','ram_pok_z'])
        bin_fields.append(['nH','plasma_beta'])

    slc_fields=['nH','pok','temperature','velocity_z','ram_pok_z']
    fields_to_draw=['star_particles','nH','temperature','pok','velocity_z']
    if mhd:
      slc_fields.append('magnetic_field_strength')
      slc_fields.append('mag_pok')
      fields_to_draw.append('magnetic_field_strength')
    scal_fields=ya.get_scalars(ds)
    slc_fields+=scal_fields

    if rank == 0:
        if not os.path.isdir(dir+'slice/'): os.mkdir(dir+'slice/')
        if not os.path.isdir(dir+'surf/'): os.mkdir(dir+'surf/')
        if not os.path.isdir(dir+'phase/'): os.mkdir(dir+'phase/')
    for i,f in enumerate(fname):
        slcfname=dir+'slice/'+kwargs['id']+f[-9:-4]+'.slice.p'
        surfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.surf.p'
        if do_phase:
            phfname=dir+'phase/'+kwargs['id']+f[-9:-4]+'.phase.p'
        else:
            phfname=f
        if rank == 0: print i,compare_files(f,slcfname),compare_files(f,surfname),compare_files(f,phfname)
        if compare_files(f,surfname) and \
           compare_files(f,phfname) and \
           compare_files(f,slcfname):
            if rank == 0: print 'all data is already there'
        else:
            if ngrids > nprocs: ds = yt.load(f,units_override=ya.unit_base)
            else: ds = yt.load(f,units_override=ya.unit_base, nprocs=nprocs)
            ya.add_yt_fields(ds,mhd=mhd,rotation=rotation,cooling=cooling)
            if not compare_files(f,surfname):
                if rank == 0: print 'projecting...'
                projection(ds,surfname)
            if not compare_files(f,slcfname):
                if rank == 0: print 'slicing...'
                slices(ds,slcfname,slc_fields)
            if not compare_files(f,phfname) and do_phase:
                if rank == 0: print 'binning...'
                le=np.array(ds.domain_left_edge)
                re=np.array(ds.domain_right_edge)
                sq=ds.box(le,re)

                phase(sq,phfname,bin_fields)

    for i,f in enumerate(fname):
        slcfname=dir+'slice/'+kwargs['id']+f[-9:-4]+'.slice.p'
        surfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.surf.p'
        starpardir='id0/'
        if os.path.isdir(dir+'starpar/'): starpardir='starpar/'
        starfname=dir+starpardir+kwargs['id']+f[-9:-4]+'.starpar.vtk'

        if do_phase:
            phfname=dir+'phase/'+kwargs['id']+f[-9:-4]+'.phase.p'
        else:
            phfname=f

        if i%nprocs == rank:
          if not compare_files(surfname,surfname+'ng'):
            print 'drawing for %s on %d' % (surfname,rank)
            plot_projection(surfname,starfname,runaway=False)
          if not compare_files(slcfname,slcfname+'ng'):
            print 'drawing for %s on %d' % (slcfname,rank)
            plot_slice(slcfname,starfname,fields_to_draw)
          if not compare_files(phfname,phfname+'ng') and do_phase:
            if do_phase:
                print 'drawing for %s on %d' % (phfname,rank)
                plot_phase(phfname,bin_fields)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
 
    parser.add_argument('-b','--base_directory',type=str,
                        default='/tigress/changgoo/',
                        help='base working directory')
    parser.add_argument('-d','--directory',type=str,default='',
                        help='working directory')
    parser.add_argument('-i','--id',type=str,
                        help='id of dataset')
    parser.add_argument('-p','--parallel',action='store_true',
                        help='parallel')
    parser.add_argument('-ph','--phase',action='store_true',
                        help='phase diagram analysis')
    parser.add_argument('-ro','--rotation',type=float,default=28.,
                        help='rotational velocity')
    parser.add_argument('-r','--range',type=str,default='',
                       help='time range, start:end:skip')
    args = parser.parse_args()
 
    if vars(args)['parallel']: yt.enable_parallelism()
    main(**vars(args))
