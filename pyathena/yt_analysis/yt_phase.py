import matplotlib as mpl
mpl.use('agg')

from ytathena import *
import glob
import argparse
import os

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,SymLogNorm,NoNorm,Normalize
import numpy as np

import yt
from shiftedColorMap import *
import string
import cPickle as pickle


aux={}
aux['nH']=dict(label=r'$n_H [{\rm cm}^{-3}]$',
               unit='cm**(-3)', limits=(1.e-6,1.e4), 
               n_bins=128, log=True)
aux['pok']=dict(label=r'$P/k_B [{\rm K} {\rm cm}^{-3}]$',
               unit='K*cm**(-3)', limits=(1.e-2,1.e8), 
               n_bins=128, log=True)
aux['temperature']=dict(label=r'$T [{\rm K}]$',
               unit='K', limits=(1.e0,1.e9), 
               n_bins=128, log=True)
aux['dvelocity_magnitude']=dict(label=r'$v [{\rm km/s}]$',
               unit='km/s', limits=(0.1,1.e4), 
               n_bins=128, log=True)
aux['velocity_z']=dict(label=r'$v_z [{\rm km/s}]$',
               unit='km/s', limits=(-1500,1500), 
               n_bins=256, log=False)
aux['magnetic_field_strength']=dict(label=r'$B [{\rm G}]$',
               unit='gauss',
               n_bins=128, log=True)
aux['mag_pok']=dict(label=r'$P_{\rm mag}/k_B [{\rm K}{\rm cm}^{-3}]$',
               unit='K*cm**(-3)', limits=(1.e-2,1.e8), 
               n_bins=128, log=True)
aux['ram_pok_z']=dict(label=r'$P_{\rm turb}/k_B [{\rm K}{\rm cm}^{-3}]$',
               unit='K*cm**(-3)', limits=(1.e-2,1.e8), 
               n_bins=128, log=True)
aux['plasma_beta']=dict(label=r'$\beta$', limits=(1.e-4,1.e16),
               n_bins=256, log=True)

bin_fields=[['nH','pok'],
            ['nH','temperature'],
            ['dvelocity_magnitude','temperature'],
            ['velocity_z','temperature'],
            ['nH','mag_pok'],
            ['nH','ram_pok_z'],
            ['nH','plasma_beta'],
           ]

fields=['cell_volume','cell_mass']
class my_pdf(object):
    def __init__(self,box):
        self.time=box.ds.current_time.in_units('Myr')
        self.left_edge=np.array(box.left_edge)
        self.right_edge=np.array(box.right_edge)
        self.pdf={}
        return

    def add_pdf(self,pdf,bin_fields):
        key=string.join(bin_fields,'-')
        print key
        self.pdf[key]={}
        self.pdf[key]['xbin']=np.array(pdf.x_bins)
        self.pdf[key]['ybin']=np.array(pdf.y_bins)
        self.pdf[key]['cell_mass']=np.array(pdf['cell_mass'])
        self.pdf[key]['cell_volume']=np.array(pdf['cell_volume'])
        return

def draw_pdf(ax,pdf,field='cell_mass',key='nH-pok'):
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

def main(**kwargs):

    dir = kwargs['base_directory']+kwargs['directory']
    fname=glob.glob(dir+'id0/'+kwargs['id']+'.????.vtk')
    fname.sort()

    ngrids=len(glob.glob(dir+'id*/'+kwargs['id']+fname[-9:]))
    comm = yt.communication_system.communicators[-1]
    nprocs = comm.size
    print ngrids,nprocs

    if yt.is_root():
        if not os.path.isdir(dir+'phase/'): os.mkdir(dir+'phase/')

    for f in fname:
        phfname=dir+'phase/'+kwargs['id']+f[-9:-4]+'.phase.p'
        if os.path.isfile(phfname):
            print '%s is already there' % phfname
        else:
            if ngrids > nprocs: ds = yt.load(f,units_override=unit_base)
            else: ds = yt.load(f,units_override=unit_base, nprocs=nprocs*8)
            le=np.array(ds.domain_left_edge)
            re=np.array(ds.domain_right_edge)
            sq=ds.box(le,re)
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
        
def plot(**kwargs):
    plt.rc('font',size=10)
    plt.rc('xtick',labelsize=8)
    plt.rc('ytick',labelsize=8)

    dir = kwargs['base_directory']+kwargs['directory']
    fname=glob.glob(dir+'id0/'+kwargs['id']+'.????.vtk')
    fname.sort()

    fig=plt.figure(0,figsize=(6,15))
    nbf=len(bin_fields)
    gs = gridspec.GridSpec(nbf,3,width_ratios=[1,1,0.03],
            wspace=0.1,hspace=0.4)
    ax=fig.add_subplot(111)
    
    for f in fname:
        phfname=dir+'phase/'+kwargs['id']+f[-9:-4]+'.phase.p'
        pdfs=pickle.load(open(phfname,'rb'))
        for j,bf in enumerate(bin_fields):
          for i in range(2):
            ax=plt.subplot(gs[j,i])
            im = draw_pdf(ax,pdfs.pdf,field=fields[i],
                          key=string.join(bf,'-'))
            im.set_cmap(plt.cm.cubehelix_r)
            im.set_clim(1.e-7,1.e-1)
            if i==1: 
                plt.setp(ax.get_yticklabels(), visible=False)
                plt.setp(ax,'ylabel','')
          cax=plt.subplot(gs[j,2])
          fig.colorbar(im,cax=cax,orientation='vertical')
        
#        fig.tight_layout()
        pngfname=dir+'phase/'+kwargs['id']+f[-9:-4]+'.phase.png'
        canvas = mpl.backends.backend_agg.FigureCanvasAgg(fig)
        canvas.print_figure(pngfname,bbox_inches='tight',num=0,dpi=150)
        fig.clf()

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
    args = parser.parse_args()
 
    if vars(args)['parallel']: yt.enable_parallelism()
    main(**vars(args))
    if yt.is_root(): plot(**vars(args))
