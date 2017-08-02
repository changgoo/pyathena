from ytathena import *
import yt
from matplotlib.colors import LogNorm,SymLogNorm,NoNorm,Normalize
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
from pyathena import read_starvtk
from slices.scatter_sp import scatter_sp,sp_legend
import pandas as pd
from shiftedColorMap import *


aux={}
aux['nH']={'title':r'$n_{\rm H} [{\rm cm}^{-3}]$',
  'cmap':'RdYlBu_r','cmin':2.e-5,'cmax':2.e2,
  'factor':1.0,'norm':LogNorm()}
aux['pok']={'title':r'$P/k_B [{\rm cm}^{-3}{\rm K}]$',
  'cmap':'gnuplot2','cmin':20,'cmax':2.e6,
  'factor':1.0,'norm':LogNorm()}
aux['temperature']={'title':r'$T [{\rm K}]$',
  'cmap':shiftedColorMap(plt.cm.Spectral_r,midpoint=2.5/6),
  'cmin':50,'cmax':5.e7,
  'factor':1.0,'norm':LogNorm()}
aux['surface_density']={'title':r'$\Sigma [{\rm M}_{\odot} {\rm pc}^{-2}]$',
  'cmap':'pink_r','cmin':0.5,'cmax':100,
  'factor':1.0,'norm':LogNorm()}
aux['velocity_z']={'title':r'$v_z [{\rm km/s}]$',
  'cmap':'RdBu_r','cmin':-200,'cmax':200,
  'factor':1.0,'norm':Normalize(),'unit':'km/s','log':False}
aux['velocity']={'title':r'$v [{\rm km/s}]$',
  'cmap':'jet','cmin':1,'cmax':1000,
  'factor':1.0,'norm':LogNorm(),'unit':'km/s'}
aux['magnetic_field_strength']={'title':r'$|B| [\mu{\rm G}]$',
  'cmap':'pink_r','cmin':0.05,'cmax':10,
  'factor':1.e6,'norm':LogNorm(),'unit':'gauss'}

fields=['nH','temperature','pok','velocity_z','magnetic_field_strength']
nf=len(fields)

def main(**kwargs):
 
    dir = kwargs['base_directory']+kwargs['directory']
    fname=dir+'id0/%s.%4.4d.vtk' % (kwargs['id'],kwargs['itime'])
    starfname=dir+'id0/%s.%4.4d.starpar.vtk' % (kwargs['id'],kwargs['itime'])

    ds = yt.load(fname,units_override=unit_base)
    ds.coordinates.x_axis[1]=0
    ds.coordinates.x_axis['y']=0
    ds.coordinates.y_axis[1]=2
    ds.coordinates.y_axis['y']=2

    sp=read_starvtk(starfname)
    ratio=ds.domain_width[2]/ds.domain_width[0]
    ix=2
    iz=ix*ratio
    fig=plt.figure(figsize=(ix*nf,iz+ix+0.1*ix))
    gs = gridspec.GridSpec(3,nf,height_ratios=[0.1*ix,iz,ix],
         hspace=0.01,wspace=0.02)

    c=ds.domain_center
    images=[]
    for i,axis in enumerate(['y','z']):
        slc= yt.SlicePlot(ds,axis,fields)
        slc_frb = slc.data_source.to_frb(slc.width[0],
                  slc.width,c,slc.width[1])
        extent=np.array(slc_frb.bounds)/1.e3
        for j,f in enumerate(fields):
            if aux[f].has_key('unit'):
                data = slc_frb[f].in_units(aux[f]['unit']).d
            else:
                data = slc_frb[f].d
            if aux[f].has_key('factor'):
                data *= aux[f]['factor']
            ax=plt.subplot(gs[i+1,j])
            im=ax.imshow(data,origin='lower',norm=aux[f]['norm'])
            im.set_clim((aux[f]['cmin'],aux[f]['cmax']))
            im.set_cmap(aux[f]['cmap'])
            im.set_extent(extent)
            images.append(im)
            if j == 0: 
		scatter_sp(sp,ax,axis=axis,runaway=False)
            elif j == 1: 
		scatter_sp(sp,ax,axis=axis)
            ax.set_xlim(extent[0],extent[1])
            ax.set_ylim(extent[2],extent[3])
         
    for j,(im,f) in enumerate(zip(images[:nf],fields)):
        cax=plt.subplot(gs[0,j])
        cbar = fig.colorbar(im,cax=cax,orientation='horizontal')
        cbar.set_label(aux[f]['title'])
        cax.xaxis.tick_top()
        cax.xaxis.set_label_position('top')
        if f is 'nH': cbar.set_ticks([1.e-4,1.e-2,1,1.e2])
        if f.startswith('velocity'): cbar.set_ticks([-100,0,100])
            
    axes=fig.axes
    sp_legend(axes[nf-1])
    plt.setp([ax.get_xticklabels() for ax in axes[:2*nf]],visible=False)
    plt.setp([ax.get_yticklabels() for ax in axes[:2*nf]],visible=False)
    for ax in fig.axes:
        ax.tick_params(axis='both',which='major',labelsize=12)
    if kwargs['label']:
        plt.setp(axes[nf:2*nf],'xlabel','x [kpc]')
        plt.setp(axes[0],'ylabel','z [kpc]')
        plt.setp(axes[nf],'ylabel','y [kpc]')
        plt.setp([ax.get_xticklabels() for ax in axes[nf:]], visible=True)
        plt.setp([ax.get_yticklabels() for ax in axes[:2*nf:nf]], visible=True)
        plt.setp([ax.xaxis.get_majorticklabels() for ax in axes[nf:2*nf]], rotation=45 )

    pngfname="%s/4slices/%s.%4.4d.4slices.png" % (kwargs['base_directory'],kwargs['id'],kwargs['itime'])
    fig.savefig(pngfname,dpi=300,bbox_inches='tight')

if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('-b','--base_directory',type=str,
                      default='/tigress/changgoo/',
                      help='base working directory')
  parser.add_argument('-d','--directory',type=str,default='',
                      help='working directory')
  parser.add_argument('-i','--id',type=str,
                      help='id of dataset')
  parser.add_argument('-t','--itime',type=int,default='',
                      help='time')
  parser.add_argument('-l','--label',action='store_true',
                      help='toggle label')
  args = parser.parse_args()

  main(**vars(args))
