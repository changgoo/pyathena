
import pyathena.yt_analysis.ytathena as ya
import yt
import glob
import argparse
import os
import pickle as pickle

import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm,SymLogNorm,NoNorm,Normalize
from pyathena import read_starvtk,texteffect,set_units
import numpy as np
import string
from .scatter_sp import scatter_sp

def plot_projection(surfname,starfname,stars=True,writefile=True,runaway=True):
    aux=ya.set_aux(os.path.basename(surfname))

    plt.rc('font',size=11)
    plt.rc('xtick',labelsize=11)
    plt.rc('ytick',labelsize=11)

    fig=plt.figure(0,figsize=(5.5,5))
    gs = gridspec.GridSpec(2,2,width_ratios=[1,0.03],wspace=0.0)

    if stars: sp=read_starvtk(starfname)
    frb=pickle.load(open(surfname,'rb'))#,encoding='latin1')
    extent=np.array(frb['bounds'])/1.e3
    x0=extent[0]
    y0=extent[2]
    Lx=extent[1]-extent[0]
    Lz=extent[3]-extent[2]
 
    if 'time' in frb:
        tMyr=frb['time']
    else:
        time,sp=read_starvtk(starfname,time_out=True)
        tMyr=time*Myr
    ax=plt.subplot(gs[:,0])
    im=ax.imshow(frb['data'],norm=LogNorm(),origin='lower')
    im.set_extent(extent)
    im.set_cmap(aux['surface_density']['cmap'])
    im.set_clim(aux['surface_density']['clim'])
    ax.text(extent[0]*0.9,extent[3]*0.9,
            't=%3d Myr' % tMyr,ha='left',va='top',**(texteffect()))

    if stars: scatter_sp(sp,ax,axis='z',runaway=runaway,type='surf')

    cax=plt.subplot(gs[0,1])
    cbar = fig.colorbar(im,cax=cax,orientation='vertical')
    cbar.set_label(aux['surface_density']['label'])

    if stars:
      cax=plt.subplot(gs[1,1])
      cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
             cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
             orientation='vertical')
      cbar.set_label(r'${\rm age [Myr]}$')
 
      norm_factor=2.
      s1=ax.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e3)/norm_factor,color='k',
        alpha=.8,label=r'$10^3 M_\odot$')
      s2=ax.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e4)/norm_factor,color='k',
        alpha=.8,label=r'$10^4 M_\odot$')
      s3=ax.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e5)/norm_factor,
        color='k',alpha=.8,label=r'$10^5 M_\odot$')

      ax.set_xlim(x0,x0+Lx)
      ax.set_ylim(y0,y0+Lz);
      legend=ax.legend((s1,s2,s3),(r'$10^3 M_\odot$',r'$10^4 M_\odot$',r'$10^5 M_\odot$'), 
                        loc=2,ncol=3,bbox_to_anchor=(0.0, 1.15),
                        fontsize='medium',frameon=True)

    ax.set_xlabel('x [kpc]')
    ax.set_ylabel('y [kpc]')

    pngfname=surfname+'ng'
    if writefile:
        plt.savefig(pngfname,bbox_inches='tight',num=0,dpi=150)
        plt.close()
    else:
        return fig
