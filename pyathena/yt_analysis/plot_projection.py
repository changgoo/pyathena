
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

    if 'time' in frb:
        tMyr=frb['time']
    else:
        time,sp=read_starvtk(starfname,time_out=True)
        tMyr=time*Myr

    if 'z' in frb: frb=frb['z']
    extent=np.array(frb['bounds'])/1.e3
    x0=extent[0]
    y0=extent[2]
    Lx=extent[1]-extent[0]
    Lz=extent[3]-extent[2]
 
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

def plot_projection_Z(surfname,starfname,axis='y',stars=True,writefile=True,runaway=True):
    aux=ya.set_aux(os.path.basename(surfname))

    plt.rc('font',size=16)
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16)


    if stars: sp=read_starvtk(starfname)
    frb=pickle.load(open(surfname,'rb'))#,encoding='latin1')
    if 'time' in frb:
        tMyr=frb['time']
    else:
        time,sp=read_starvtk(starfname,time_out=True)
        tMyr=time*Myr

    extent=np.array(frb[axis]['bounds'])/1.e3
    x0=extent[0]
    y0=extent[2]
    Lx=extent[1]-extent[0]
    Lz=extent[3]-extent[2]
 
    ix=3
    iz=ix*Lz/Lx
    fig=plt.figure(0,figsize=(ix*2+0.5,iz))
    gs = gridspec.GridSpec(2,3,width_ratios=[1,1,0.03],wspace=0.0)
    ax1=plt.subplot(gs[:,0])
    im1=ax1.imshow(frb['y']['data'],norm=LogNorm(),origin='lower')
    im1.set_extent(extent)
    im1.set_cmap(aux['surface_density']['cmap'])
    clim=aux['surface_density']['clim']
    cmin=clim[0]
    cmax=clim[1]
    clim=(cmin*0.1,cmax)
    im1.set_clim(clim)
    ax1.text(extent[0]*0.9,extent[3]*0.9,
            't=%3d Myr' % tMyr,ha='left',va='top',**(texteffect()))
    if stars: scatter_sp(sp,ax1,axis='y',runaway=runaway,type='surf')

    ax2=plt.subplot(gs[:,1])
    im2=ax2.imshow(frb['x']['data'],norm=LogNorm(),origin='lower')
    im2.set_extent(extent)
    im2.set_cmap(aux['surface_density']['cmap'])
    im2.set_clim(clim)
    if stars: scatter_sp(sp,ax2,axis='x',runaway=runaway,type='surf')


    cax=plt.subplot(gs[0,2])
    cbar = fig.colorbar(im1,cax=cax,orientation='vertical')
    cbar.set_label(aux['surface_density']['label'])

    if stars:
      cax=plt.subplot(gs[1,2])
      cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
             cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
             orientation='vertical')
      cbar.set_label(r'${\rm age [Myr]}$')
 
      norm_factor=2.
      s1=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e3)/norm_factor,color='k',
        alpha=.8,label=r'$10^3 M_\odot$')
      s2=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e4)/norm_factor,color='k',
        alpha=.8,label=r'$10^4 M_\odot$')
      s3=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e5)/norm_factor,
        color='k',alpha=.8,label=r'$10^5 M_\odot$')

      ax1.set_xlim(x0,x0+Lx)
      ax1.set_ylim(y0,y0+Lz);
      ax2.set_xlim(x0,x0+Lx)
      ax2.set_ylim(y0,y0+Lz);
      legend=ax1.legend((s1,s2,s3),(r'$10^3 M_\odot$',r'$10^4 M_\odot$',r'$10^5 M_\odot$'),
                        loc='lower left',fontsize='medium',frameon=True)

    ax1.set_xlabel('x [kpc]')
    ax1.set_ylabel('z [kpc]')

    ax2.set_xlabel('y [kpc]')
    plt.setp(ax2.get_yticklabels(),visible=False)
    pngfname=surfname+'ng'
    if writefile:
        plt.savefig(pngfname,bbox_inches='tight',num=0,dpi=150)
        plt.close()
    else:
        return fig

def plot_projection_icm(surfname,scalfname,starfname,
  stars=True,writefile=True,runaway=True):
    aux=ya.set_aux(os.path.basename(surfname))

    plt.rc('font',size=16)
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16)


    if stars: sp=read_starvtk(starfname)
    frb=pickle.load(open(surfname,'rb'))#,encoding='latin1')
    icm=pickle.load(open(scalfname,'rb'))#,encoding='latin1')
    if 'time' in frb:
        tMyr=frb['time']
    else:
        time,sp=read_starvtk(starfname,time_out=True)
        tMyr=time*Myr

    extent=np.array(frb['y']['bounds'])/1.e3
    x0=extent[0]
    y0=extent[2]
    Lx=extent[1]-extent[0]
    Lz=extent[3]-extent[2]
 
    ix=3
    iz=ix*Lz/Lx

    cism=plt.cm.bone_r
    cicm=plt.cm.Reds

    fig=plt.figure(0,figsize=(ix*2+0.5,iz))
    gs = gridspec.GridSpec(3,3,width_ratios=[1,1,0.03],wspace=0.0)
    ax1=plt.subplot(gs[:,0])
    im1=ax1.imshow(frb['y']['data'],norm=LogNorm(),origin='lower')
    im1.set_extent(extent)
    im1.set_cmap(cism)
    clim=aux['surface_density']['clim']
    cmin=clim[0]
    cmax=clim[1]
    clim=(cmin*0.1,cmax)
    im1.set_clim(clim)

    im11=ax1.imshow(icm['y']['data'],norm=LogNorm(),origin='lower',alpha=0.7)
    im11.set_extent(extent)
    im11.set_cmap(cicm)
    im11.set_clim((1.e-2,1))

    ax1.text(extent[0]*0.9,extent[3]*0.9,
            't=%3d Myr' % tMyr,ha='left',va='top',**(texteffect(20)))
    if stars: scatter_sp(sp,ax1,axis='y',runaway=runaway,type='surf')

    ax2=plt.subplot(gs[:,1])
    im2=ax2.imshow(frb['x']['data'],norm=LogNorm(),origin='lower')
    im2.set_extent(extent)
    im2.set_cmap(cism)
    im2.set_clim(clim)

    im21=ax2.imshow(icm['x']['data'],norm=LogNorm(),origin='lower',alpha=0.7)
    im21.set_extent(extent)
    im21.set_cmap(cicm)
    im21.set_clim((1.e-2,1))


    if stars: scatter_sp(sp,ax2,axis='x',runaway=runaway,type='surf')


    cax=plt.subplot(gs[0,2])
    cbar = fig.colorbar(im1,cax=cax,orientation='vertical')
    cbar.set_label(aux['surface_density']['label'])

    cax=plt.subplot(gs[1,2])
    cbar = fig.colorbar(im11,cax=cax,orientation='vertical')
    cbar.set_label(r'$C_{\rm ICM}$')


    if stars:
      cax=plt.subplot(gs[2,2])
      cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
             cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
             orientation='vertical')
      cbar.set_label(r'${\rm age [Myr]}$')
 
      norm_factor=2.
      s1=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e3)/norm_factor,color='k',
        alpha=.8,label=r'$10^3 M_\odot$')
      s2=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e4)/norm_factor,color='k',
        alpha=.8,label=r'$10^4 M_\odot$')
      s3=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e5)/norm_factor,
        color='k',alpha=.8,label=r'$10^5 M_\odot$')

      ax1.set_xlim(x0,x0+Lx)
      ax1.set_ylim(y0,y0+Lz);
      ax2.set_xlim(x0,x0+Lx)
      ax2.set_ylim(y0,y0+Lz);
      legend=ax1.legend((s1,s2,s3),(r'$10^3 M_\odot$',r'$10^4 M_\odot$',r'$10^5 M_\odot$'),
                        loc='lower left',fontsize='medium',frameon=True)

    ax1.set_xlabel('x [kpc]')
    ax1.set_ylabel('z [kpc]')

    ax2.set_xlabel('y [kpc]')
    plt.setp(ax2.get_yticklabels(),visible=False)
    pngfname=surfname+'ng'
    if writefile:
        plt.savefig(pngfname,bbox_inches='tight',num=0,dpi=150)
        plt.close()
    else:
        return fig
