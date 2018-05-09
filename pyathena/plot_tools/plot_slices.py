
import glob
import os

import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm,SymLogNorm,NoNorm,Normalize
from pyathena import read_starvtk,texteffect,set_units
import numpy as np
import string
from .scatter_sp import scatter_sp
import cPickle as pickle

unit=set_units(muH=1.4271)
to_Myr=unit['time'].to('Myr').value

def slice(slcfname,starfname,fields_to_draw,zoom=1.,aux={},\
               writefile=True,tstamp=True,stars=True,field_label=True,norm_factor=2):
    plt.rc('font',size=14)
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)

    slc_data=pickle.load(open(slcfname,'rb'))
    x0=slc_data['yextent'][0]
    y0=slc_data['yextent'][3]
    Lx=slc_data['yextent'][1]-slc_data['yextent'][0]
    Lz=slc_data['yextent'][3]-slc_data['yextent'][2]
    Lz=Lz/zoom
    ix=2
    iz=ix*Lz/Lx
    nf=len(fields_to_draw)
    fig=plt.figure(1,figsize=(ix*nf,iz+ix*2))
    gs = gridspec.GridSpec(2,nf,height_ratios=[iz,ix])
    gs.update(left=0.10,right=0.90,wspace=0,hspace=0)

    sp=read_starvtk(starfname)
    if 'time' in slc_data:
        tMyr=slc_data['time']
    else:
        time,sp=read_starvtk(starfname,time_out=True)
        tMyr=time*Myr
    images=[]
    for i,axis in enumerate(['y','z']):
        for j,f in enumerate(fields_to_draw):
            data=slc_data[axis][f]
            ax=plt.subplot(gs[i,j])
            im=ax.imshow(data,origin='lower')
            if f in aux:
                if 'norm' in aux[f]: im.set_norm(aux[f]['norm']) 
                if 'cmap' in aux[f]: im.set_cmap(aux[f]['cmap']) 
                if 'clim' in aux[f]: im.set_clim(aux[f]['clim']) 
            extent=slc_data[axis+'extent']
            im.set_extent(extent)
            images.append(im)
            if stars:
                if j == 0: 
                    scatter_sp(sp,ax,axis=axis,runaway=False,norm_factor=norm_factor)
                elif j == 1: 
                    scatter_sp(sp,ax,axis=axis,norm_factor=norm_factor)
            ax.set_xlim(extent[0],extent[1])
            ax.set_ylim(extent[2],extent[3])

    gs2 = gridspec.GridSpec(nf+2+stars,1)
    gs2.update(left=0.91,right=0.93,hspace=0.05)
    for j,(im,f) in enumerate(zip(images,fields_to_draw)):
        cax=plt.subplot(gs2[j+stars])
        cbar = fig.colorbar(im,cax=cax,orientation='vertical')
        if f in aux:
            if 'label' in aux[f]: cbar.set_label(aux[f]['label'])
            if 'cticks' in aux[f]: cbar.set_ticks(aux[f]['cticks'])

    if stars:
        cax=plt.subplot(gs2[0])
        cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
             cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
             orientation='vertical')
        cbar.set_label(r'${\rm age [Myr]}$')

    axes=fig.axes[:2*nf]
    if field_label: 
      for ax,f in zip(axes[:nf],fields_to_draw):
        if f in aux:
            if 'label' in aux[f]:
                lab=aux[f]['label']
                label=lab[:lab.rfind(r'\;')]+'$'
                ax.text(0.5,0.95,label,size=20,horizontalalignment='center',
                        transform = ax.transAxes,**(texteffect()))
 
    if stars: 
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
                         scatterpoints = 1, loc='lower left',fontsize='medium',frameon=True)

    plt.setp([ax.get_xticklabels() for ax in axes[:2*nf]],visible=False)
    plt.setp([ax.get_yticklabels() for ax in axes[:2*nf]],visible=False)
    plt.setp(axes[:nf],'ylim',(slc_data['yextent'][2]/zoom,slc_data['yextent'][3]/zoom))

    plt.setp(axes[nf:2*nf],'xlabel','x [kpc]')
    plt.setp(axes[0],'ylabel','z [kpc]')
    if tstamp: 
        plt.setp(axes[0],'title','t=%3d Myr' % tMyr)
    plt.setp(axes[nf],'ylabel','y [kpc]')
    plt.setp([ax.get_xticklabels() for ax in axes[nf:]], visible=True)
    plt.setp([ax.get_yticklabels() for ax in axes[:2*nf:nf]], visible=True)
    plt.setp([ax.xaxis.get_majorticklabels() for ax in axes[nf:2*nf]], rotation=45 )

    pngfname=slcfname[:-1]+'png'
    #canvas = mpl.backends.backend_agg.FigureCanvasAgg(fig)
    #canvas.print_figure(pngfname,num=1,dpi=150,bbox_inches='tight')
    if writefile:
        plt.savefig(pngfname,bbox_inches='tight',num=0,dpi=150)
        plt.close()
    else:
        return fig

def slice2(slcfname,starfname,fields_to_draw,zoom=1.,aux={},\
               writefile=True,tstamp=True,stars=True,field_label=True,norm_factor=2):
    plt.rc('font',size=14)
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)

    slc_data=pickle.load(open(slcfname,'rb'))
    x0=slc_data['yextent'][0]
    y0=slc_data['yextent'][3]
    Lx=slc_data['yextent'][1]-slc_data['yextent'][0]
    Lz=slc_data['yextent'][3]-slc_data['yextent'][2]
    Lz=Lz/zoom
    ix=2
    iz=ix*Lz/Lx
    nf=len(fields_to_draw)
    fig=plt.figure(1,figsize=(ix*nf,iz+ix*1.2))
    gs = gridspec.GridSpec(2,nf,height_ratios=[iz,ix])
    gs.update(top=0.95,left=0.10,right=0.95,wspace=0.05,hspace=0)

    if stars:
        sp=read_starvtk(starfname)

    if 'time' in slc_data:
        tMyr=slc_data['time']
    else:
        time,sp=read_starvtk(starfname,time_out=True)
        tMyr=time*Myr
    images=[]
    for i,axis in enumerate(['y','z']):
        for j,f in enumerate(fields_to_draw):
            ax=plt.subplot(gs[i,j])
            if f is 'star_particles': 
                scatter_sp(sp,ax,axis=axis,norm_factor=norm_factor,type='surf')
                if axis is 'y': 
                    ax.set_xlim(x0,x0+Lx)
                    ax.set_ylim(y0,y0+Lz);
                if axis is 'z': 
                    ax.set_xlim(x0,x0+Lx)
                    ax.set_ylim(x0,x0+Lx)
                ax.set_aspect(1.0)
            else:
                data=slc_data[axis][f]
                im=ax.imshow(data,origin='lower',interpolation='bilinear')
                if f in aux:
                    if 'norm' in aux[f]: im.set_norm(aux[f]['norm']) 
                    if 'cmap' in aux[f]: im.set_cmap(aux[f]['cmap']) 
                    if 'clim' in aux[f]: im.set_clim(aux[f]['clim']) 

                extent=slc_data[axis+'extent']
                im.set_extent(extent)
                images.append(im)
                ax.set_xlim(extent[0],extent[1])
                ax.set_ylim(extent[2],extent[3])

    for j,(im,f) in enumerate(zip(images,fields_to_draw[1:])):
        ax=plt.subplot(gs[0,j+1])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", "3%", pad="1%") 
        cbar = fig.colorbar(im,cax=cax,orientation='horizontal')
        if f in aux:
            if 'label' in aux[f]: cbar.set_label(aux[f]['label'])
            if 'cticks' in aux[f]: cbar.set_ticks(aux[f]['cticks'])
        cax.xaxis.tick_top()
        cax.xaxis.set_label_position('top')

    ax=plt.subplot(gs[0,0])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", "3%", pad="1%") 
    cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
           cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
           orientation='horizontal')
    cax.xaxis.tick_top()
    cax.xaxis.set_label_position('top')
    cbar.set_label(r'${\rm age [Myr]}$')

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
                     scatterpoints = 1, loc='lower left',fontsize='medium',frameon=True)

    axes=fig.axes
    plt.setp([ax.get_xticklabels() for ax in axes[:2*nf]],visible=False)
    plt.setp([ax.get_yticklabels() for ax in axes[:2*nf]],visible=False)
    plt.setp(axes[:nf],'ylim',(slc_data['yextent'][2]/zoom,slc_data['yextent'][3]/zoom))

    plt.setp(axes[nf:2*nf],'xlabel','x [kpc]')
    plt.setp(axes[0],'ylabel','z [kpc]')
    if tstamp: 
        ax=axes[0]
        ax.text(0.5,0.95,'t=%3d Myr' % tMyr,size=16,horizontalalignment='center',
                transform = ax.transAxes,**(texteffect()))
    plt.setp(axes[nf],'ylabel','y [kpc]')
    plt.setp([ax.get_xticklabels() for ax in axes[nf:]], visible=True)
    plt.setp([ax.get_yticklabels() for ax in axes[:2*nf:nf]], visible=True)
    plt.setp([ax.xaxis.get_majorticklabels() for ax in axes[nf:2*nf]], rotation=45 )

    pngfname=slcfname+'ng'
    #canvas = mpl.backends.backend_agg.FigureCanvasAgg(fig)
    #canvas.print_figure(pngfname,num=1,dpi=150,bbox_inches='tight')
    if writefile:
        plt.savefig(pngfname,bbox_inches='tight',num=0,dpi=150)
        plt.close()
    else:
        return fig
