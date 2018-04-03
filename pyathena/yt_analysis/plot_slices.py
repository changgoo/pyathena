import matplotlib as mpl
#mpl.use('Agg')

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

fields=['cell_volume','cell_mass']
unit=set_units(muH=1.4271)
Myr=unit['time'].to('Myr').value

def projection(sp,axis):
    if axis == 0 or axis == 'z':
        spx=sp['x1']
        spy=sp['x2']
        spz=sp['x3']
    elif axis == 1 or axis == 'y':
        spx=sp['x1']
        spy=sp['x3']
        spz=sp['x2']
    elif axis == 2 or axis == 'x':
        spx=sp['x2']
        spy=sp['x3']
        spz=sp['x1']
    return spx,spy,spz

def projection_v(sp,axis):
    if axis == 0 or axis == 'z':
        spx=sp['v1']
        spy=sp['v2']
        spz=sp['v3']
    elif axis == 1 or axis == 'y':
        spx=sp['v1']
        spy=sp['v3']
        spz=sp['v2']
    elif axis == 2 or axis == 'x':
        spx=sp['v2']
        spy=sp['v3']
        spz=sp['v1']
    return spx,spy,spz

def scatter_sp(sp,ax,axis=0,thickness=10,norm_factor=4., \
  type='slice',kpc=True,runaway=True):
    unit=set_units(muH=1.4271)
    Msun=unit['mass'].to('Msun').value
    Myr=unit['time'].to('Myr').value
    #print len(sp)
    if len(sp) >0:
      runaways=(sp['mass'] == 0.0)
      sp_runaway=sp[runaways]
      sp_normal=sp[-runaways]
      #print len(sp_runaway)
      if len(sp_runaway) > 0 and runaway:
        spx,spy,spz=projection(sp_runaway,axis)
        spvx,spvy,spvz=projection_v(sp_runaway,axis)
        if kpc:
            spx = spx/1.e3
            spy = spy/1.e3
        if type == 'slice': 
            islab=np.where(abs(spz) < thickness)
            #ax.scatter(spx.iloc[islab],spy.iloc[islab],marker='.',c='k',alpha=1.0)
            #ax.quiver(spx.iloc[islab],spy.iloc[islab],
            #      spvx.iloc[islab],spvy.iloc[islab],color='k',alpha=1.0)
        ax.scatter(spx,spy,marker='o',color='k',alpha=0.5,s=10.0/norm_factor)
        #ax.quiver(spx,spy,spvx,spvy,color='w',alpha=0.5)
 
      if len(sp_normal) > 0: 
        spx,spy,spz=projection(sp_normal,axis)
        if kpc:
            spx = spx/1.e3
            spy = spy/1.e3
        if type == 'slice': xbool=abs(spz) < thickness
        spm=np.sqrt(sp_normal['mass']*Msun)/norm_factor
        spa=sp_normal['age']*Myr
        iyoung=np.where(spa < 40.)
        #print len(iyoung[0])
        if type == 'slice':
            islab=np.where(xbool*(spa<40))
            ax.scatter(spx.iloc[islab],spy.iloc[islab],marker='o',\
                s=spm.iloc[islab],c=spa.iloc[islab],\
                vmax=40,vmin=0,cmap=plt.cm.cool_r,alpha=1.0)
        ax.scatter(spx.iloc[iyoung],spy.iloc[iyoung],marker='o',\
            s=spm.iloc[iyoung],c=spa.iloc[iyoung],\
            vmax=40,vmin=0,cmap=plt.cm.cool_r,alpha=0.5)


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

def slice1(slcfname,vtkfname,fields_to_draw,zoom=1.,\
               writefile=True,tstamp=True,stars=True,field_label=True):
    aux=ya.set_aux(os.path.basename(slcfname))
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

    sp=read_starvtk(vtkfname[:-3]+'starpar.vtk')
    if 'time' in slc_data:
        tMyr=slc_data['time']
    else:
        time,sp=read_starvtk(vtkfname[:-3]+'starpar.vtk',time_out=True)
        tMyr=time*Myr
    images=[]
    for i,axis in enumerate(['y','z']):
        for j,f in enumerate(fields_to_draw):
            data=slc_data[axis][f]
            ax=plt.subplot(gs[i,j])
            im=ax.imshow(data,origin='lower')
            if aux[f]['log']: im.set_norm(LogNorm()) 
            extent=slc_data[axis+'extent']
            im.set_extent(extent)
            im.set_cmap(aux[f]['cmap'])
            im.set_clim(aux[f]['clim'])
            images.append(im)
            if stars:
              if j == 0: 
                scatter_sp(sp,ax,axis=axis,runaway=False,norm_factor=4.)
              elif j == 1: 
                scatter_sp(sp,ax,axis=axis,norm_factor=4.)
            ax.set_xlim(extent[0],extent[1])
            ax.set_ylim(extent[2],extent[3])

    gs2 = gridspec.GridSpec(nf+2+stars,1)
    gs2.update(left=0.91,right=0.93,hspace=0.05)
    for j,(im,f) in enumerate(zip(images,fields_to_draw)):
        cax=plt.subplot(gs2[j+stars])
        cbar = fig.colorbar(im,cax=cax,orientation='vertical')
        cbar.set_label(aux[f]['label'])
        #cax.xaxis.tick_top()
        #cax.xaxis.set_label_position('top')
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
        lab=aux[f]['label']
        label=lab[:lab.rfind(r'\;')]+'$'
        ax.text(0.5,0.95,label,size=20,horizontalalignment='center',
                transform = ax.transAxes,**(texteffect()))
 
    if stars: legend=sp_legend(axes[-1],top=False,norm_factor=4.)

    plt.setp([ax.get_xticklabels() for ax in axes[:2*nf]],visible=False)
    plt.setp([ax.get_yticklabels() for ax in axes[:2*nf]],visible=False)
    plt.setp(axes[:nf],'ylim',(slc_data['yextent'][2]/zoom,slc_data['yextent'][3]/zoom))

    plt.setp(axes[nf:2*nf],'xlabel','x [kpc]')
    plt.setp(axes[0],'ylabel','z [kpc]')
    if tstamp: 
#      axes[0].text(x0*0.9,y0*0.9,
#                   't=%3d Myr' % tMyr,ha='left',va='top',**(texteffect(16)))
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

def slice2(slcfname,starfname,fields_to_draw,zoom=1.,\
               writefile=True,tstamp=True,stars=True,field_label=True):
    aux=ya.set_aux(os.path.basename(slcfname))
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
    norm_factor=2.

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
              if aux[f]['log']: im.set_norm(LogNorm()) 
              extent=slc_data[axis+'extent']
              im.set_extent(extent)
              im.set_cmap(aux[f]['cmap'])
              im.set_clim(aux[f]['clim'])
              images.append(im)
              ax.set_xlim(extent[0],extent[1])
              ax.set_ylim(extent[2],extent[3])

    for j,(im,f) in enumerate(zip(images,fields_to_draw[1:])):
        ax=plt.subplot(gs[0,j+1])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", "3%", pad="1%") 
#        cax=plt.subplot(gs[0,j])
        cbar = fig.colorbar(im,cax=cax,orientation='horizontal')
        cbar.set_label(aux[f]['label'])
        cax.xaxis.tick_top()
        cax.xaxis.set_label_position('top')
        if 'cticks' in aux[f]: cbar.set_ticks(aux[f]['cticks'])

    ax=plt.subplot(gs[0,0])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", "3%", pad="1%") 
#    cax=plt.subplot(gs[0,0])
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
#      axes[0].text(x0*0.9,y0*0.9,
#                   't=%3d Myr' % tMyr,ha='left',va='top',**(texteffect(16)))
#      plt.setp(axes[0],'title','t=%3d Myr' % tMyr)
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

def slice3(slcfname,fields_to_draw,axis='y',extent=None,\
               writefile=True,tstamp=True,field_label=True):
    aux=ya.set_aux(os.path.basename(slcfname))

    slc_data=pickle.load(open(slcfname,'rb'))
    x0=slc_data[axis+'extent'][0]
    y0=slc_data[axis+'extent'][2]
    print(x0,y0)
    Lx=slc_data[axis+'extent'][1]-slc_data[axis+'extent'][0]
    Ly=slc_data[axis+'extent'][3]-slc_data[axis+'extent'][2]
    if extent is None: extent=[x0,x0+Lx,y0,y0+Ly]
    x0=extent[0]
    y0=extent[2]
    lx=extent[1]-extent[0]
    ly=extent[3]-extent[2]
    print(extent,lx,ly)
    ix=2
    iz=ix*ly/lx
    nf=len(fields_to_draw)
    fig=plt.figure(1,figsize=(ix*nf,iz+ix*1.2))
    gs = gridspec.GridSpec(1,nf)
    gs.update(top=0.95,left=0.10,right=0.95,wspace=0.05,hspace=0)
    norm_factor=2.

    starname=slcfname.replace('slice/','starpar/').replace('slice.p','starpar.vtk')
    sp=read_starvtk(starname)
    if 'time' in slc_data:
        tMyr=slc_data['time']
    else:
        time,sp=read_starvtk(starname,time_out=True)
        tMyr=time*Myr
    images=[]
    star_axis=-1
    for j,f in enumerate(fields_to_draw):
        ax=plt.subplot(gs[0,j])
        if f is 'star_particles': 
          scatter_sp(sp,ax,axis=axis,norm_factor=norm_factor,type='surf')
          ax.set_xlim(x0,x0+lx)
          ax.set_ylim(y0,y0+ly);
          ax.set_aspect(1.0)
          star_axis=j
        else:
          data=slc_data[axis][f]
          im=ax.imshow(data,origin='lower',interpolation='bilinear')#interpolation='nearest',resample=True)
          if aux[f]['log']: im.set_norm(LogNorm()) 
          im.set_extent(slc_data[axis+'extent'])
          im.set_cmap(aux[f]['cmap'])
          im.set_clim(aux[f]['clim'])
          images.append(im)
          ax.set_xlim(extent[0],extent[1])
          ax.set_ylim(extent[2],extent[3])
          ax.set_aspect(1.0)

    for j,(im,f) in enumerate(zip(images,fields_to_draw)):
        if f != 'star_particles':
            ax=plt.subplot(gs[0,j])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("top", "3%", pad="1%") 
##           cax=plt.subplot(gs[0,j])
            cbar = fig.colorbar(im,cax=cax,orientation='horizontal')
            cbar.set_label(aux[f]['label'])
            cax.xaxis.tick_top()
            cax.xaxis.set_label_position('top')
            if 'cticks' in aux[f]: cbar.set_ticks(aux[f]['cticks'])

    if star_axis != -1:
        ax=plt.subplot(gs[0,star_axis])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", "3%", pad="1%") 
##       cax=plt.subplot(gs[0,0])
        cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
               cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
               orientation='horizontal')
        cax.xaxis.tick_top()
        cax.xaxis.set_label_position('top')
        cbar.set_label(r'${\rm age [Myr]}$')
 
        s1=ax.scatter(Lx*2,Ly*2,
          s=np.sqrt(1.e3)/norm_factor,color='k',
          alpha=.8,label=r'$10^3 M_\odot$')
        s2=ax.scatter(Lx*2,Ly*2,
          s=np.sqrt(1.e4)/norm_factor,color='k',
          alpha=.8,label=r'$10^4 M_\odot$')
        s3=ax.scatter(Lx*2,Ly*2,
          s=np.sqrt(1.e5)/norm_factor,
          color='k',alpha=.8,label=r'$10^5 M_\odot$')
 
        ax.set_xlim(x0,x0+lx)
        ax.set_ylim(y0,y0+ly);
        legend=ax.legend((s1,s2,s3),(r'$10^3 M_\odot$',r'$10^4 M_\odot$',r'$10^5 M_\odot$'), 
                         scatterpoints = 1, loc='lower left',fontsize='medium',frameon=True)

    axes=fig.axes
    plt.setp([ax.get_xticklabels() for ax in axes[:nf]],visible=False)
    plt.setp([ax.get_yticklabels() for ax in axes[:nf]],visible=False)
    plt.setp(axes[:nf],'xlim',(x0,x0+lx))
    plt.setp(axes[:nf],'ylim',(y0,y0+ly))

    plt.setp(axes[0],'xlabel','x [kpc]')
    plt.setp(axes[0],'ylabel','z [kpc]')
    if tstamp: 
      ax=axes[0]
      ax.text(0.5,0.95,'t=%3d Myr' % tMyr,size=16,horizontalalignment='center',
              transform = ax.transAxes,**(texteffect()))
#      axes[0].text(x0*0.9,y0*0.9,
#                   't=%3d Myr' % tMyr,ha='left',va='top',**(texteffect(16)))
#      plt.setp(axes[0],'title','t=%3d Myr' % tMyr)
    plt.setp([ax.get_xticklabels() for ax in axes[:nf:nf]], visible=True)
    plt.setp([ax.get_yticklabels() for ax in axes[:nf:nf]], visible=True)
    plt.setp([ax.xaxis.get_majorticklabels() for ax in axes[:nf:nf]], rotation=45 )

    pngfname=slcfname+'ng'
    #canvas = mpl.backends.backend_agg.FigureCanvasAgg(fig)
    #canvas.print_figure(pngfname,num=1,dpi=150,bbox_inches='tight')
    if writefile:
        plt.savefig(pngfname,bbox_inches='tight',num=0,dpi=150)
        plt.close()
    else:
        return fig
