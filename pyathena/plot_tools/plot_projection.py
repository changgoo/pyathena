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

def plot_projection(surfname,starfname,stars=True,writefile=True,runaway=True,aux={},norm_factor=2.,active=False):

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
    im=ax.imshow(frb['data'],origin='lower')
    im.set_extent(extent)
    if 'norm' in aux: im.set_norm(aux['norm'])
    if 'cmap' in aux: im.set_cmap(aux['cmap'])
    if 'clim' in aux: im.set_clim(aux['clim'])
    ax.text(extent[0]*0.9,extent[3]*0.9,
            't=%3d Myr' % tMyr,ha='left',va='top',**(texteffect()))

    if stars: scatter_sp(sp,ax,axis='z',runaway=runaway,type='surf',norm_factor=norm_factor,active=active)

    cax=plt.subplot(gs[0,1])
    cbar = fig.colorbar(im,cax=cax,orientation='vertical')
    if 'label' in aux: cbar.set_label(aux['label'])

    if stars:
      cax=plt.subplot(gs[1,1])
      cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
             cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
             orientation='vertical')
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

def plot_projection_Z(surfname,starfname,stars=True,writefile=True,runaway=True,norm_factor=2.,aux={}):

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

    extent=np.array(frb['y']['bounds'])/1.e3
    x0=extent[0]
    y0=extent[2]
    Lx=extent[1]-extent[0]
    Lz=extent[3]-extent[2]
 
    ix=2
    iz=ix*Lz/Lx
    fig=plt.figure(0,figsize=(ix*2+0.5,iz))
    gs = gridspec.GridSpec(2,3,width_ratios=[1,1,0.1],wspace=0.0)
    ax1=plt.subplot(gs[:,0])
    im1=ax1.imshow(frb['y']['data'],norm=LogNorm(),origin='lower')
    im1.set_extent(extent)
    if 'cmap' in aux: im1.set_cmap(aux['cmap'])
    if 'clim' in aux: im1.set_clim(aux['clim'])
    ax1.text(extent[0]*0.9,extent[3]*0.9,
            't=%3d Myr' % tMyr,ha='left',va='top',**(texteffect()))
    if stars: scatter_sp(sp,ax1,axis='y',runaway=runaway,type='surf',norm_factor=norm_factor)

    extent=np.array(frb['x']['bounds'])/1.e3
    ax2=plt.subplot(gs[:,1])
    im2=ax2.imshow(frb['x']['data'],norm=LogNorm(),origin='lower')
    im2.set_extent(extent)
    if 'cmap' in aux: im2.set_cmap(aux['cmap'])
    if 'clim' in aux: im2.set_clim(aux['clim'])
    if stars: scatter_sp(sp,ax2,axis='x',runaway=runaway,type='surf')

    cax=plt.subplot(gs[0,2])
    cbar = fig.colorbar(im1,cax=cax,orientation='vertical')
    if 'label' in aux: cbar.set_label(aux['label'])

    if stars:
      cax=plt.subplot(gs[1,2])
      cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
             cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
             orientation='vertical')
      cbar.set_label(r'${\rm age [Myr]}$')
 
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
  stars=True,writefile=True,runaway=True,norm_factor=2.,aux={}):

    plt.rc('font',size=20)
    plt.rc('xtick',labelsize=20)
    plt.rc('ytick',labelsize=20)

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
    cicm._init()
    x=np.arange(cicm.N)
    alphas=0.4*(np.tanh((x-80)/30.)+1)
    #alphas = np.linspace(0.5, 0.5, cicm.N)
    cicm._lut[:-3,-1] = alphas
    cicm._lut[-3,-1] = alphas.min()
    cicm._lut[-2,-1] = alphas.max()

    if 'clim' in aux: clim=aux['clim']
    clim_icm=(0.0,0.5)
    norm_icm=Normalize()

    fig=plt.figure(0,figsize=(ix*2+0.5,iz))
    gs = gridspec.GridSpec(3,3,width_ratios=[1,1,0.05],wspace=0.0)
    ax1=plt.subplot(gs[:,0])
    im1=ax1.imshow(frb['y']['data'],norm=LogNorm(),origin='lower')
    im1.set_extent(extent)
    im1.set_cmap(cism)
    im1.set_clim(clim)

    im11=ax1.imshow(icm['y']['data'],norm=norm_icm,origin='lower')
    im11.set_extent(extent)
    im11.set_cmap(cicm)
    im11.set_clim(clim_icm)

    ax1.text(0.,extent[3]*0.9,
            't=%3d Myr' % tMyr,ha='center',va='top',**(texteffect(20)))
    if stars: scatter_sp(sp,ax1,axis='y',runaway=runaway,type='surf',norm_factor=norm_factor)

    ax2=plt.subplot(gs[:,1])
    im2=ax2.imshow(frb['x']['data'],norm=LogNorm(),origin='lower')
    im2.set_extent(extent)
    im2.set_cmap(cism)
    im2.set_clim(clim)

    im21=ax2.imshow(icm['x']['data'],norm=norm_icm,origin='lower')
    im21.set_extent(extent)
    im21.set_cmap(cicm)
    im21.set_clim(clim_icm)


    if stars: scatter_sp(sp,ax2,axis='x',runaway=runaway,type='surf',norm_factor=norm_factor)


    cax=plt.subplot(gs[0,2])
    cbar = fig.colorbar(im1,cax=cax,orientation='vertical')
    if 'label' in aux: cbar.set_label(aux['label'])

    cax=plt.subplot(gs[1,2])
    cbar = fig.colorbar(im11,cax=cax,orientation='vertical')
    cbar.set_label(r'$C_{\rm ICM}$')


    if stars:
      cax=plt.subplot(gs[2,2])
      cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
             cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
             orientation='vertical')
      cbar.set_label(r'${\rm age [Myr]}$')
 
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
      ax2.set_xlim(x0,x0+Lx)
      ax1.set_ylim(y0,y0+Lz)
      ax2.set_ylim(y0,y0+Lz)
      legend=ax1.legend((s1,s2,s3),(r'$10^3 M_\odot$',r'$10^4 M_\odot$',r'$10^5 M_\odot$'),
                        loc='lower left',fontsize='medium',frameon=True)

    ax1.set_xlabel('x [kpc]')
    ax1.set_ylabel('z [kpc]')

    ax2.set_xlabel('y [kpc]')
    plt.setp(ax2.get_yticklabels(),visible=False)
    pngfname=scalfname+'ng'
    if writefile:
        plt.savefig(pngfname,bbox_inches='tight',num=0,dpi=150)
        plt.close()
    else:
        return fig

def plot_projection_all(surfnames,axis='y',pngfname=None,runaway=True,icm_max=1.0,iscal=-1,norm_factor=2.):

    plt.rc('font',size=30)
    plt.rc('xtick',labelsize=30)
    plt.rc('ytick',labelsize=30)

    nsurf=len(surfnames)
    setup=True
    for i,surfname in enumerate(surfnames):
        scalfnames=glob.glob(surfname.replace('surf.p','scal?.p'))
        scalfnames.sort()
        starfnames=glob.glob(surfname.replace('surf/','starpar/').replace('surf.p','starpar.vtk'))+glob.glob(surfname.replace('surf/','id0/').replace('surf.p','starpar.vtk'))
        nstar=len(starfnames)
        nscal=len(scalfnames)
       
        if nstar > 0:
            starfname=starfnames[0]
            sp=read_starvtk(starfname)
        
        frb=pickle.load(open(surfname,'rb'))
        if nscal > 0:
            icm=pickle.load(open(scalfnames[iscal],'rb'))
            if icm[axis]['data'].max() > (icm_max*1.1):
                print scalfnames[iscal],icm[axis]['data'].min(),icm[axis]['data'].max()

        if 'time' in frb:
            tMyr=frb['time']
        else:
            time,sp=read_starvtk(starfname,time_out=True)
            tMyr=time*Myr
 
        if setup:
            extent=np.array(frb[axis]['bounds'])/1.e3
            x0=extent[0]
            y0=extent[2]
            Lx=extent[1]-extent[0]
            Lz=extent[3]-extent[2]
  
            if axis == 'z': ix=6
            else: ix=3
            iz=ix*Lz/Lx

            cism=plt.cm.bone_r
            cicm=plt.cm.Reds
            cicm._init()
            x=np.arange(cicm.N)
            alphas=0.4*(np.tanh((x-100)/50.)+1)
            #alphas = np.linspace(0.5, 0.5, cicm.N)
            cicm._lut[:-3,-1] = alphas
            cicm._lut[-3,-1] = alphas.min()
            cicm._lut[-2,-1] = alphas.max()
  
            if 'clim' in aux: clim=aux['clim']
            else: clim=(1.e-3,10)
            clim_icm=(0.0,0.5)
            norm_icm=Normalize()

            fig=plt.figure(0,figsize=(ix*nsurf+0.5,iz))
            width_list=[1]*nsurf+[0.05]
            gs = gridspec.GridSpec(3,nsurf+1,width_ratios=width_list,wspace=0.0)
            setup=False

        ax1=plt.subplot(gs[:,i])
        im1=ax1.imshow(frb[axis]['data'],norm=LogNorm(),origin='lower')
        im1.set_extent(extent)
        im1.set_cmap(cism)
        im1.set_clim(clim)

        if nscal > 0:
            icm[axis]['data'] /= icm_max
            im11=ax1.imshow(icm[axis]['data'],norm=norm_icm,origin='lower')
            im11.set_extent(extent)
            im11.set_cmap(cicm)
            im11.set_clim(clim_icm)

        if nstar > 0: 
            scatter_sp(sp,ax1,axis=axis,runaway=runaway,type='surf',norm_factor=norm_factor)

        if i==0:
            ax1.set_title('t=%3dMyr' % tMyr,**(texteffect(30)))
        else:
            ax1.set_title('%3dMyr' % tMyr,**(texteffect(30)))

    axes=fig.axes[:nsurf]
    ax1=axes[-1]
    cax=plt.subplot(gs[0,nsurf])
    cbar = fig.colorbar(im1,cax=cax,orientation='vertical')
    if 'label' in aux: cbar.set_label(aux['label'])

    cax=plt.subplot(gs[1,nsurf])
    cbar = fig.colorbar(im11,cax=cax,orientation='vertical')
    cbar.set_label(r'$C_{\rm ICM}$')

    if nstar > 0:
      cax=plt.subplot(gs[2,nsurf])
      cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
             cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
             orientation='vertical')
      cbar.set_label(r'${\rm age [Myr]}$')
 
      s1=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e3)/norm_factor,color='k',
        alpha=.8,label=r'$10^3 M_\odot$')
      s2=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e4)/norm_factor,color='k',
        alpha=.8,label=r'$10^4 M_\odot$')
      s3=ax1.scatter(Lx*2,Lz*2,
        s=np.sqrt(1.e5)/norm_factor,
        color='k',alpha=.8,label=r'$10^5 M_\odot$')

      starlabels=(r'$10^3 M_\odot$',r'$10^4 M_\odot$',r'$10^5 M_\odot$') 
      if axis=='z':
          legend=ax1.legend((s1,s2,s3),starlabels,
                            loc='upper left',ncol=3,bbox_to_anchor=(0.0, 1.15),
                            fontsize='medium',frameon=True)
      else:
          legend=ax1.legend((s1,s2,s3),starlabels,
                            loc='lower right',fontsize='medium',frameon=True)


    plt.setp(axes,'xlim',(x0,x0+Lx))
    plt.setp(axes,'ylim',(y0,y0+Lz))

    ax1=axes[0]
    if axis=='z':
        #plt.setp(axes,'xlabel','x [kpc]')
        ax1.set_xlabel('x [kpc]')
        ax1.set_ylabel('y [kpc]')
    elif axis=='y':
        #plt.setp(axes,'xlabel','x [kpc]')
        ax1.set_xlabel('x [kpc]')
        ax1.set_ylabel('z [kpc]')
    elif axis=='x':
        #plt.setp(axes,'xlabel','y [kpc]')
        ax1.set_ylabel('y [kpc]')
        ax1.set_ylabel('z [kpc]')


    plt.setp([ax.get_yticklabels() for ax in axes[1:]],visible=False)
    plt.setp([ax.get_xticklabels() for ax in axes[1:]],visible=False)

    if pngfname is not None:
        plt.savefig(pngfname,bbox_inches='tight',num=0,dpi=150)
        plt.close()
    else:
        return fig
