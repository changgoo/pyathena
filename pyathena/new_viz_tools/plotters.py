import numpy as np

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm,LogNorm,Normalize
import matplotlib.gridspec as gridspec

import cmasher.cm as cma
import cmocean.cm as cmo
from matplotlib.pyplot import cm

from .tigunits import *

plt.rcParams['figure.dpi']=200

def set_axis_color(ax,xy,color):
    '''
        Change spines, label, tick colors
    '''
    if 'x' in xy:
        ax.tick_params(axis='x',colors=color, which='both')
        ax.xaxis.label.set_color(color)
        ax.spines['bottom'].set_color(color)
        ax.spines['top'].set_color(color)
    if 'y' in xy:
        ax.spines['left'].set_color(color)
        ax.spines['right'].set_color(color)
        ax.tick_params(axis='y',colors=color, which='both')
        ax.yaxis.label.set_color(color)

def add_alpha_to_cmap(cmap,expo=0.5,reversed=False):
    '''
        Add a list of alpha to color map
    '''
    my_cmap = cmap(np.arange(cmap.N))

    alpha=np.linspace(0, 1, cmap.N)**expo

    if reversed:
        my_cmap[:,-1] = alpha[::-1]
    else:
        my_cmap[:,-1] = alpha

    my_cmap = mcolors.ListedColormap(my_cmap)
    
    return my_cmap

def plot_sfr(sim,ax):
    hzp=sim.hzp
    ds=sim.ds

    ax.plot(hzp.tMyr,hzp.sfr10)
    ax.axvline(ds.domain['time']*to_Myr,ls='--',color='C1')
    ax.plot(ds.domain['time']*to_Myr,np.interp(ds.domain['time']*to_Myr,hzp.tMyr,hzp.sfr10),'o',color='C1')
    ax.set_xlabel('time [Myr], current time = {:5.1f} Myr'.format(ds.domain['time']*to_Myr))
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_yscale('log')
    ax.set_ylim(5.e-4,5.e-2)
    ax.set_ylabel(r'$\Sigma_{\rm SFR}\,[{\rm M_\odot\,kpc^{-2}\,yr^{-1}}]$')

def plot_Pmid(sim,ax):
    hzp=sim.hzp
    ds=sim.ds

    ax.plot(hzp.tMyr,hzp.Pmid_2p,label=r'$P_{\rm mid}$')
    ax.plot(hzp.tMyr,hzp.W_2p,label=r'$\mathcal{W}$')
    ax.plot(hzp.tMyr,hzp.sfr10*1.e3*eta_conv,label=r'$\eta\Sigma_{\rm SFR}$')
    ax.axvline(ds.domain['time']*to_Myr,ls='--',color='C7')
    #ax.plot(ds.domain['time']*to_Myr,np.interp(ds.domain['time']*to_Myr,hzp.tMyr,hzp.Pmid_2p),'o',color='C1')
    ax.set_xlabel('time [Myr], current time = {:5.1f} Myr'.format(ds.domain['time']*to_Myr))
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    #ax.set_yscale('log')
    #ax.set_ylim(5.e3,1.e5)
    ax.set_ylim(0,5.e4)
    ax.set_ylabel(r'$P/k_B\,[{\rm cm^{-3}\,K}]$')
    ax.legend(loc='upper right')

def plot_outflux(sim,ax,ph,labels=False):
    hzp=sim.hzp
    ds=sim.ds

    if ph == '2p': 
        ax.plot(hzp.tMyr,hzp.massflux_out_05_d_2p,label=r'$\mathcal{F}_M (cool)$',color='C2',lw=3)
        ax.plot(hzp.tMyr,hzp.massflux_out_10_d_2p,label=r'$\mathcal{F}_M (cool)$',color='C2',lw=2)
    if ph == 'h': 
        ax.plot(hzp.tMyr,hzp.massflux_out_05_d_h,label=r'$\mathcal{F}_M (hot)$',color='C3',lw=3)
        ax.plot(hzp.tMyr,hzp.massflux_out_10_d_h,label=r'$\mathcal{F}_M (hot)$',color='C3',lw=2)
    ax.fill_between(hzp.tMyr,hzp.sfr10,label=r'$\eta\Sigma_{\rm SFR}$',color='C0',alpha=0.5,lw=0)
    ax.axvline(ds.domain['time']*to_Myr,ls='--',color='C6',lw=2)
    #ax.plot(ds.domain['time']*to_Myr,np.interp(ds.domain['time']*to_Myr,hzp.tMyr,hzp.Pmid_2p),'o',color='C1')
    if labels:
        ax.set_xlabel('time [Myr], current time = {:5.1f} Myr'.format(ds.domain['time']*to_Myr))
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
    ax.set_yscale('log')
    ax.set_ylim(5.e-4,8.e-2)
    ax.set_xlim(ds.domain['time']*to_Myr-100,ds.domain['time']*to_Myr+100)
    ax.set_ylabel(r'$[{\rm M_\odot\,kpc^{-2}\,yr^{-1}}]$')

def plot_clusters(sim,ax,xaxis,yaxis,yoffset,Ly):
    spx = sim.clpos[xaxis]/1.e3
    spy = sim.clpos[yaxis]/1.e3
    star=ax.scatter(spx,spy+yoffset,marker='o',s=sim.clmass,c=sim.clage,
                    vmax=40,vmin=0,linewidths=0,cmap=sim.clcmap)
    if yoffset !=0: 
        star=ax.scatter(spx,spy+yoffset+Ly,marker='o',s=sim.clmass,c=sim.clage,
                        vmax=40,vmin=0,linewidths=0,cmap=sim.clcmap)

def plot_runaways(sim,ax,xaxis,yaxis,yoffset,Ly):
    rspx = sim.runpos[xaxis]/1.e3
    rspy = sim.runpos[yaxis]/1.e3
    rstar=ax.scatter(rspx,rspy+yoffset,marker='o',s=2,color='k',linewidths=0)
    if yoffset != 0: 
        rstar=ax.scatter(rspx,rspy+yoffset+Ly,marker='o',s=2,color='k',linewidths=0)

def plot_image(sim,ax,cutaxis='z',scafield='nH',vecfield='',
               yshift=0,log=True,colorbars=True,
               clusters=False,runaways=False,
               im_kw=dict(),st_kw=dict()):
    if cutaxis != 'z': yshift=0
    scal=np.roll(sim.slc[cutaxis][scafield],sim.joffset*yshift,axis=0)
    
    if log: scal=np.log10(scal)
    xmin,xmax,ymin,ymax=sim.slc[cutaxis+'extent']
    Ly=ymax-ymin
    im=ax.imshow(scal,extent=[xmin,xmax,ymin,ymax],origin='lower',**im_kw)
    cbar_collection=[im]

    axarr=np.array(['x','y','z'])
    xaxis,yaxis=axarr[axarr != cutaxis]
    if vecfield+'_x' in sim.slc[cutaxis]:
        Ny,Nx=scal.shape
        xfc=np.linspace(xmin,xmax,Nx+1)
        yfc=np.linspace(ymin,ymax,Ny+1)
        xcc=0.5*(xfc[1:]+xfc[:-1])
        ycc=0.5*(yfc[1:]+yfc[:-1])

        vx = np.roll(sim.slc[cutaxis][vecfield+'_'+xaxis],sim.joffset*yshift,axis=0)
        vy = np.roll(sim.slc[cutaxis][vecfield+'_'+yaxis],sim.joffset*yshift,axis=0)
        vmag = np.sqrt(vx**2+vy**2)
        st=ax.streamplot(xcc,ycc,vx,vy,color=vmag,**st_kw)
        cbar_collection.append(st.lines)
    if clusters & (len(sim.sp)>0):
        plot_clusters(sim,ax,xaxis,yaxis,sim.yoffset*yshift/1.e3,Ly*yshift)
    if runaways & (len(sim.sp)>0):
        plot_runaways(sim,ax,xaxis,yaxis,sim.yoffset*yshift/1.e3,Ly*yshift)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_aspect('equal')
    #ax.axis('off')
    if colorbars:
        for coll,loc in zip(cbar_collection,[1,4]):
            if cutaxis == 'z':
                cax = inset_axes(ax,width="30%",height="3%",loc=loc)
                cbar=plt.colorbar(coll,cax=cax,orientation='horizontal')
            else:
                cax = inset_axes(ax,width="10%",height="10%",loc=loc)
                cbar=plt.colorbar(coll,cax=cax,orientation='vertical')
                cax.yaxis.tick_left()
                cax.yaxis.set_label_position('left')
            if (loc == 4) & (cutaxis == 'z'):
                cax.xaxis.tick_top()
                cax.xaxis.set_label_position('top')
            set_axis_color(cax,'xy','k')
            cbar.outline.set_edgecolor('k')
            if scafield in ['nH','pok','mag_pok','ram_pok_z','magnetic_field_strength']:
                set_axis_color(cax,'xy','w')
                cbar.outline.set_edgecolor('w')
    
def plot_surf(sim,ax,scafield='surf',cutaxis='z',vecfield='',
               yshift=0,colorbars=True,log=True,
               clusters=False,runaways=False,
               im_kw=dict(),st_kw=dict()):
    if cutaxis != 'z': yshift=0

    if scafield == 'surf':
        scal=np.roll(sim.surf[cutaxis]['data'],sim.joffset*yshift,axis=0)
    if scafield == 'dproj':
        scal=np.roll(sim.dproj[cutaxis]['data'],sim.joffset*yshift,axis=0)
    if scafield == 'Tproj':
        scal=np.roll(sim.Tproj[cutaxis]['data'],sim.joffset*yshift,axis=0)
    if (scafield == 'Bproj') and (cutaxis == 'y'):
        Bmag = np.sqrt(sim.Bproj['x']['tot']**2+sim.Bproj['y']['tot']**2+sim.Bproj['z']['tot']**2)
        scal=np.roll(Bmag,sim.joffset*yshift,axis=0)
    
    if log: scal=np.log10(scal)
    xmin,xmax,ymin,ymax=sim.surf[cutaxis]['bounds']/1.e3
    Ly=ymax-ymin

    im=ax.imshow(scal,extent=[xmin,xmax,ymin,ymax],origin='lower',**im_kw)
    axarr=np.array(['x','y','z'])
    xaxis,yaxis=axarr[axarr != cutaxis]

    cbar_collection=[im]
    if vecfield+'_x' in sim.slc[cutaxis]:
        Ny,Nx=scal.shape
        xfc=np.linspace(xmin,xmax,Nx+1)
        yfc=np.linspace(ymin,ymax,Ny+1)
        xcc=0.5*(xfc[1:]+xfc[:-1])
        ycc=0.5*(yfc[1:]+yfc[:-1])

        vx = np.roll(sim.slc[cutaxis][vecfield+'_'+xaxis],sim.joffset*yshift,axis=0)
        vy = np.roll(sim.slc[cutaxis][vecfield+'_'+yaxis],sim.joffset*yshift,axis=0)
        vmag = np.sqrt(vx**2+vy**2)
        st=ax.streamplot(xcc,ycc,vx,vy,color=vmag,**st_kw)
        cbar_collection.append(st.lines)
    if (cutaxis == 'y') and hasattr(sim,vecfield):
        Ny,Nx=scal.shape
        xfc=np.linspace(xmin,xmax,Nx+1)
        yfc=np.linspace(ymin,ymax,Ny+1)
        xcc=0.5*(xfc[1:]+xfc[:-1])
        ycc=0.5*(yfc[1:]+yfc[:-1])

        vx = np.roll(sim.Bproj['x']['mean'],sim.joffset*yshift,axis=0)
        vy = np.roll(sim.Bproj['z']['mean'],sim.joffset*yshift,axis=0)
        vmag = np.sqrt(vx**2+vy**2)
        st=ax.streamplot(xcc,ycc,vx,vy,color=vmag,**st_kw)
        cbar_collection.append(st.lines)

    if clusters & (len(sim.sp)>0):
        plot_clusters(sim,ax,xaxis,yaxis,sim.yoffset*yshift/1.e3,Ly*yshift)
    if runaways & (len(sim.sp)>0):
        plot_runaways(sim,ax,xaxis,yaxis,sim.yoffset*yshift/1.e3,Ly*yshift)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    if colorbars:
        for coll,loc in zip(cbar_collection,[1,4]):
            if cutaxis == 'z':
                cax = inset_axes(ax,width="30%",height="3%",loc=loc)
                cbar=plt.colorbar(coll,cax=cax,orientation='horizontal')
            else:
                cax = inset_axes(ax,width="10%",height="10%",loc=loc)
                cbar=plt.colorbar(coll,cax=cax,orientation='vertical')
                cax.yaxis.tick_left()
                cax.yaxis.set_label_position('left')
            if (loc == 4) & (cutaxis == 'z'):
                cax.xaxis.tick_top()
                cax.xaxis.set_label_position('top')
            set_axis_color(cax,'xy','w')
            cbar.outline.set_edgecolor('w')
            if scafield == 'Tproj' and cutaxis=='z':
                set_axis_color(cax,'xy','k')
                cbar.outline.set_edgecolor('k')

    ax.set_aspect('equal')
