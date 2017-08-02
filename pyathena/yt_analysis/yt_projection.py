import matplotlib as mpl
mpl.use('agg')

from ytathena import *
import yt
import glob
import argparse
import os
import cPickle as pickle

import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,SymLogNorm,NoNorm,Normalize
from pyathena import read_starvtk
from slices.scatter_sp import scatter_sp,sp_legend

def do_projection(ds,surfname):
    proj = ds.proj('density',axis='z')
    w=ds.domain_width[0]
    res=(ds.domain_dimensions[0],ds.domain_dimensions[1])
    frb = proj.to_frb(width=w,resolution=res)
    surf=np.array(frb['density'].in_units('Msun/pc**2'))
    if yt.is_root():
        pickle.dump({'data':surf,'bounds':frb.bounds},
          open(surfname,'wb'),pickle.HIGHEST_PROTOCOL)

def main(**kwargs):
 
    dir = kwargs['base_directory']+kwargs['directory']
    fname=glob.glob(dir+'id0/'+kwargs['id']+'.????.vtk')
    fname.sort()

    ngrids=len(glob.glob(dir+'id*/'+kwargs['id']+fname[-9:]))
    comm = yt.communication_system.communicators[-1]
    nprocs = comm.size
    print ngrids,nprocs

    if yt.is_root():
        if not os.path.isdir(dir+'surf/'): os.mkdir(dir+'surf/')
    for f in fname:
        surfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.surf.p'
        if os.path.isfile(surfname):
            print '%s is already there' % surfname
        else:
            if ngrids > nprocs: ds = yt.load(f,units_override=unit_base)
            else: ds = yt.load(f,units_override=unit_base, nprocs=nprocs*8)
            do_projection(ds,surfname)

def plot(**kwargs):

    plt.rc('font',size=11)
    plt.rc('xtick',labelsize=11)
    plt.rc('ytick',labelsize=11)

    dir = kwargs['base_directory']+kwargs['directory']
    fname=glob.glob(dir+'id0/'+kwargs['id']+'.????.vtk')
    fname.sort()

    fig=plt.figure(0,figsize=(5.5,5))
    gs = gridspec.GridSpec(2,2,width_ratios=[1,0.03],wspace=0.0)
    for f in fname:
        surfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.surf.p'
        sp= read_starvtk(f[:-3]+'starpar.vtk')
        frb=pickle.load(open(surfname,'rb'))
        ax=plt.subplot(gs[:,0])
        im=ax.imshow(frb['data'],norm=LogNorm(),origin='lower')
        im.set_extent(np.array(frb['bounds'])/1.e3)
        im.set_cmap(plt.cm.pink_r)
        im.set_clim(1.e-1,1.e2)

	scatter_sp(sp,ax,axis='z',runaway=True,type='surf')
        sp_legend(ax,top=True)

        cax=plt.subplot(gs[0,1])
        cbar = fig.colorbar(im,cax=cax,orientation='vertical')
        cbar.set_label(r'$\Sigma [M_\odot {\rm pc}^{-2}]$')

        cax=plt.subplot(gs[1,1])
        cbar = colorbar.ColorbarBase(cax, ticks=[0,20,40],
               cmap=plt.cm.cool_r, norm=Normalize(vmin=0,vmax=40), 
               orientation='vertical')
        cbar.set_label(r'${\rm age [Myr]}$')

        ax.set_xlabel('x [kpc]')
        ax.set_ylabel('y [kpc]')

        pngfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.surf.png'
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
