import matplotlib as mpl
mpl.use('agg')

import glob
import sys,os,shutil
import pandas as pd
import numpy as np

sys.path.insert(0,'../')
from pyathena import get_params
from pyathena.plot_tools import plot_slices,plot_projection,set_aux
import pyathena.plot_tools.movie as movie
from pyathena.utils import compare_files
import pickle as p

narg=len(sys.argv)
system='tigress'
base='/tigress/changgoo/'
if narg > 1:
    system = sys.argv[1]
    if system == 'pleiades':
        base='/u/ckim14/'
    elif system =='tigress':
        base='/tigress/changgoo/'
    elif system =='tigress_rps':
        base='/tigress/changgoo/'
    elif system =='tigress_arm':
        base='/tigress/changgoo/ARM/'
    elif system =='rusty':
        base='/mnt/ceph/users/ckim/'
    elif system =='cori':
        base='/global/cscratch1/sd/changgoo/'
    else:
        print('{} is not supported'.format(system))
        sys.exit()
if narg > 2:
    dirs=glob.glob('{}/{}'.format(base,sys.argv[2]))
else:
    dirs=glob.glob('{}/*'.format(base))
ids=[]
for dd in dirs:
    if os.path.isdir(dd):
        if os.path.isdir(dd+'/slice/'):
            ids.append(os.path.basename(dd))

if narg > 3:
    overwrite = eval(sys.argv[3])
else:
    overwrite = False

for pid in ids:
    print(pid)
    par = get_params('{}{}/{}.par'.format(base,pid,pid))
    if 'pattern' in par:
        vy0 = par['Omega']*(1.0-par['pattern'])*par['R0']
        print('v_y,circ = {}'.format(vy0))
    else:
        vy0 = 0.0
    slc_files=glob.glob('{}{}/slice/{}.????.slice.p'.format(base,pid,pid))
    slc_files.sort()
    nf=len(slc_files)
    aux=set_aux.set_aux(pid)
    aux_surf=aux['surface_density']
    if system == 'tigress_rps':
        field_list=['star_particles','surface_density','specific_scalar3_proj','nH','specific_scalar3','temperature','pok','ram_pok_z']
    else:
        field_list=['star_particles','surface_density','nH','temperature','pok','velocity_z']
    slcdata=p.load(open(slc_files[0],'rb'), encoding='latin1')
    if 'magnetic_field_strength' in slcdata['x']:
        if system == 'tigress_rps':
            field_list += ['mag_pok']
        else:
            field_list += ['magnetic_field_strength']
    for slcname in slc_files:
        print(slcname)
        starname=slcname.replace('slice.p','starpar.vtk').replace('slice','starpar')
        projname=slcname.replace('slice','surf')
        if not compare_files(slcname,slcname.replace('.p','_proj.png')) or overwrite:
            #plot_slices.slice2(slcname,starname,field_list,aux=aux,vy0=vy0)
            plot_slices.slice_proj(slcname,projname,starname,field_list,aux=aux,vy0=vy0)
        #if not compare_files(projname,projname+'ng') or overwrite:
            plot_projection.plot_projection(projname,starname,
              scale_func=np.cbrt,runaway=False,aux=aux_surf,vy0=vy0)

    if system.startswith('tigress'):
        basedir1='{}{}/'.format(base,pid)
        basedir2='/tigress/changgoo/public_html/TIGRESS_figures/movies/'
        if system == 'tigress_arm': basedir2 += 'ARM/'
        ffig = os.path.join(basedir1,'slice/*.slice_proj.png')
        fmp4 = os.path.join(basedir1,'slice/{}_slice_proj.mp4'.format(pid))
        movie.make_movie(ffig, fmp4)
        shutil.copy(fmp4,basedir2)
        ffig = os.path.join(basedir1,'surf/*.surf.png')
        fmp4 = os.path.join(basedir1,'slice/{}_surf.mp4'.format(pid))
        movie.make_movie(ffig, fmp4)
        shutil.copy(fmp4,basedir2)
