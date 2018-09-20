import matplotlib as mpl
mpl.use('agg')

import glob
import sys,os,shutil
import pandas as pd

sys.path.insert(0,'../')
from pyathena.plot_tools import plot_slices,plot_projection,set_aux
from pyathena.utils import compare_files

narg=len(sys.argv)
system='tigress'
base='/tigress/changgoo/'
if narg > 1:
    system = sys.argv[1]
    if system == 'pleiades':
        base='/u/ckim14/'
    elif system =='tigress':
        base='/tigress/changgoo/'
    else:
        print '{} is not supported'.format(system)
        sys.exit()
if narg > 2:
    dirs=glob.glob('{}/{}'.format(base,sys.argv[2]))
else:
    dirs=glob.glob('{}/*'.format(base))
ids=[]
do_pickling = False
for dd in dirs:
    if os.path.isdir(dd):
        if os.path.isdir(dd+'/id0/'):
            do_pickling = True
        else:
            do_pickling = False
        if os.path.isdir(dd+'/hst/'):
            ids.append(os.path.basename(dd))

for pid in ids:
    print pid
    slc_files=glob.glob('{}{}/slice/{}.????.slice.p'.format(base,pid,pid))
    nf=len(slc_files)
    aux=set_aux.set_aux(pid)
    aux_surf=aux['surface_density']
    field_list=['star_particles','nH','temperature','pok',
           'velocity_z','magnetic_field_strength']
    for itime in range(nf):
        slcname='{}{}/slice/{}.{:04d}.slice.p'.format(base,pid,pid,itime)
        starname='{}{}/starpar/{}.{:04d}.starpar.vtk'.format(base,pid,pid,itime)
        projname='{}{}/surf/{}.{:04d}.surf.p'.format(base,pid,pid,itime)
        if not compare_files(slcname,slcname+'ng'):
            plot_slices.slice2(slcname,starname,field_list,aux=aux)
        if not compare_files(projname,projname+'ng'):
            plot_projection.plot_projection(projname,starname,runaway=False,aux=aux_surf)
