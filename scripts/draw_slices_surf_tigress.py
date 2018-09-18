import matplotlib as mpl
mpl.use('agg')

import glob,os
from pyathena.plot_tools import plot_slices,plot_projection,set_aux
from pyathena.utils import compare_files

base='/tigress/changgoo/'
dirs=glob.glob('{}/*acc*/hst'.format(base))
ids=[]
for dd in dirs:
    ids.append(dd.split('/')[-2])
print ids
ids=['R4_4pc_rst_newacc']

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
