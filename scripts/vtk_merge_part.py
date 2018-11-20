
from __future__ import print_function
from pyathena.parse_par import *
import pyathena as pa
from pyathena.utils import compare_files
import subprocess
import shutil
import numpy as np
import string
import glob
import argparse
import os


def main(**kwargs):
    base = kwargs['base_directory']
    id = kwargs['id']
    if len(kwargs['directory']) > 0:
        dir = kwargs['directory']
    else:
        dir = id
    join_vtk= kwargs['join_vtk']
    if kwargs['var'] != '':
        var=kwargs['var']
        join_vtk=os.path.dirname(join_vtk)+'/'+\
                 os.path.basename(join_vtk).replace('join_vtk','join_vtk_var')
    itime = kwargs['itime']
    istart = kwargs['istart']
    iend = kwargs['iend']+1
    fname='{}{}/id0/{}.{:04d}.vtk'.format(base,dir,id,itime)

    if kwargs['new_id'] != '':
        newid=kwargs['new_id']
    else:
        newid=id

    if kwargs['new_base_directory'] != '':
        newbase='%s%s/' % (kwargs['new_base_directory'],newid)
    else:
        newbase='%s%s/merged/' % (base,dir)

    if not os.path.isdir(newbase): os.mkdir(newbase)

    parfile='%s%s/%s.par' % (base,dir,id)
    par=get_params(parfile)

    NGrids=[int(par['NGrid_x1']),\
            int(par['NGrid_x2']),\
            int(par['NGrid_x3'])]
    Nslab=NGrids[2]
    Nproc=np.prod(NGrids)
    Nproc_h=NGrids[0]*NGrids[1]
    gid=np.arange(Nproc)

    print(Nproc,NGrids)

    print(fname)
    fpath,fbase,fstep,fext,mpi=pa.parse_filename(fname)

    grids=gid[istart:iend]
    baseid=newid
    command=[join_vtk]
    outfile='%s/%s.%s.vtk' % (newbase,baseid,fstep)
    if kwargs['var'] != '':
        command.append('-o %s' % outfile.replace('vtk',kwargs['var'][0]+'.vtk'))
        command.append('-v %s' % var)
    else:
        command.append('-o %s' % outfile)
    zmin=1.e10
    zmax=-1.e10
    for gidx in grids:
        if gidx == 0: 
            vtkfile='%s%s/id%d/%s.%s.%s' % (base,dir,gidx,id,fstep,fext)
        else:
            vtkfile='%s%s/id%d/%s-id%d.%s.%s' % (base,dir,gidx,id,gidx,fstep,fext)
        ds=pa.AthenaDataSet(vtkfile,serial=True)
        zmin=min(ds.domain['left_edge'][2],zmin)
        zmax=max(ds.domain['right_edge'][2],zmax)
        command.append(vtkfile)
    #print command
    print(string.join(command))
    print('id={} to {} corresponds to z={} to {}'.format(istart,iend,zmin,zmax))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
 
    parser.add_argument('-b','--base_directory',type=str,
                        default='/tigress/changgoo/',
                        help='base working directory')
    parser.add_argument('-nb','--new_base_directory',type=str,
                        default='',
                        help='new base directory')
    parser.add_argument('-d','--directory',type=str,default='',
                        help='working directory')
    parser.add_argument('-i','--id',type=str,
                        help='id of dataset')
    parser.add_argument('-ni','--new_id',type=str,
                        default='',
                        help='new id')
    parser.add_argument('-v','--var',type=str,
                        default='',
                        help='variable')
    parser.add_argument('-j','--join_vtk',type=str,
                        default='/tigress/changgoo/join_vtk/join_vtk',
                        help='path to join_vtk excutable')
    parser.add_argument('-is','--istart',type=int,default=0)
    parser.add_argument('-ie','--iend',type=int,default=-1)
    parser.add_argument('-it','--itime',type=int,default=0)
    args = parser.parse_args()
 
    main(**vars(args))

