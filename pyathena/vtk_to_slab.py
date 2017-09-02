from parse_par import parse_par
from vtk_reader import parse_filename
from utils import compare_files
import subprocess
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

    files=glob.glob('%s%s/id0/%s.????.vtk' % (base,dir,id))
    files.sort()

    if kwargs['range'] != '':
        sp=kwargs['range'].split(',')
        start = eval(sp[0])
        end = eval(sp[1])
        fskip = eval(sp[2])
    else:
        start = 0
        end = len(files)
        fskip = 1

    files=files[start:end:fskip]
    if kwargs['new_base_directory'] != '':
        newbase='%s%s/' % (kwargs['new_base_directory'],dir)
    else:
        newbase='%s%s/slab/' % (base,dir)

    rstfiles=glob.glob('%s%s/id0/%s.????.rst' % (base,dir,id))
    rstfiles+=glob.glob('%s%s/rst/%s.????.rst' % (base,dir,id))

    par,block,field=parse_par(rstfiles[0])

    NGrids=[eval(par['domain1']['NGrid_x1'][0]),\
            eval(par['domain1']['NGrid_x2'][0]),\
            eval(par['domain1']['NGrid_x3'][0])]
    Nslab=NGrids[2]
    Nproc=np.prod(NGrids)
    Nproc_h=NGrids[0]*NGrids[1]
    gid=np.arange(Nproc)

    print Nproc,NGrids

    if not os.path.isdir(newbase): os.mkdir(newbase)

    for f in files:
        print f
        fpath,fbase,fstep,fext,mpi=parse_filename(f)

        for islab in range(Nslab):
            print '%d of %d' % (islab, Nslab)
            grids=gid[gid/Nproc_h == islab]
            if islab == 0: baseid=id
            else: baseid='%s-id%d' %(id,islab)
            if not os.path.isdir('%s/id%d' % (newbase,islab)):
                os.mkdir('%s/id%d' % (newbase,islab))
            command=[join_vtk]
            outfile='%s/id%d/%s.%s.vtk' % (newbase,islab,baseid,fstep)
            command.append('-o %s' % outfile)
            for gidx in grids:
                if gidx == 0: 
                    vtkfile='%s%s/id%d/%s.%s.%s' % (base,dir,gidx,id,fstep,fext)
                else:
                    vtkfile='%s%s/id%d/%s-id%d.%s.%s' % (base,dir,gidx,id,gidx,fstep,fext)
                command.append(vtkfile)
            #print command
            if not compare_files(vtkfile,outfile) or kwargs['overwrite']:
                subprocess.call(string.join(command),shell=True)
            else:
                print '%s is newer than %s' % (outfile, vtkfile)
# delete originals
        file_originals=glob.glob('%s/id*/%s-id*.%s.%s' % (fpath,fbase,fstep,fext))
        if len(file_originals) > 0: 
            for f in file_originals: os.remove(f)
            os.remove('%s/id0/%s.%s.%s' % (fpath,fbase,fstep,fext))
# copy starpar.vtk
        src_starpar_name='%s/id0/%s.%s.starpar.vtk' % (fpath,fbase,fstep)
        dst_name='%s/id0/' % (newbase)
        if os.path.isfile(src_starpar_name): 
            command=['mv',src_starpar_name,dst_name]
            subprocess.call(string.join(command),shell=True)

# move zprof
        src_zprof_names=glob.glob('%s/id0/%s.%s.*.zprof' % (fpath,fbase,fstep))
        dst_name='%s/id0/' % (newbase)
        for f in src_zprof_names:
            if os.path.isfile(f):
                command=['mv',f,dst_name]
                subprocess.call(string.join(command),shell=True)
# copy history
    src_hst_name='%s/id0/%s.hst' % (fpath,fbase)
    dst_name='%s/id0/' % (newbase)
    if os.path.isfile(src_hst_name):
        command=['cp',src_hst_name,dst_name]
        subprocess.call(string.join(command),shell=True)

    src_hst_name='%s/id0/%s.sn' % (fpath,fbase)
    dst_name='%s/id0/' % (newbase)
    if os.path.isfile(src_hst_name):
        command=['cp',src_hst_name,dst_name]
        subprocess.call(string.join(command),shell=True)

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
    parser.add_argument('-j','--join_vtk',type=str,
                        default='/tigress/changgoo/join_vtk/join_vtk',
                        help='path to join_vtk excutable')
    parser.add_argument('-o','--overwrite',action='store_true',
                        help='overwrite')
    parser.add_argument('-r','--range',type=str,default='',
                       help='time range, start:end:skip')
    args = parser.parse_args()
 
    main(**vars(args))
