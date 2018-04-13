import numpy as np
import os
import argparse

import rst_handler as rh

def main(**kwargs):
    f_lowres=kwargs['file']
    dir=kwargs['dir']
    id=kwargs['id']
    itime=int(f_lowres[-8:-4])
    print f_lowres,dir,id,itime

    par=rh.parse_par(f_lowres)
    dm=par['domain1']
    Nx=np.array([dm['Nx1'],dm['Nx2'],dm['Nx3']])
    Ng=np.array([dm['NGrid_x1'],dm['NGrid_x2'],dm['NGrid_x3']])
    Nb=Nx/Ng

    grids,NG=rh.calculate_grid(Nx,Nb)

    rstdata_low=rh.read(f_lowres,grids,NG,verbose=True)
    ns=0
    for f in rstdata_low:
        if f.startswith('SCALAR'): ns+=1
    if kwargs['noscalar']: ns=0
    print 'nscalars:',ns
    
    pardata_low=rh.parse_misc_info(f_lowres)
    par=pardata_low['par']

    new_Nx=Nx
    new_NB=Nx
    new_grids,new_NG=rh.calculate_grid(new_Nx,new_NB)

    par=par.replace('NGrid_x1      = %d' % NG[0],'NGrid_x1      = %d' % new_NG[0])
    par=par.replace('NGrid_x2      = %d' % NG[1],'NGrid_x2      = %d' % new_NG[1])
    par=par.replace('NGrid_x3      = %d' % NG[2],'NGrid_x3      = %d' % new_NG[2])
    par=par.replace('AutoWithNProc = %d' % NG[0]*NG[1]*NG[2],'AutoWithNProc = 0')
    pardata_low['par']=par

    rh.write_allfile(pardata_low,rstdata_low,new_grids,\
        dname=dir,id=id,itime=itime,verbose=True,scalar=ns)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('-f','--file',type=str,
                      help='file name for original low resolution')
  parser.add_argument('-i','--id',type=str,
                      help='id of new dataset')
  parser.add_argument('-d','--dir',type=str,
                      help='dir of new dataset')
  parser.add_argument('-c','--crop',action='store_true',help='cropping')
  parser.add_argument('-ns','--noscalar',action='store_true',help='noscalar')
  args = parser.parse_args()
  main(**vars(args))
