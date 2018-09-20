import numpy as np
import os
import glob
import argparse

import rst_handler as rh

def degrade(**kwargs):
    f_orig=kwargs['file']
    dir=kwargs['dir']
    id=kwargs['id']
    itime=int(f_orig[-8:-4])
    print f_orig,dir,id,itime

    par=rh.parse_par(f_orig)
    dm=par['domain1']
    Nx=np.array([dm['Nx1'],dm['Nx2'],dm['Nx3']])
    Ng=np.array([dm['NGrid_x1'],dm['NGrid_x2'],dm['NGrid_x3']])
    Nb=Nx/Ng

    grids,NG=rh.calculate_grid(Nx,Nb)

    rstdata=rh.read(f_orig,grids,NG,verbose=True)
    ns=0
    for f in rstdata:
        if f.startswith('SCALAR'): ns+=1
    if kwargs['noscalar']: ns=0
    print 'nscalars:',ns
    
    rstdata_degrade=rh.degrade(rstdata,scalar=ns)
    if kwargs['hydro']:
        fc_varnames=['1-FIELD','2-FIELD','3-FIELD']
        for f in fc_varnames:
            if f in rstdata_degrade: 
                print 'removing {}'.format(f)
                del rstdata_degrade[f]
    grids_deg,NG_deg=rh.calculate_grid(Nx/2,[32,32,64])

    pardata=rh.parse_misc_info(f_orig)
    par=pardata['par']
    par=par.replace('Nx1           = %d' % Nx[0],'Nx1           = %d' % (Nx[0]/2))
    par=par.replace('Nx2           = %d' % Nx[1],'Nx2           = %d' % (Nx[1]/2))
    par=par.replace('Nx3           = %d' % Nx[2],'Nx3           = %d' % (Nx[2]/2))
    par=par.replace('NGrid_x1      = %d' % NG[0],'NGrid_x1      = %d' % NG_deg[0])
    par=par.replace('NGrid_x2      = %d' % NG[1],'NGrid_x2      = %d' % NG_deg[1])
    par=par.replace('NGrid_x3      = %d' % NG[2],'NGrid_x3      = %d' % NG_deg[2])
    par=par.replace('AutoWithNProc = %d' % NG[0]*NG[1]*NG[2],'AutoWithNProc = 0')
    pardata['par']=par

    rh.write_allfile(pardata,rstdata_degrade,grids_deg,\
        dname=dir,id=id,itime=itime,verbose=True,scalar=ns)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
 
    parser.add_argument('-f','--file',type=str,
                        help='file name for original low resolution')
    parser.add_argument('-i','--id',type=str,
                        help='id of new dataset')
    parser.add_argument('-d','--dir',type=str,
                        help='dir of new dataset')
    parser.add_argument('-ns','--noscalar',action='store_true',help='noscalar')
    parser.add_argument('-hd','--hydro',action='store_true',help='hydro')
    args = parser.parse_args()
    degrade(**vars(args))
