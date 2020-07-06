import numpy as np
import os
import glob
import argparse

import rst_handler as rh

def refine(**kwargs):
    f_lowres=kwargs['file']
    pdir=os.path.dirname(f_lowres)+'/'
    pid=os.path.basename(f_lowres).split('.')[0]
    dir=kwargs['dir']
    id=kwargs['id']
    itime=int(f_lowres[-8:-4])
    print(f_lowres,pdir,pid)
    print(dir,id,itime)

    files=glob.glob('{}*{}'.format(f_lowres[:-9],f_lowres[-9:]))
    nf=len(files)

    par=rh.parse_par(f_lowres)
    sghost=kwargs['starghost']
    ns=int(par['configure']['nscalars'])
    if kwargs['nscalars'] >= 0: ns=kwargs['nscalars']
    print(sghost,ns)

    dm=par['domain1']
    Nx=np.array([dm['Nx1'],dm['Nx2'],dm['Nx3']])
    Ng=np.array([dm['NGrid_x1'],dm['NGrid_x2'],dm['NGrid_x3']])
    Nb=Nx/Ng
    if not (Ng.prod() == nf): 
        print("something wrong...",Ng,nf)
        return

    grids,NG=rh.calculate_grid(Nx,Nb)
    if kwargs['split']: # just to split restart files 
        grids_refine,NG_refine=rh.calculate_grid(Nx,[Nb[0]/2,Nb[1]/2,Nb[2]/2])
    else: # refine (2x higher resolution)
        grids_refine,NG_refine=rh.calculate_grid(Nx*2,[Nb[0],Nb[1],Nb[2]])

    pardata=rh.parse_misc_info(f_lowres)

    par=pardata['par'].decode()
    print(par)
    if not kwargs['split']:
        par=par.replace('Nx1           = %d' % Nx[0],'Nx1           = %d' % (Nx[0]*2))
        par=par.replace('Nx2           = %d' % Nx[1],'Nx2           = %d' % (Nx[1]*2))
        par=par.replace('Nx3           = %d' % Nx[2],'Nx3           = %d' % (Nx[2]*2))
    par=par.replace('NGrid_x1      = %d' % NG[0],'NGrid_x1      = %d' % NG_refine[0])
    par=par.replace('NGrid_x2      = %d' % NG[1],'NGrid_x2      = %d' % NG_refine[1])
    par=par.replace('NGrid_x3      = %d' % NG[2],'NGrid_x3      = %d' % NG_refine[2])
    par=par.replace('AutoWithNProc = %d' % NG[0]*NG[1]*NG[2],'AutoWithNProc = 0')
    pardata['par']=par.encode()

    print(par[par.rfind('<domain1'):par.rfind('<problem')])

    for g_orig in grids:
        if g_orig['id'] == 0:
            fname_orig='{}{}.{:04d}.rst'.format(pdir,pid,itime)
        else:
            fname_orig='{}{}-id{}.{:04d}.rst'.format(pdir,pid,g_orig['id'],itime)
        print(fname_orig)
        fm,rstdata=rh.read_rst_grid(fname_orig,starghost=sghost)

        if kwargs['split']:
            is_refine=g_orig['is']
            ie_refine=is_refine+g_orig['Nx']
        else:
            is_refine=g_orig['is']*2
            ie_refine=is_refine+g_orig['Nx']*2

        grids_part=[]
        for g_new in grids_refine:
            if ((g_new['is']-is_refine+1)*(g_new['is']-ie_refine+1) <=0).all():
                #print g
                grids_part.append(g_new)
        if kwargs['split']: 
            rstdata_refine=rstdata
        else:
            rstdata_refine=rh.refine(rstdata,scalar=ns)
        if kwargs['hydro']:
            fc_varnames=['1-FIELD','2-FIELD','3-FIELD']
            for f in fc_varnames:
                if f in rstdata_refine: 
                    print('removing {}'.format(f))
                    del rstdata_refine[f]

        rh.write_allfile(pardata,rstdata_refine,grids_part,grid_disp=is_refine,\
                      id=id,dname=dir,\
                      itime=itime,verbose=True,scalar=ns)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
 
    parser.add_argument('-f','--file',type=str,
                        help='file name for original low resolution')
    parser.add_argument('-i','--id',type=str,
                        help='id of new dataset')
    parser.add_argument('-d','--dir',type=str,
                        help='dir of new dataset')
    parser.add_argument('-ns','--nscalars',type=int,help='number of scalar',default=-1)
    parser.add_argument('-hd','--hydro',action='store_true',help='hydro')
    parser.add_argument('-sp','--split',action='store_true',help='split')
    parser.add_argument('-sg','--starghost',action='store_false',
                        default=True,help='starghost')
    args = parser.parse_args()

    refine(**vars(args))
