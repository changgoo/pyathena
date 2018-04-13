import numpy as np
import os
import argparse

import rst_handler as rh

def crop(data,z1,z2):
    cropped={}
#cropping
    for f in data:
        if f == '3-FIELD':
            cropped[f]=data[f][z1:z2+1,:,:].copy()
        else:
            cropped[f]=data[f][z1:z2,:,:].copy()
        print f,data[f].shape,cropped[f].shape

    return cropped

def main(**kwargs):
    f_lowres=kwargs['file']
    id=kwargs['id']
    print f_lowres,id

    par=rh.parse_par(f_lowres)
    dm=par['domain1']
    Nx=np.array([dm['Nx1'],dm['Nx2'],dm['Nx3']])
    Ng=np.array([dm['NGrid_x1'],dm['NGrid_x2'],dm['NGrid_x3']])
#    Ng=np.array([1,1,1])
    Nb=Nx/Ng

    grids,NG=rh.calculate_grid(Nx,Nb)

    rstdata_low=rh.read(f_lowres,grids,NG,verbose=True)
    ns=0
    for f in rstdata_low:
        if f.startswith('SCALAR'): ns+=1
    if kwargs['noscalar']: ns=0
    print 'nscalars:',ns
    
    #cropping
    if kwargs['crop']: 
        z1=Nx[2]/2-Nx[2]/4
        z2=Nx[2]/2+Nx[2]/4
        x3_orig=[dm['x3min'],dm['x3max']]
        dx=(dm['x3max']-dm['x3min'])/Nx[2]
        xc_pos=np.arange(x3_orig[0],x3_orig[1],dx)+0.5*dx
        x3_new=[xc_pos[z1],xc_pos[z2-1]]
        print z1,z2,xc_pos[0],xc_pos[-1],xc_pos[z1],xc_pos[z2-1]
        rstdata_low=crop(rstdata_low,z1,z2)
    if kwargs['refine']: 
        rstdata_high=rh.refine(rstdata_low,scalar=ns)
    else:
        rstdata_high=rstdata_low
    for i,d in enumerate([rstdata_low,rstdata_high]):
        Bx=d['1-FIELD']
        By=d['2-FIELD']
        Bz=d['3-FIELD']
        dBx = np.diff(Bx,axis=2)
        dBy = np.diff(By,axis=1)
        dBz = np.diff(Bz,axis=0)
        dB=dBx+dBy+dBz
        print np.abs(dB).max(),dB.std()

    pardata_low=rh.parse_misc_info(f_lowres)
    par=pardata_low['par']

    new_Nx=np.array(rstdata_high['DENSITY'].shape)[::-1]
    new_NB=np.array([64,64,64])
    new_grids,new_NG=rh.calculate_grid(new_Nx,new_NB)

    par=par.replace('Nx1           = %d' % Nx[0],'Nx1           = %d' % new_Nx[0])
    par=par.replace('Nx2           = %d' % Nx[1],'Nx2           = %d' % new_Nx[1])
    par=par.replace('Nx3           = %d' % Nx[2],'Nx3           = %d' % new_Nx[2])
    if kwargs['crop']: 
        par=par.replace('x3min         = %d' % x3_orig[0],'x3min         = %d' % x3_new[0])
        par=par.replace('x3max         = %d' % x3_orig[1],'x3max         = %d' % x3_new[1])
    par=par.replace('NGrid_x1      = %d' % NG[0],'NGrid_x1      = %d' % new_NG[0])
    par=par.replace('NGrid_x2      = %d' % NG[1],'NGrid_x2      = %d' % new_NG[1])
    par=par.replace('NGrid_x3      = %d' % NG[2],'NGrid_x3      = %d' % new_NG[2])
    par=par.replace('AutoWithNProc = %d' % NG[0]*NG[1]*NG[2],'AutoWithNProc = 0')
    pardata_low['par']=par

    print par[par.rfind('<domain1'):par.rfind('<problem')]

    rh.write_allfile(pardata_low,rstdata_high,new_grids,\
#    dname='/nobackup/ckim14/%s/rst/' % id,id=id,verbose=True,scalar=ns)
#    dname='/u/ckim14/rst/',id=id,verbose=True,scalar=ns)
#    dname='/tigress/PERSEUS/changgoo/%s/rst/' % id,id=id,verbose=True,scalar=ns)
    dname='/scratch/gpfs/changgoo/%s/rst/' % id,id=id,verbose=True,scalar=ns)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('-f','--file',type=str,
                      help='file name for original low resolution')
  parser.add_argument('-i','--id',type=str,
                      help='id of dataset')
  parser.add_argument('-c','--crop',action='store_true',help='cropping')
  parser.add_argument('-r','--refine',action='store_true',help='refining')
  parser.add_argument('-ns','--noscalar',action='store_true',help='noscalar')
  args = parser.parse_args()
  main(**vars(args))
