#!/usr/bin/env python

import sys
import pyathena as pa
import os,glob
import argparse
import sys
import pandas as pd
from pyathena import preprocessing
from pyathena.set_plt import *

def check_files(directory_name,problem_id):
    rstfiles=glob.glob('{}/id0/{}.????.rst'.format(directory_name,problem_id))
    rstfiles+=glob.glob('{}/rst/{}.????.rst'.format(directory_name,problem_id))
    parfile='{}/{}.par'.format(directory_name,problem_id)
    print('==================================================')
    print('checking files...')
    if os.path.isfile(parfile):
        print('par file is already there for {}!'.format(problem_id))
    else:
        if len(rstfiles):
            pa.write_par_from_rst(rstfiles[0],parfile)
        else:
            print('cannot find restart files!')

    hstfiles=glob.glob('{}/id0/{}.hst'.format(directory_name,problem_id))
    hstfiles+=glob.glob('{}/hst/{}.hst'.format(directory_name,problem_id))

    for h in hstfiles:
        print('found history file at {}'.format(h))

    hstfiles=glob.glob('{}/hst/{}.hst_zp.p'.format(directory_name,problem_id))
    for h_zp in hstfiles:
        print('found history file from zprof at {}'.format(h_zp))

    print('==================================================')

    return parfile,h,h_zp

def draw_history(hzpfile,outdir,problem_id):
    if os.path.isfile(hzpfile):
        print('drawing history ...')
        h_zp=pd.read_pickle(hzpfile)
        sfrmean=h_zp['sfr10'].mean()
        snrmean=h_zp['snr10'].mean()
        Pmidmean=h_zp['Pmid_2p'].mean()
        metadata=[
                ('surf',labels['surf']+label_units['surf'],'linear',[],None,None),
                ('sfr10',labels['sfr']+label_units['sfr'],'log',['sfr40','sfr100'],sfrmean*0.05,None),
                ('snr10',labels['snr']+label_units['snr'],'log',['snr40','snr100'],snrmean*0.05,None),
                ('H',labels['scaleH']+label_units['scaleH'],'linear',['H_2p','H_c','H_u','H_w'],None,None),
                ('sigma_eff_2p',labels['sigma_z']+label_units['velocity'],'linear',['v3_2p','vA_2p','cs_2p'],None,None),
                ('Pmid_2p',labels['pressure_mid_2p']+label_units['pressure'],'log',
                 ['Pturb_mid_2p','Pth_mid_2p','Pimag_mean_mid_2p','Pimag_turb_mid_2p'],Pmidmean*0.01,None),
                ('Pmid_2p',labels['pressure']+label_units['pressure'],'linear',['PDE2','W_2p'],None,None),
                ('massflux_bd_d',labels['massflux']+label_units['massflux'],'linear',
                 ['massflux_bd_d_h','massflux_bd_d_2p'],None,None),
                ('snr10',labels['metalflux']+label_units['massflux'],'linear',
                 ['massflux_bd_s1','massflux_bd_s2','massflux_bd_s3'],None,None),
                ('sfr10',labels['massflux']+label_units['massflux'],'log',
                 ['massflux_bd_d_h','massflux_out_5_h','massflux_out_10_h'],sfrmean*0.005,None),
                ('mf_c',labels['massfrac'],'linear',
                 ['mf_u','mf_w'],None,None),
                 ]
        figfname='{}figures/{}-history.png'.format(outdir,problem_id)
        preprocessing.draw_history(h_zp,metadata,figfname)
        print("output: {}".format(figfname))

def print_info(parfile,h,h_zp,pid=None,outdir='./'):
    if pid is None: problem_id=os.path.basename(parfile).replace('.par','')
    else: problem_id=pid
    params=pa.get_params(parfile)
    with open('{}/mds/{}.md'.format(outdir,problem_id),'w') as fp:
        fp.write('## {}\n'.format(problem_id))

        # gas info
        fp.write('### INITIAL CONDITIONS\n')
        fp.write('* Surface Density of Gas (M_sun/pc^2): {surf}\n'.format(**params))
        fp.write('* Surface Density of Stars (M_sun/pc^2): {SurfS}\n'.format(**params))
        fp.write('* Scale Heights of Stars (pc): {zstar}\n'.format(**params))
        fp.write('* Midplane density of DM (M_sun/pc^3): {rhodm}\n'.format(**params))
        fp.write('* Galactocentric Distance (pc): {R0}\n'.format(**params))
        fp.write('* Galactic Rotation (km/s/pc): {Omega}\n'.format(**params))
        fp.write('* Nscalars (Metal tracers): {nscalars}\n'.format(**params))
 
        # domain info
        fp.write('\n### DOMAIN\n')
        for i in range(3): 
            params['Lx{:d}'.format(i+1)]=params['x{:d}max'.format(i+1)]-params['x{:d}min'.format(i+1)]
            params['dx{:d}'.format(i+1)]=params['Lx{:d}'.format(i+1)]/params['Nx{:d}'.format(i+1)]
            params['Lx{:d}'.format(i+1)]=int(params['Lx{:d}'.format(i+1)])
            params['dx{:d}'.format(i+1)]=int(params['dx{:d}'.format(i+1)])
            params['Nx{:d}'.format(i+1)]=int(params['Nx{:d}'.format(i+1)])
        fp.write('* Lx: {Lx1} x {Lx2} x {Lx3}\n'.format(**params))
        fp.write('* Nx: {Nx1} x {Nx2} x {Nx3}\n'.format(**params))
        fp.write('* dx: {dx1} x {dx2} x {dx3}\n'.format(**params))
 
        # perturbation info
        fp.write('\n### SNe\n')
        params['runaways']=False
        if 'fbin' in params:
            if params['fbin']>0: params['runaways']=True
 
        params['random_SN_driving']=False
        if 'snrate_random' in params:
            if (params['iprob'] == 5) & (params['snrate_random'] > 0): params['random_SN_driving']=True
 
        #params['turb_driving']=False
        #if params.has_key('dedt'):
        #    if (params['dedt'] >0) & ((params['iprob'] ==4) | (params['iprob'] == 5)): params['turb_driving']=True
        #    tdrive = params['driving_time']
 
        params['Ia_SN']=False
        if 'Ia_amp' in params:
            if (params['Ia_amp'] == 1): params['Ia_SN']=True
        for k in ['runaways','Ia_SN','random_SN_driving']:
            fp.write('* {}: {}\n'.format(k,params[k]))
 
        fp.write('\n### Movies\n')
        fp.write('* [Slice](../movies/{}_slice_proj.mp4)\n\n'.format(problem_id))
        fp.write('* [Surface Density](../movies/{}_surf.mp4)\n\n'.format(problem_id))
 
        fp.write('\n### History\n')
        fp.write('* `{}`\n\n'.format(h))
        fp.write('* `{}`\n\n'.format(h_zp))
        fp.write('![{} history](../figures/{}-history.png)\n'.format(problem_id,problem_id))
    print('output: {}/mds/{}.md'.format(outdir,problem_id))

def create_indexhtml(outdir):

    from dominate import document
    from dominate import tags

    mds = glob.glob('{}mds/*.md'.format(outdir))
    mds.sort()
    with document(title='TIGRESS Models') as doc:
        tags.h1('TIGRESS Models Quick View')
        for path in mds:
            title = os.path.basename(path).replace('.md','')
            tags.li(tags.a(title, href='./mds/{}.md'.format(title)))

    with open('{}index.html'.format(outdir), 'w') as f:
        f.write(doc.render())

if __name__ == '__main__':
    narg=len(sys.argv)
    system='tigress'
    base='/tigress/changgoo/'
    if narg > 1:
        system = sys.argv[1]
        if system == 'pleiades':
            base='/u/ckim14/'
        elif system =='tigress':
            base='/tigress/changgoo/'
        elif system =='tigress_arm':
            base='/tigress/changgoo/ARM/'
            vtkbase='/tigress/woongkim/TWO/'
        elif system =='perseus':
            base='/perseus/scratch/gpfs/changgoo/'
        elif system =='cori':
            base='/global/cscratch1/sd/changgoo/'
        elif system =='rusty':
            base='/mnt/ceph/users/ckim/'
        else:
            print('{} is not supported'.format(system))
            sys.exit()

    if narg > 2:
        pid = sys.argv[2]
    else:
        print('problem id must be specified')
        sys.exit()

    directory_name = base+pid
    problem_dir = pid
    if os.path.isdir(directory_name+'/slab'):
        directory_name += '/slab/'
        problem_dir += '/slab/'

    preprocessing.doall(base,pid,problem_dir=problem_dir,do_pickling=False)
    parfile,h,h_zp=check_files(directory_name,pid)

    if system == 'tigress':
        outdir='/tigress/changgoo/public_html/TIGRESS_Models/'
    elif system == 'rusty':
        outdir='/mnt/home/ckim/public_www/TIGRESS_Models/'
    else:
        print('{} is not supported for model output'.format(system))
    if not os.path.isdir(outdir): os.mkdir(outdir)
    if not os.path.isdir(outdir+'figures/'): os.mkdir(outdir+'figures/')
    if not os.path.isdir(outdir+'movies/'): os.mkdir(outdir+'movies/')
    if not os.path.isdir(outdir+'mds/'): os.mkdir(outdir+'mds/')
    draw_history(h_zp,outdir,pid)
    print_info(parfile,h,h_zp,pid=pid,outdir=outdir)
    create_indexhtml(outdir)

    if system == 'tigress':
        print('http://tigress-web.princeton.edu/~changgoo/TIGRESS_Models/index.html')
    elif system == 'rusty':
        print('https://users.flatironinstitute.org/~ckim/TIGRESS_Models/index.html')
