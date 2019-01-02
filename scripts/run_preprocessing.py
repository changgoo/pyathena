import matplotlib as mpl
mpl.use('agg')

import glob
import sys,os,shutil
import pandas as pd

sys.path.insert(0,'../')
from pyathena import preprocessing
from pyathena.set_plt import *

narg=len(sys.argv)
system='tigress'
base='/tigress/changgoo/'
if narg > 1:
    system = sys.argv[1]
    if system == 'pleiades':
        base='/u/ckim14/'
    elif system =='tigress':
        base='/tigress/changgoo/'
    elif system =='perseus':
        base='/perseus/scratch/gpfs/changgoo/'
    elif system =='cori':
        base='/global/cscratch1/sd/changgoo/'
    elif system =='rusty':
        base='/mnt/ceph/users/ckim/'
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
            ids.append(os.path.basename(dd))
            do_pickling = True
        else:
            do_pickling = False
            if os.path.isdir(dd+'/hst/'):
                ids.append(os.path.basename(dd))

if narg > 3:
    do_pickling = eval(sys.argv[3])

print system,base,dd,ids,do_pickling
for problem_id in ids: 
    if (system == 'cori') | (system == 'rusty'):
        preprocessing.doall(base,problem_id,problem_dir='{}/slab/'.format(problem_id),\
                            do_pickling=do_pickling,use_yt=False)
    else:
        preprocessing.doall(base,problem_id,do_pickling=do_pickling,use_yt=False)

for problem_id in ids:
    if (system == 'cori') | (system == 'rusty'):
        problem_dir='{}/slab/'.format(problem_id)
    else:
        problem_dir='{}/'.format(problem_id)
    print 'drawing {} ...'.format(problem_id)
    h_zp=pd.read_pickle('{}{}hst/{}.hst_zp.p'.format(base,problem_dir,problem_id))
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
            ('sfr10',labels['massflux']+label_units['massflux'],'linear',
             ['massflux_bd_d_h','massflux_out_5_h','massflux_out_10_h'],None,None),
            ('mf_c',labels['massfrac'],'linear',
             ['mf_u','mf_w'],None,None),
             ]
    figdir='{}{}/figures/'.format(base,problem_dir)
    if not os.path.isdir(figdir): os.mkdir(figdir)
    figfname='{}{}/figures/{}-history.png'.format(base,problem_dir,problem_id)
    preprocessing.draw_history(h_zp,metadata,figfname)
    if system == 'tigress':
        figdir_tigress='{}/public_html/TIGRESS_figures/history/'.format(base)
        if len(figfname) > 0: 
            os.chmod(figfname,0644)
            shutil.copy2(figfname,figdir_tigress)
