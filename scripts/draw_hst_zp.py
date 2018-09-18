import matplotlib as mpl
mpl.use('agg')

from pyathena import preprocessing
from pyathena.set_plt import *
import glob,os
import pandas as pd

base='/u/ckim14/'
dirs=glob.glob('{}/*_newacc*/id0'.format(base))
ids=[]
for dd in dirs:
    ids.append(dd.split('/')[-2])
print ids

for problem_id in ids:
    h_zp=pd.read_pickle('{}{}/hst/{}.hst_zp.p'.format(base,problem_id,problem_id))
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
             ]
    figdir='{}{}/figures/'.format(base,problem_id)
    if not os.path.isdir(figdir): os.mkdir(figdir)
    preprocessing.draw_history(h_zp,metadata,figfname='{}{}/figures/{}-history.png'.format(base,problem_id,problem_id))

