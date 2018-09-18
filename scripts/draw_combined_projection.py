import matplotlib as mpl
mpl.use('agg')

import numpy as np
import matplotlib.pyplot as plt
from pyathena.plot_tools.combined_projection import *

def draw_one(base,pid,i):
    ds=pa.AthenaDataSet('{}/{}/id0/{}.{:04d}.vtk'.format(base,pid,pid,i))
    dproj,Tproj,Bproj,extent=get_proj_ims(ds,slab_range=[5,25],test=True)
    Nz,Nx=dproj.shape
    ratio=Nz/float(Nx)
    fig=plt.figure(0,figsize=(15,15*ratio))
    draw_merged_proj(fig,dproj,Tproj,Bproj,extent)
    fig.savefig('{}{}/proj_figures/{}.merged_proj.{:04d}.png'.format(base,pid,pid,i),dpi=100)

import glob
base = '/u/ckim14/'
pid = 'R4_2pc_newacc'
#draw_one(base,pid,100)
files = glob.glob('{}/{}/id0/{}.*.vtk'.format(base,pid,pid))
files.sort()
for i in range(270,len(files)):
    print '{} of {}'.format(i, len(files))
    draw_one(base,pid,i)
