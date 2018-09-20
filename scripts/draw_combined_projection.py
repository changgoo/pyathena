import matplotlib as mpl
mpl.use('agg')

import glob,sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,'../')
from pyathena.plot_tools.combined_projection import *

def draw_one(base,pid,i):
    ds=pa.AthenaDataSet('{}/{}/id0/{}.{:04d}.vtk'.format(base,pid,pid,i))
    dproj,Tproj,Bproj,extent=get_proj_ims(ds,slab_range=[5,25],test=True)
    Nz,Nx=dproj.shape
    ratio=Nz/float(Nx)
    fig=plt.figure(0,figsize=(15,15*ratio))
    draw_merged_proj(fig,dproj,Tproj,Bproj,extent)
    fig.savefig('{}{}/proj_figures/{}.merged_proj.{:04d}.png'.format(base,pid,pid,i),dpi=100)

narg=len(sys.argv)
system='tigress'
if narg < 3:
    print 'base directory and problem id should be given'
    sys.exit()

base=sys.argv[1]
pid=sys.argv[2]
files = glob.glob('{}/{}/id0/{}.*.vtk'.format(base,pid,pid))
files.sort()
istr=int(files[0].split('.')[-2])
iend=int(files[-1].split('.')[-2])
print istr,iend

if narg == 4:
    istr=max(istr,int(sys.argv[3]))
    iend=istr+1

if narg == 5:
    iend=min(iend,int(sys.argv[4]))

if narg < 6:
    istride=1
else:
    istride=int(sys.argv[5])

for i in range(istr,iend,istride):
    print '{} of {}'.format(i, len(files))
    draw_one(base,pid,i)
