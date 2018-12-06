from __future__ import print_function
import matplotlib as mpl
mpl.use('agg')

import glob
import sys,os,shutil

sys.path.insert(0,'../')
import pyathena.yt_analysis.phase_plots_with_yt as ph
import yt

narg=len(sys.argv)
system='tigress'
base='/tigress/changgoo/'

# first argument to set system
system = sys.argv[1]
if system == 'pleiades':
    base='/u/ckim14/'
elif system =='tigress':
    base='/tigress/changgoo/'
elif system =='perseus':
    base='/perseus/scratch/gpfs/changgoo/'
elif system =='cori':
    base='/global/cscratch1/sd/changgoo/'
else:
    print('{} is not supported'.format(system))
    sys.exit()

# second argument to set problem ids
dirs=glob.glob('{}/{}'.format(base,sys.argv[2]))

# third and fourth arguments to set time range
t1 = eval(sys.argv[3])
t2 = eval(sys.argv[4])

overwrite=False
if narg > 5:
    overwrite=eval(sys.argv[5])

do_pickling = False
ids=[]
for dd in dirs:
    if os.path.isdir(dd):
        if os.path.isdir(dd+'/id0/'):
            ids.append(os.path.basename(dd))
            do_pickling = True
        else:
            do_pickling = False
            if os.path.isdir(dd+'/phase/'):
                ids.append(os.path.basename(dd))

for pid in ids: 
    outdir='{}{}/phase/'.format(base,pid)
    if not os.path.isdir(outdir): os.mkdir(outdir)
    print(pid)
    if pid.startswith('R8') or pid.startswith('MHD') or pid.startswith('LGR8'):
        params=ph.phase_parameters(Sigma_factor=1.)
    elif pid.startswith('R4') or pid.startswith('LGR4'):
        params=ph.phase_parameters(Sigma_factor=10.)
        print(params.extrema['nH'])
    elif pid.startswith('R2') or pid.startswith('LGR2'):
        params=ph.phase_parameters(Sigma_factor=100.)
    elif pid.startswith('R16'):
        params=ph.phase_parameters(Sigma_factor=0.1)
    else:
        raise ValueError('need to set scale factor for {}'.format(pid))

    for i in range(t1,t2):
        filename='{}/{}/id0/{}.{:04d}.vtk'.format(base,pid,pid,i)
        if os.path.isfile(filename): 
            print('*** Creating joint PDFs for {} ***'.format(os.path.basename(filename)))
            ph.phase_by_slab(filename,zmax=yt.YTQuantity(2048,'pc'),params=params,
                             verbose=60,overwrite=overwrite) 
            for bf in params.bin_fields:
                print('*** Drawing for itime={} fields={bf[0]}-{bf[1]}'.format(i,bf=bf))
                for sidx in [[1,2],[8,9]]:
                    fig=ph.plot_joint_PDFs_slab(base,pid,i,bf,slab_index=sidx)
                    pngfname='{}{}/phase/{}.{:04d}-{b[0]}-{b[1]}-slab{s[0]}-{s[1]}.png'
                    pngfname=pngfname.format(base,pid,pid,i,b=bf,s=sidx)
                    fig.savefig(pngfname,bbox_inches='tight',dpi=150)
