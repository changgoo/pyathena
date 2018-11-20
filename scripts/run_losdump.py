import pyathena.synthetic_observations as syn
from pyathena.synthetic_observations.los_dump_all import los_dump

import numpy as np
import os,glob

base='/tigress/changgoo/'
id='R8_8pc_rst'

Nside=128
center=[0,0,0]
smax=3500
for itime in range(300,301,1):
    fname='%s%s/id0/%s.%4.4d.vtk' % (base,id,id,itime)
    if os.path.isfile(fname):
        print(fname)
        ds,domain=syn.setup_domain(fname,vel=True)
        deltas=domain['dx'][2]
        syn.make_directory(domain,smax,Nside=Nside,center=center)
        #for f in ['density','temperature','velocity1','velocity2','velocity3']:
        for f in ['density','magnetic_field1','magnetic_field2','mangetic_field3']:
            los_dump(ds,domain,deltas,smax,[f],Nside,center,False)
