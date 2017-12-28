#from startup import *
import sys
import subprocess

sysconfig={}
sysconfig['ccalin004']={'base':'/mnt/ceph/users/ckim/','source':'/mnt/home/ckim/Sources/pyathena-TIGRESS/'}
sysconfig['tigress.princeton.edu']={'base':'/tigress/changgoo/','source':'/tigress/changgoo/pyathena-TIGRESS/'}
sysconfig['C02TW028HTDH']={'base':'/Users/ckim/','source':'/Users/ckim/Sources/pyathena-TIGRESS/'}

uname=subprocess.check_output(['uname','-a'])
sysname=uname.split(' ')[1]

if sysname in sysconfig:
    print '### setting up for %s system' % sysname
    sys.path.insert(0,sysconfig[sysname]['source'])
    base=sysconfig[sysname]['base']
    print '### base directory path is %s' % base
else:
    print '### no system is matched with %s' % sysname

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
