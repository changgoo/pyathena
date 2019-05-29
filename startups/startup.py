import sys,os,glob
import subprocess

sysconfig={}
sysconfig['ccalin004']={'base':'/mnt/ceph/users/ckim/','source':'/mnt/home/ckim/Sources/pyathena-TIGRESS/'}
sysconfig['tigressdata.princeton.edu']={'base':'/tigress/changgoo/','source':'/tigress/changgoo/pyathena-TIGRESS/'}
sysconfig['tigressdata2.princeton.edu']={'base':'/tigress/changgoo/','source':'/tigress/changgoo/pyathena-TIGRESS/'}
sysconfig['cori']={'base':'/global/cscratch1/sd/changgoo/','source':'/global/u2/c/changgoo/pyathena/'}
sysconfig['princeton-macbook']={'base':'/Users/ckim/Research/TIGRESS/','source':'/Users/ckim/Sources/pyathena-TIGRESS/'}

sysname=os.uname()[1]

if sysname in sysconfig:
    print('### setting up for %s system' % sysname)
    #sys.path.insert(0,sysconfig[sysname]['source'])
    base=sysconfig[sysname]['base']
    sourcedir=sysconfig[sysname]['source']
    print('### base directory path is %s' % base)
else:
    print('### no system is matched with %s' % sysname)

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LogNorm,Normalize,SymLogNorm
import cmocean

import pandas as pd
import xarray as xr

from IPython import get_ipython
ipython = get_ipython()

# If in ipython, load autoreload extension
if 'ipython' in globals():
    print('\nWelcome to IPython!\n')
    ipython.magic('load_ext autoreload')
    ipython.magic('autoreload 2')

# Display all cell outputs in notebook
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = 'all'
