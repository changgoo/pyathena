import numpy as np
import os,sys
import glob
import argparse

sys.path.insert(0,'../')
import pyathena.rst_handler as rh

if len(sys.argv) == 1:
    print "please specify filename"
else:
    filename=sys.argv[1]
    for sghost in [False,True]:
        fm,data=rh.read_rst_grid(filename,verbose=True,starghost=sghost)
