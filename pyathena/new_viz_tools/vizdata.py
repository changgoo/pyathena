import pyathena as pa
import pandas as pd
import numpy as np
import glob,os

from .plotters import *

class data(object):
    def __init__(self,pid,itime,base='/tigress/changgoo/',pdir=None,vtkdir=None,Omega=0.028):

        if pdir is None:
            pdir=pid
        if vtkdir is None:
            vtkdir=base+pdir
        self.pid=pid
        self.pdir=pdir+'/'
        self.base=base
        self.vtkdir=vtkdir
        self.itime=itime
        self.par=pa.get_params('{}{}/{}.par'.format(base,pdir,pid))
        if 'Omega' in self.par:
            self.Omega=self.par['Omega']
        else:
            self.Omega=Omega
        self.get_list_itimes()

    def get_list_itimes(self):
        vtkfiles=glob.glob('{}{}/surf/{}.????.surf.p'.format(self.base,self.pdir,self.pid))
        vtkfiles.sort() 
        self.itimes = []
        for f in vtkfiles:
            fbase=os.path.basename(f)
            self.itimes.append(int(fbase.split('.')[-3]))

    def readhst(self):
        pid=self.pid
        pdir=self.pdir
        base=self.base
        itime=self.itime
        self.hzp=pd.read_pickle('{}{}/hst/{}.hst_zp.p'.format(base,pdir,pid))
        self.sn=pd.read_pickle('{}{}/hst/{}.sn.p'.format(base,pdir,pid))

    def readall(self):
        pid=self.pid
        pdir=self.pdir
        base=self.base
        itime=self.itime
        
        self.readhst()
        #self.ds=pa.AthenaDataSet('{}/id0/{}.{:04d}.vtk'.format(self.vtkdir,pid,itime))
        self.sp=pa.read_starvtk('{}{}/starpar/{}.{:04d}.starpar.vtk'.format(base,pdir,pid,itime))
        self.slc=pd.read_pickle('{}{}/slice/{}.{:04d}.slice.p'.format(base,pdir,pid,itime))
        self.surf=pd.read_pickle('{}{}/surf/{}.{:04d}.surf.p'.format(base,pdir,pid,itime))
        for proj in ['ddproj','dTproj','Bproj']:
            proj_file='{}{}/proj/{}.{:04d}.{}.p'.format(base,pdir,pid,itime,proj)
            if os.path.isfile(proj_file):
                if proj == 'Bproj':
                    setattr(self,proj,pd.read_pickle(proj_file))
                else:
                    setattr(self,proj[1:],pd.read_pickle(proj_file))
        
        
        sp=self.sp
        if len(sp)>0:
            if not 'mage' in sp: sp['mage']=sp['age']
            icluster = (sp['mass']>0) & (sp['mage']*to_Myr <40)
            cl=sp[icluster]
            runaway=sp[sp['mass']==0]
            self.clpos=dict(x=cl['x1'],y=cl['x2'],z=cl['x3'])
            self.clmass=np.cbrt(cl['mass']*to_Msun)*5
            self.clage=cl['mage']*to_Myr
            self.runpos=dict(x=runaway['x1'],y=runaway['x2'],z=runaway['x3'])
        self.clmass_keys=np.cbrt([1.e3,1.e4,1.e5])*5
        self.clcmap=add_alpha_to_cmap(cma.ember_r,expo=1,reversed=True)
        
        time = self.slc['time']
        Lx = self.par['x1max']-self.par['x1min']
        Ly = self.par['x2max']-self.par['x2min']
        dy = Lx/self.par['Nx2']
        self.yoffset=-np.mod(time*Lx,Ly)
        self.joffset=int(self.yoffset/dy)
