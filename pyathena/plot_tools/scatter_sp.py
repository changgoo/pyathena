import numpy as np
import matplotlib.pyplot as plt
from pyathena import read_starvtk,texteffect,set_units
def projection(sp,axis):
    if axis == 0 or axis == 'z':
        spx=sp['x1']
        spy=sp['x2']
        spz=sp['x3']
    elif axis == 1 or axis == 'y':
        spx=sp['x1']
        spy=sp['x3']
        spz=sp['x2']
    elif axis == 2 or axis == 'x':
        spx=sp['x2']
        spy=sp['x3']
        spz=sp['x1']
    return spx,spy,spz

def projection_v(sp,axis):
    if axis == 0 or axis == 'z':
        spx=sp['v1']
        spy=sp['v2']
        spz=sp['v3']
    elif axis == 1 or axis == 'y':
        spx=sp['v1']
        spy=sp['v3']
        spz=sp['v2']
    elif axis == 2 or axis == 'x':
        spx=sp['v2']
        spy=sp['v3']
        spz=sp['v1']
    return spx,spy,spz

def scatter_sp(sp,ax,axis=0,thickness=10,norm_factor=4., \
  type='slice',kpc=True,runaway=True,active=True):
    unit=set_units(muH=1.4271)
    Msun=unit['mass'].to('Msun').value
    Myr=unit['time'].to('Myr').value
    #print len(sp)
    if 'flag' in sp and active:
      sp=sp[sp['flag'] > -2]
    if len(sp) >0:
      runaways=(sp['mass'] == 0.0)
      sp_runaway=sp[runaways]
      sp_normal=sp[-runaways]
      #print len(sp_runaway)
      if len(sp_runaway) > 0 and runaway:
        spx,spy,spz=projection(sp_runaway,axis)
        spvx,spvy,spvz=projection_v(sp_runaway,axis)
        if kpc:
            spx = spx/1.e3
            spy = spy/1.e3
        if type == 'slice': 
            islab=np.where(abs(spz) < thickness)
            #ax.scatter(spx.iloc[islab],spy.iloc[islab],marker='.',c='k',alpha=1.0)
            #ax.quiver(spx.iloc[islab],spy.iloc[islab],
            #      spvx.iloc[islab],spvy.iloc[islab],color='k',alpha=1.0)
        ax.scatter(spx,spy,marker='o',color='k',alpha=1.0,s=10.0/norm_factor)
        #ax.quiver(spx,spy,spvx,spvy,color='w',alpha=0.5)
 
      if len(sp_normal) > 0: 
        spx,spy,spz=projection(sp_normal,axis)
        if kpc:
            spx = spx/1.e3
            spy = spy/1.e3
        if type == 'slice': xbool=abs(spz) < thickness
        spm=np.sqrt(sp_normal['mass']*Msun)/norm_factor
        spa=sp_normal['age']*Myr
        iyoung=np.where(spa < 40.)
        #print len(iyoung[0])
        if type == 'slice':
            islab=np.where(xbool*(spa<40))
            ax.scatter(spx.iloc[islab],spy.iloc[islab],marker='o',\
                s=spm.iloc[islab],c=spa.iloc[islab],\
                vmax=40,vmin=0,cmap=plt.cm.cool_r,alpha=1.0)
        ax.scatter(spx.iloc[iyoung],spy.iloc[iyoung],marker='o',\
            s=spm.iloc[iyoung],c=spa.iloc[iyoung],\
            vmax=40,vmin=0,cmap=plt.cm.cool_r,alpha=0.7)
