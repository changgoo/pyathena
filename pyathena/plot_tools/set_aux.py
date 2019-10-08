from matplotlib import cm
from matplotlib.colors import LogNorm,Normalize,SymLogNorm
from .shiftedColorMap import *

def set_aux(model='solar',verbose=False):
    aux={}
    aux['density']=dict(label=r'$n_H\;[{\rm cm}^{-3}]$', \
        unit='cm**(-3)', limits=(1.e-6,1.e6), \
        cmap=cm.Spectral_r,clim=(2.e-5,2.e2), \
        cticks=(1.e-4,1.e-2,1,1.e2), \
        n_bins=128, norm=LogNorm())
    aux['nH']=dict(label=r'$n_H\;[{\rm cm}^{-3}]$', \
        unit='cm**(-3)', limits=(1.e-6,1.e6), \
        cmap=cm.Spectral_r,clim=(2.e-5,2.e2), \
        cticks=(1.e-4,1.e-2,1,1.e2), \
        n_bins=128, norm=LogNorm())
    aux['pok']=dict(label=r'$P/k_B\;[{\rm K}\,{\rm cm}^{-3}]$', \
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        cmap=cm.plasma,clim=(10,5.e5), \
        cticks=(1.e2,1.e3,1.e4,1.e5), \
        n_bins=128, norm=LogNorm())
    aux['temperature']=dict(label=r'$T\;[{\rm K}]$', \
        unit='K', limits=(1.e0,1.e9), \
        cmap=shiftedColorMap(cm.RdYlBu_r,midpoint=3/7.), \
        clim=(10,1.e8), \
        cticks=(1.e2,1.e4,1.e6,1.e8), \
        n_bins=128, norm=LogNorm())
    aux['surface_density']=dict( \
        label=r'$\Sigma\;[{\rm M}_{\odot} {\rm pc}^{-2}]$', \
        cmap=cm.pink_r,clim=(0.1,100),norm=LogNorm())
    aux['velocity_z']=dict(label=r'$v_z\;[{\rm km/s}]$', \
        unit='km/s', limits=(-1500,1500), \
        cmap=cm.RdBu_r,clim=(-200,200), \
        cticks=(-100,0,100), \
        n_bins=256, norm=Normalize())
    aux['velocity_x']=dict(label=r'$v_x\;[{\rm km/s}]$', \
        unit='km/s', limits=(-1500,1500), \
        cmap=cm.RdBu_r,clim=(-200,200), \
        cticks=(-100,0,100), \
        n_bins=256, norm=Normalize())
    aux['velocity_y']=dict(label=r'$v_y\;[{\rm km/s}]$', \
        unit='km/s', limits=(-1500,1500), \
        cmap=cm.RdBu_r,clim=(-200,200), \
        cticks=(-100,0,100), \
        n_bins=256, norm=Normalize())
    aux['magnetic_field_strength']=dict(label=r'$B\;[\mu{\rm G}]$', \
        unit='uG', \
        cmap=cm.viridis,clim=(0.01,10),factor=1, \
        n_bins=128, norm=LogNorm())
    aux['mag_pok']=dict(label=r'$P_{\rm mag}/k_B\;[{\rm K}{\rm cm}^{-3}]$',\
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        cmap=cm.plasma,clim=(10,5.e5), \
        n_bins=128, norm=LogNorm())
    aux['ram_pok_z']=dict(\
        label=r'$P_{\rm turb}/k_B\;[{\rm K}{\rm cm}^{-3}]$', \
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        cmap=cm.plasma,clim=(10,5.e5), \
        n_bins=128, norm=LogNorm())
    aux['plasma_beta']=dict(label=r'$\beta$', limits=(1.e-4,1.e16), \
        n_bins=256, norm=LogNorm())
    aux['specific_scalar3']=dict(label=r'$C_{\rm icm}$', limits=(-0.1,1.1), clim=(-0.1,1.10),\
        cticks=(0,0.5,1),n_bins=256, norm=Normalize(), cmap=cm.Reds)
    aux['specific_scalar3_proj']=dict(label=r'$C_{\rm icm}$', limits=(-0.1,1.1), clim=(-0.1,1.10),\
        cticks=(0,0.5,1),n_bins=256, norm=Normalize(), cmap=cm.Reds)

    aux['star_particles']=dict(label=r'${\rm age [Myr]}$', \
        unit='Myr', limits=(0,40), \
        cmap=cm.cool_r,clim=(0,40), \
        cticks=(0,20,40), \
        n_bins=256, norm=LogNorm())

    if model.startswith('R4') or model.startswith('LGR4'):
        if verbose: print('auxilary information is set for R4')
        aux['nH']['clim']=(2.e-4,2.e3)
        aux['nH']['cticks']=(1.e-2,1,1.e2)
        aux['pok']['clim']=(10,5.e6)
        aux['pok']['cticks']=(1.e2,1.e4,1.e6)
        aux['surface_density']['clim']=(1,1000)
        aux['velocity_z']['clim']=(-300,300)
        aux['velocity_z']['cticks']=(-200,0,200)
        aux['magnetic_field_strength']['clim']=(0.01,100)

    if model.startswith('R2') or model.startswith('LGR2'):
        if verbose: print('auxilary information is set for R2')
        aux['nH']['clim']=(2.e-4,2.e3)
        aux['nH']['cticks']=(1.e-2,1,1.e2)
        aux['pok']['clim']=(10,5.e6)
        aux['pok']['cticks']=(1.e2,1.e4,1.e6)
        aux['velocity_z']['clim']=(-300,300)
        aux['velocity_z']['cticks']=(-200,0,200)
        aux['surface_density']['clim']=(1,1000)
        aux['magnetic_field_strength']['clim']=(0.01,100)

    if model is 'multi_SN':
        aux['nH']['clim']=(2.e-5,2.e2)
        aux['pok']['clim']=(50,1.e5)

    return aux
