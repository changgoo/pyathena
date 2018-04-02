from yt import add_field
from yt import YTQuantity
from yt.utilities.physical_constants import \
    mh, \
    me, \
    sigma_thompson, \
    clight, \
    kboltz, \
    G
import pickle as p
import numpy as np
import pyathena as pa

# basic qunatities with renormalization
def _ndensity(field, data):
        return data["gas","density"]/(1.4271*mh)

def _ram_pok_z(field,data):
        return data["gas","density"]*data["gas","velocity_z"]**2/kboltz

# thermodynamics quantities
def _pok(field, data):
        return data["gas","pressure"]/kboltz

def _cs(field, data):
        return np.sqrt(data["gas","pressure"]/data["gas","density"])

def _T1(field, data):
        return data["gas","pressure"]/data["gas","density"]*mh/kboltz

def _mu(field, data):
        cf=pa.coolftn()
        T1=data["gas","T1"].d
        temp=cf.get_temp(T1)
        return temp/T1

def _temperature(field,data):
        return data["gas","T1"]*data["gas","mu"]

# rotation
Omega=YTQuantity(28,"km/s/kpc")
def _dvelocity(field,data):
        return data["gas","velocity_y"]+data["gas","x"]*Omega

def _dvelocity_mag(field,data):
        return np.sqrt(data["gas","velocity_x"]**2+data["gas","dvelocity_y"]**2+data["gas","velocity_z"]**2)

def _dkinetic_energy(field,data):
    return 0.5*data['gas','dvelocity_magnitude']**2*data['gas','density']

# magnetic fields
def _mag_pok(field,data):
        return data["gas","magnetic_pressure"]/kboltz

# metals
def _metal(field,data):
        return data["athena","specific_scalar[0]"]*data["gas","density"]

def _metal_cl(field,data):
        return data["athena","specific_scalar[1]"]*data["gas","density"]

def _metal_run(field,data):
        return data["athena","specific_scalar[2]"]*data["gas","density"]

unit_base={"length_unit": (1.0,"pc"), 
           "time_unit": (1.0,"s*pc/km"), 
           "mass_unit": (2.38858753789e-24,"g/cm**3*pc**3"), 
           "velocity_unit": (1.0,"km/s"),
           "magnetic_unit": (5.4786746797e-07,"gauss")}

import matplotlib.pyplot as plt
from .shiftedColorMap import *

def get_scalars(ds):
    scal_fields=[]
    for f in ds.field_list: 
        code,field=f
        if field.startswith('specific_scalar'):
            scal_fields.append(field)

    return scal_fields


def add_yt_fields(ds,cooling=True,mhd=True,rotation=True):
    ds.add_field(("gas","nH"),function=_ndensity, \
      units='cm**(-3)',display_name=r'$n_{\rm H}$')
    ds.add_field(("gas","ram_pok_z"),function=_ram_pok_z, \
      units='K*cm**(-3)',display_name=r'$P_{\rm turb}/k_{\rm B}$')
    if cooling:
        ds.add_field(("gas","pok"),function=_pok, \
          units='K*cm**(-3)',display_name=r'$P/k_{\rm B}$')
        ds.add_field(("gas","cs"),function=_cs, \
          units='km*s**(-1)',display_name=r'$c_s$')
        ds.add_field(("gas","T1"),function=_T1, \
          units='K',display_name=r'$T_1$')
        ds.add_field(("gas","mu"),function=_mu, \
          units='',display_name=r'$\mu$',force_override=True)
        ds.add_field(("gas","temperature"),function=_temperature, \
          units='K',display_name=r'$T$',force_override=True)
    if rotation:
        ds.add_field(("gas","dvelocity_y"),function=_dvelocity, \
          units='km/s',display_name=r'$\delta v_y$',force_override=True)
        ds.add_field(("gas","dvelocity_magnitude"),function=_dvelocity_mag, \
          units='km/s',display_name=r'$v$',force_override=True)
        ds.add_field(("gas","dkinetic_energy"),function=_dkinetic_energy, \
          units='erg/cm**3',display_name=r'$E_k$',force_override=True)
    if mhd:
        ds.add_field(("gas","mag_pok"),function=_mag_pok, \
          units='K*cm**(-3)',display_name=r'$P_{\rm mag}/k_{\rm B}$')
    scal_fields=get_scalars(ds)
    if len(scal_fields)>0:
        ds.add_field(("gas","metal0"),function=_metal, \
          units='g*cm**(-3)',display_name=r'$\rho_{\rm metal}$')
    if len(scal_fields)>1:
        ds.add_field(("gas","metal1"),function=_metal_cl, \
          units='g*cm**(-3)',display_name=r'$\rho_{\rm metal,cl}$')
    if len(scal_fields)>2:
        ds.add_field(("gas","metal2"),function=_metal_run, \
          units='g*cm**(-3)',display_name=r'$\rho_{\rm metal,run}$')

def set_aux(model='solar',verbose=False):
    aux={}
    aux['density']=dict(label=r'$n_H\;[{\rm cm}^{-3}]$', \
        unit='cm**(-3)', limits=(1.e-6,1.e6), \
        cmap=plt.cm.Spectral_r,clim=(2.e-5,2.e2), \
        cticks=(1.e-4,1.e-2,1,1.e2), \
        n_bins=128, log=True)
    aux['nH']=dict(label=r'$n_H\;[{\rm cm}^{-3}]$', \
        unit='cm**(-3)', limits=(1.e-6,1.e6), \
        cmap=plt.cm.Spectral_r,clim=(2.e-5,2.e2), \
        cticks=(1.e-4,1.e-2,1,1.e2), \
        n_bins=128, log=True)
    aux['pok']=dict(label=r'$P/k_B\;[{\rm K}\,{\rm cm}^{-3}]$', \
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        cmap=plt.cm.gnuplot2,clim=(10,5.e5), \
        cticks=(1.e2,1.e3,1.e4,1.e5), \
        n_bins=128, log=True)
    aux['temperature']=dict(label=r'$T\;[{\rm K}]$', \
        unit='K', limits=(1.e0,1.e9), \
        cmap=shiftedColorMap(plt.cm.RdYlBu_r,midpoint=3/7.), \
        clim=(10,1.e8), \
        cticks=(1.e2,1.e4,1.e6,1.e8), \
        n_bins=128, log=True)
    aux['surface_density']=dict( \
        label=r'$\Sigma\;[{\rm M}_{\odot} {\rm pc}^{-2}]$', \
        cmap=plt.cm.pink_r,clim=(0.1,100),log=True)
    aux['dvelocity_magnitude']=dict(label=r'$v\;[{\rm km/s}]$', \
        unit='km/s', limits=(0.1,1.e4), \
        cmap=plt.cm.jet,clim=(1,1000), \
        n_bins=128, log=True)
    aux['velocity_z']=dict(label=r'$v_z\;[{\rm km/s}]$', \
        unit='km/s', limits=(-1500,1500), \
        cmap=plt.cm.RdBu_r,clim=(-200,200), \
        cticks=(-100,0,100), \
        n_bins=256, log=False)
    aux['magnetic_field_strength']=dict(label=r'$B\;[\mu{\rm G}]$', \
        unit='uG', \
        cmap=plt.cm.viridis,clim=(0.01,10),factor=1, \
        n_bins=128, log=True)
    aux['mag_pok']=dict(label=r'$P_{\rm mag}/k_B\;[{\rm K}{\rm cm}^{-3}]$',\
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        cmap=plt.cm.gnuplot2,clim=(10,5.e5), \
        n_bins=128, log=True)
    aux['ram_pok_z']=dict(\
        label=r'$P_{\rm turb}/k_B\;[{\rm K}{\rm cm}^{-3}]$', \
        unit='K*cm**(-3)', limits=(1.e-2,1.e8), \
        cmap=plt.cm.gnuplot2,clim=(10,5.e5), \
        n_bins=128, log=True)
    aux['plasma_beta']=dict(label=r'$\beta$', limits=(1.e-4,1.e16), \
        n_bins=256, log=True)
    aux['star_particles']=dict(label=r'${\rm age [Myr]}$', \
        unit='Myr', limits=(0,40), \
        cmap=plt.cm.cool_r,clim=(0,40), \
        cticks=(0,20,40), \
        n_bins=256, log=False)
    aux['specific_scalar[0]']=dict(label=r'$Z$', limits=(0,2), \
        n_bins=256, log=False)
    aux['specific_scalar[1]']=dict(label=r'$Z_{\rm cl}$', limits=(0,2), \
        n_bins=256, log=False)
    aux['specific_scalar[2]']=dict(label=r'$Z_{\rm run}$', limits=(0,2), \
        n_bins=256, log=False)
    aux['specific_scalar[3]']=dict(label=r'$Z_{\rm Ia}$', limits=(0,2), \
        n_bins=256, log=False)
    aux['specific_scalar[4]']=dict(label=r'$C_{\rm ICM}$', limits=(0,2), \
        n_bins=256, log=False)

    if model.startswith('R4'):
        if verbose: print('auxilary information is set for R4')
        aux['nH']['clim']=(2.e-4,2.e3)
        aux['nH']['cticks']=(1.e-2,1,1.e2)
        aux['pok']['clim']=(10,5.e6)
        aux['pok']['cticks']=(1.e2,1.e4,1.e6)
        aux['surface_density']['clim']=(1,1000)
        aux['velocity_z']['clim']=(-300,300)
        aux['velocity_z']['cticks']=(-200,0,200)
        aux['magnetic_field_strength']['clim']=(0.01,100)


    elif model.startswith('R2'):
        if verbose: print('auxilary information is set for R2')
        aux['nH']['clim']=(2.e-4,2.e3)
        aux['nH']['cticks']=(1.e-2,1,1.e2)
        aux['pok']['clim']=(10,5.e6)
        aux['pok']['cticks']=(1.e2,1.e4,1.e6)
        aux['velocity_z']['clim']=(-300,300)
        aux['velocity_z']['cticks']=(-200,0,200)
        aux['surface_density']['clim']=(1,1000)
        aux['magnetic_field_strength']['clim']=(0.01,100)

    elif model is 'multi_SN':
        aux['nH']['clim']=(2.e-5,2.e2)
        aux['pok']['clim']=(50,1.e5)
    else:
        if verbose: print('auxilary information is set for Solar nbhd.')
    return aux

def check_aux(fields):
    aux=set_aux()
    for f in fields:
        if f not in aux:
            print("auxiliary information for %s is missing",f)
            print(aux[f])
