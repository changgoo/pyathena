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

def _turb_pok(field,data):
        return data["gas","density"]*data["gas","velocity_magnitude"]**2/kboltz

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
def _metallicity(field,data):
        return data["athena","specific_scalar[0]"]

def _metal(field,data):
        return data["athena","specific_scalar[0]"]*data["gas","density"]

def _metal_cl(field,data):
        return data["athena","specific_scalar[1]"]*data["gas","density"]

def _metal_run(field,data):
        return data["athena","specific_scalar[2]"]*data["gas","density"]

def _radius(field, data):
    return np.sqrt(data['x']**2+data['y']**2+data['z']**2)

def _velocity_r(field, data):
    return data['velocity_x']*data['x']/data['radius']+\
           data['velocity_y']*data['y']/data['radius']+\
           data['velocity_z']*data['z']/data['radius']

def _momentum_r(field, data):
    return data["cell_mass"]*data["velocity_r"]

def _total_kinetic_energy(field, data):
    return 0.5*data["cell_mass"]*data["velocity_magnitude"]**2

def _total_magnetic_energy(field, data):
    return data["magnetic_field_magnitude"]**2*data["cell_volume"]/8.0/np.pi

unit_base={"length_unit": (1.0,"pc"), 
           "time_unit": (1.0,"s*pc/km"), 
           "mass_unit": (2.38858753789e-24,"g/cm**3*pc**3"), 
           "velocity_unit": (1.0,"km/s"),
           "magnetic_unit": (5.4786746797e-07,"gauss")}

def get_scalars(ds):
    scal_fields=[]
    for f in ds.field_list: 
        code,field=f
        if field.startswith('specific_scalar'):
            scal_fields.append(field)

    return scal_fields


def add_yt_fields(ds,cooling=True,mhd=True,rotation=True):
    ds.add_field(("gas","nH"),function=_ndensity,sampling_type='cell', \
      units='cm**(-3)',display_name=r'$n_{\rm H}$')
    ds.add_field(("gas","ram_pok_z"),function=_ram_pok_z,sampling_type='cell', \
      units='K*cm**(-3)',display_name=r'$P_{\rm turb,z}/k_{\rm B}$')
    ds.add_field(("gas","turb_pok"),function=_turb_pok,sampling_type='cell', \
      units='K*cm**(-3)',display_name=r'$P_{\rm turb}/k_{\rm B}$')
    ds.add_field(("gas","radius"), function=_radius, \
      sampling_type='cell',units='pc', \
      display_name=r'$r$',force_override=True)
    ds.add_field(("gas","total_kinetic_energy"), function=_total_kinetic_energy, \
      sampling_type='cell',units='erg', \
      display_name=r'$E_{\rm kin}$',force_override=True)
    ds.add_field(("gas","velocity_r"), function=_velocity_r, \
      sampling_type='cell',units='km/s', display_name=r'$v_r$')
    ds.add_field(("gas","momentum_r"), function=_momentum_r, \
      sampling_type='cell',units='Msun*km/s', display_name=r'$p_r$')
    if cooling:
        ds.add_field(("gas","pok"),function=_pok,sampling_type='cell', \
          units='K*cm**(-3)',display_name=r'$P/k_{\rm B}$')
        ds.add_field(("gas","cs"),function=_cs,sampling_type='cell', \
          units='km*s**(-1)',display_name=r'$c_s$')
        ds.add_field(("gas","T1"),function=_T1,sampling_type='cell', \
          units='K',display_name=r'$T_1$')
        ds.add_field(("gas","mu"),function=_mu,sampling_type='cell', \
          units='',display_name=r'$\mu$',force_override=True)
        ds.add_field(("gas","temperature"),function=_temperature,sampling_type='cell', \
          units='K',display_name=r'$T$',force_override=True)
    if rotation:
        ds.add_field(("gas","dvelocity_y"),function=_dvelocity,sampling_type='cell', \
          units='km/s',display_name=r'$\delta v_y$',force_override=True)
        ds.add_field(("gas","dvelocity_magnitude"),function=_dvelocity_mag,sampling_type='cell', \
          units='km/s',display_name=r'$v$',force_override=True)
        ds.add_field(("gas","dkinetic_energy"),function=_dkinetic_energy,sampling_type='cell', \
          units='erg/cm**3',display_name=r'$E_k$',force_override=True)
    if mhd:
        ds.add_field(("gas","mag_pok"),function=_mag_pok,sampling_type='cell', \
          units='K*cm**(-3)',display_name=r'$P_{\rm mag}/k_{\rm B}$')
        ds.add_field(("gas","total_magnetic_energy"), function=_total_magnetic_energy, \
          sampling_type='cell', \
          units='erg', display_name=r'$E_{\rm mag}$',force_override=True)

    scal_fields=get_scalars(ds)
    if len(scal_fields)>0:
        ds.add_field(("gas","metallicity"),function=_metallicity,sampling_type='cell', \
          units='dimensionless',display_name=r'$Z$')
        ds.add_field(("gas","metal0"),function=_metal,sampling_type='cell', \
          units='g*cm**(-3)',display_name=r'$\rho_{\rm metal}$')
    if len(scal_fields)>1:
        ds.add_field(("gas","metal1"),function=_metal_cl,sampling_type='cell', \
          units='g*cm**(-3)',display_name=r'$\rho_{\rm metal,cl}$')
    if len(scal_fields)>2:
        ds.add_field(("gas","metal2"),function=_metal_run,sampling_type='cell', \
          units='g*cm**(-3)',display_name=r'$\rho_{\rm metal,run}$')
