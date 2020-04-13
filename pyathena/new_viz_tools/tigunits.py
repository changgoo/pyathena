import pyathena as pa
import astropy.constants as ac
import astropy.units as au

coolftn=pa.coolftn()
units=pa.set_units(muH=1.4271)
to_Myr=units['time'].to('Myr').value
to_Msun=units['mass'].to('Msun').value
to_uG=units['magnetic_field'].value
to_pok=units['pressure'].cgs.value/ac.k_B.cgs.value
to_T1=(units['pressure']/units['density']*ac.m_p/ac.k_B).cgs.value

eta_conv=(ac.M_sun/ac.kpc**2/au.yr*au.km/au.s/ac.k_B).cgs.value
