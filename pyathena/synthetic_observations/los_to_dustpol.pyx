# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True
import numpy as np
cimport numpy as np # access to Numpy from Cython layer
from libc.math cimport exp 


cdef void calc_stokes(double[::1] nH, double[::1] Bx, double[::1] By, double[::1] Bz, double[::1] pol,
                         int Nlos,double deltas, double Bnu, double sigma, double p0,int attenuation):
    """
    inputs: 
        nH in units of cm^-3
        Bfield (X,Y,Z) following the Healpix convention
    parameters:
        Bnu(T_dust,nu0): planck function (default T_dust=18K, nu0=353GHz; Planck XX 2015)
        sigma: absorption crossection in units of cm^-2 at 353GHz (default 1.2e-26 cm^-2; Planck XX 2015)
        deltas: length of lins segments in units of pc
        p0: intrinsic polarization fraction (default 20%; Planck XX 2015)
        attenuation: if true, self-attenuation is considered (default False)
    output:
        I, Q, U stokes' parameters in units of the input Bnu units. (default is MJy/sr)
    """
    cdef:
        double ds
        int i
        double Bperp2, B2, cos2phi, sin2phi, cosgam2
        double dtau, I, Q, U,tau
    ds=deltas*3.085677581467192e+18
    I=0.
    Q=0.
    U=0.
    
    for i in range(Nlos):
        Bperp2=Bx[i]*Bx[i]+By[i]*By[i]
        B2=Bperp2+Bz[i]*Bz[i]
        cos2phi=(Bx[i]**2-By[i]**2)/Bperp2
        sin2phi=Bx[i]*By[i]/Bperp2
        cosgam2=Bperp2/B2

        dtau=sigma*nH[i]*ds
     #   print nH.sum()*ds.cgs,tau[-1], np.exp(-tau[-1])

        I+=Bnu*(1.0-p0*(cosgam2-2./3.0))*dtau#*np.exp(-tau)
        Q+=p0*Bnu*cos2phi*cosgam2*dtau#*np.exp(-tau)
        U+=p0*Bnu*sin2phi*cosgam2*dtau#*np.exp(-tau)
        if attenuation == 1:
            tau=0
            for j in range(i):
                tau+=dtau
            I=I*exp(-tau)
            Q=Q*exp(-tau)
            U=U*exp(-tau)
    pol[0]=I
    pol[1]=Q
    pol[2]=U

cpdef los_to_dustpol(nH, Bx, By, Bz, deltas, args):
    cdef:
        double[::1] _nH= np.array(nH, np.float64)
        double[::1] _Bx= np.array(Bx, np.float64)
        double[::1] _By= np.array(By, np.float64)
        double[::1] _Bz= np.array(Bz, np.float64)
        double[::1] _pol=np.empty(3, np.float64)
        int Nlos = len(nH)
        double ds = float(deltas)
        double Bnu = float(args['Bnu'])
        double sigma = float(args['sigma'])
        double p0 = float(args['p0'])
        int attenuation = int(args['attenuation'])
    
    calc_stokes(_nH, _Bx, _By, _Bz, _pol, Nlos, ds, Bnu, sigma, p0, attenuation)
    
    return np.array(_pol)
