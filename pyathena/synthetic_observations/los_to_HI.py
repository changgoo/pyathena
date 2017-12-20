import numpy as np

def los_to_HI(dens,temp,vel,vchannel,deltas=1.,WF=False):
    """
        inputs:
            dens: number density of hydrogen in units of 1/cm^3
            temp: temperature in units of K
            vel: line-of-sight velocity in units of km/s
        parameters:
            vmax: maximum range of velocity channel (+- vmax) in units of km/s
            dvch: velocity channel resolution in units of km/s
            deltas: length of line segments in units of pc
        outputs: a dictionary
            TB: the brightness temperature
            tau: optical depth
            vchannel: velocity channel
            TBthin: the brithgtness temperature assuming optically-thin case.
    """
    
    Tlos=temp[:,:,np.newaxis]
    vlos=vel[:,:,np.newaxis]
    nlos=dens[:,:,np.newaxis]

    if WF: Tspin=Tspin_WF(Tlos,nlos)
    else: Tspin=Tlos
    
    ds=deltas*3.085677581467192e+18

    v_L=0.21394414*np.sqrt(Tlos) # in units of km/s
    phi_v=0.00019827867/v_L*np.exp(-(1.6651092223153954*
          (vchannel[np.newaxis,np.newaxis,:]-vlos)/v_L)**2) # time
    kappa_v=2.6137475e-15*nlos/Tspin*phi_v # area/volume = 1/length
    tau_los=kappa_v*ds # dimensionless

    tau_cumul=tau_los.cumsum(axis=1)
    TB=np.nansum(Tspin*(1-np.exp(-tau_los))*np.exp(-tau_cumul),axis=1) # same unit with Tspin
    TBthin=np.nansum(Tspin*tau_los,axis=1) # same unit with Tspin
    tau_v=np.nansum(kappa_v*ds,axis=1) # dimensionless
    return TB,tau_v,TBthin

def k10h(T2):
    """
       input: T/100K
       output: collisional excitation rate of hydrogen in c.g.s.
    """
    k10_1=1.19e-10*T2**(0.74-0.20*np.log(T2))
    k10_2=2.24e-10*T2**0.207*np.exp(-0.876/T2)
    k10=k10_1
    idx=k10_2 > k10_1
    k10[idx]=k10_2[idx]
    return k10

def k10e(T2):
    """
       input: T/100K
       output: collisional excitation rate of electrion in c.g.s.
    """
    Temp=T2*1.e2
    k10=-9.607+0.5*np.log10(Temp)*np.exp(-(np.log10(Temp))**4.5/1.8e3)
    return 10.**k10

def Tspin_WF(Temp,nH,nalpha=1.e-6):
    """
       spin temperature including the WF effect for a given nalpha
    """
    T2=Temp/(1.e2)
    Tk=Temp
    TL=Temp
    TA=3.77
    T21=0.068168759
    A10=2.8843e-15

    k10=k10h(T2)

    yc=(T21*k10*nH/A10/Temp)
    yalpha=(7.76e11*nalpha/TL/np.sqrt(Tk))
    Ts=(TA+yc*Temp+yalpha*TL)/(1.0+yc+yalpha)

    return Ts


