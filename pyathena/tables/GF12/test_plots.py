from pyathena.set_plt import *
from .GF12 import *

def plot_ion_frac(tbl,ion_name):
    element=tbl.elements.loc[ion_name]
    ion_frac=tbl.ion_frac
    for i in range(element['number']+1):
        plt.loglog(ion_frac['temp'],ion_frac['{}{}'.format(ion_name,i)],
                   label='{}{}'.format(ion_name,int_to_roman(i+1)))
    plt.legend()
    plt.ylabel(r'{} $X_i$'.format(pt.elements.symbol(ion_name).name))
    plt.xlabel(r'Temperature [K]')
    
def plot_cie_cooling_ion(tbl,ion_name,ion_by_ion=True,noplot=False):
    element=tbl.elements.loc[ion_name]
    nion=element['number']+1
    A=element['abundance']
    total_cooling=tbl.get_cie_cooling(ion_name)
    if not noplot:
        l,=plt.loglog(tbl.temp,total_cooling*A,label=ion_name)
        if ion_by_ion and '{}0'.format(ion_name) in tbl.ion_frac:
            for i in range(nion):
                cie_cooling_ion=A*tbl.cie_cooling_per_ion[ion_name][:,i]*tbl.ion_frac['{}{}'.format(ion_name,i)]
                plt.loglog(tbl.temp,cie_cooling_ion,color=l.get_color(),ls=':')
            plt.ylabel(r'$\Lambda_{{CIE,{}}} [{{\rm erg s^{{-1}} cm^{{3}}}}]$'.format(ion_name))
            plt.ylim(total_cooling.min()*0.1*A)
    return total_cooling*A
    
def plot_cie_cooling_element(tbl,elements=['H','He','C','N','O','Ne','Mg','Si','S','Fe'],ion_by_ion=True):
    total_cooling=np.zeros_like(tbl.temp)
    for ion_name in elements:
        total_cooling+=plot_cie_cooling_ion(tbl,ion_name,ion_by_ion=ion_by_ion)
    plt.plot(tbl.temp,total_cooling,color='k',lw=3,alpha=0.5)
    if not ion_by_ion: plt.legend()
    plt.ylabel(r'$\Lambda_{{CIE}} [{{\rm erg s^{{-1}} cm^{{3}}}}]$')

if __name__ == '__main__':
    tbl=GF12_table()
    plot_ion_frac(tbl,'O')
    plt.show()
    plot_cie_cooling_ion(tbl,'O');
    plt.show()
    plot_cie_cooling_element(tbl,ion_by_ion=False)
    plt.ylim(1.e-26,1.e-21)
    plt.show()
    total_cooling=tbl.get_total_cie_cooling()
    plt.loglog(tbl.temp,total_cooling)
    total_cooling=tbl.get_total_cie_cooling(elements=tbl.elements.index)
    plt.loglog(tbl.temp,total_cooling)
    plt.show()
