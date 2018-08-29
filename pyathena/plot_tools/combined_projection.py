import pyathena as pa
import matplotlib.pyplot as plt
from pyathena.plot_tools.scatter_sp import scatter_sp

def get_proj_ims(ds,slab_range=None,test=False):
    if slab_range is None:
        ns_start=ds.NGrids[2]/4+1
        ns_end=ds.NGrids[2]/4*3+1
    else:
        ns_start,ns_end = slab_range

    Lx=ds.domain['Lx'][1]
    x2min=ds.domain['left_edge'][1]-Lx
    x2max=ds.domain['right_edge'][1]+Lx

    x3min=ds.domain['right_edge'][2]
    x3max=ds.domain['left_edge'][2]
    for ns in range(ns_start,ns_end):
        g=ds._get_slab_grid(slab=ns,verbose=True)
        x3min=min(x3min,g[0]['left_edge'][2])
        x3max=max(x3max,g[0]['right_edge'][2])
    extent=(x2min,x2max,x3min,x3max)
    
    dlist=[]
    for ns in range(ns_start,ns_end):
        dlist.append(ds.read_all_data('density',slab=ns))
    den=np.concatenate(dlist)

    dlist=[]
    for ns in range(ns_start,ns_end):
        dlist.append(ds.read_all_data('T1',slab=ns))
    T1=np.concatenate(dlist)

    dlist=[]
    for ns in range(ns_start,ns_end):
        dlist.append(ds.read_all_data('magnetic_field',slab=ns))
    B=np.concatenate(dlist)

    Bmag=np.sqrt((B**2).sum(axis=-1))

    coolftn=pa.coolftn()
    temp=coolftn.get_temp(T1)

    Nz,Ny,Nx=den.shape

    unit=pa.set_units(muH=1.4271)
    unit['magnetic_field']

    dproj=den.mean(axis=2)
    Tproj=(den*temp).mean(axis=2)/dproj
    Bproj=(den*Bmag).mean(axis=2)*unit['magnetic_field'].value/dproj
    #dproj=np.roll(dproj,Nx/2,axis=1)
    #Tproj=np.roll(Tproj,Nx/2,axis=1)
    #Bproj=np.roll(Bproj,Nx/2,axis=1)
    #Tproj=(temp).mean(axis=2)
    dproj=np.concatenate([dproj,dproj,dproj],axis=1)
    Tproj=np.concatenate([Tproj,Tproj,Tproj],axis=1)
    Bproj=np.concatenate([Bproj,Bproj,Bproj],axis=1)

    return dproj,Tproj,Bproj,extent


def draw_merged_proj(fig,dproj,Tproj,Bproj,extent):
    Nz,Nx=dproj.shape
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False, aspect=1)

    drange=[1.e-4,100]
    Trange=[100,1.e6]
    Brange=[1.e-1,100]
    ix=np.arange(Nx)
    a1=(1-np.tanh((ix-1.0*Nx/3.)/float(0.1*Nx)))*0.5
    a3=(0.5*np.tanh((ix-2.0*Nx/3.)/float(0.1*Nx)))+0.5
    a2=1-(a1+a3)
    #plt.plot(ix,a1)
    #plt.plot(ix,a2)
    #plt.plot(ix,1-(a1+a2))
    # Create an alpha channel of linearly increasing values moving to the right.
    alpha1 = np.zeros(dproj.shape)
    alpha1=a1[np.newaxis,:]*np.ones(Nz)[:,np.newaxis]

    # Normalize the colors b/w 0 and 1, we'll then pass an MxNx4 array to imshow
    from matplotlib.colors import Normalize
    cmap=plt.cm.bone_r
    color1 = LogNorm(drange[0],drange[1], clip=True)(dproj)
    color1 = cmap(color1)
    color1[..., -1] = alpha1

    # Create an alpha channel of linearly increasing values moving to the right.
    alpha2 = np.ones(dproj.shape)
    alpha2=a2[np.newaxis,:]*np.ones(Nz)[:,np.newaxis]

    cmap=plt.cm.RdYlBu_r
    color2 = LogNorm(Trange[0],Trange[1], clip=True)(Tproj)
    color2 = cmap(color2)
    color2[..., -1] = alpha2

    # Create an alpha channel of linearly increasing values moving to the right.
    alpha3 = np.ones(dproj.shape)
    alpha3=a3[np.newaxis,:]*np.ones(Nz)[:,np.newaxis]

    cmap=plt.cm.cubehelix_r
    color3 = LogNorm(Brange[0],Brange[1], clip=True)(Bproj)
    color3 = cmap(color3)
    color3[..., -1] = alpha3

    # Create the figure and image
    # Note that the absolute values may be slightly different
    im1=ax.imshow(color1,origin='lower',extent=extent)
    im2=ax.imshow(color2,origin='lower',extent=extent)
    im3=ax.imshow(color3,origin='lower',extent=extent)

    # Create scale bar
    ax.plot([extent[0]+200,extent[0]+700],[extent[2]+350,extent[2]+350],lw=10,color='k')
    ax.text(extent[0]+350,extent[2]+400,'500 pc',ha='center',va='bottom',color='k',weight=300,**pa.texteffect(25))
    #ax.set_xlim(0,Nx-1)
    ax.set_axis_off()

    # Create colorbar
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cax1 = fig.add_axes([0.05, 0.03, 0.2, 0.02])
    cmap1 = plt.cm.bone_r
    norm1 = LogNorm(drange[0],drange[1])
    cb1 = mpl.colorbar.ColorbarBase(cax1, cmap=cmap1, norm=norm1, orientation='horizontal',
                                   ticks=[1.e-3,1.e-1,10.],)
    cb1.ax.set_title(r'Number Density [cm$^{-3}$]',fontsize=18)
    cb1.ax.tick_params(labelsize=20,direction='out') 

    cax2 = fig.add_axes([0.4, 0.03, 0.2, 0.02])
    cmap2 = plt.cm.RdYlBu_r
    norm2 = LogNorm(Trange[0],Trange[1])
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap=cmap2, norm=norm2, orientation='horizontal',)
    cb2.ax.set_title(r'Temperature [K]',fontsize=18)
    cb2.ax.tick_params(labelsize=20,direction='out') 

    cax3 = fig.add_axes([0.75, 0.03, 0.2, 0.02])
    cmap3 = plt.cm.cubehelix_r
    norm3 = LogNorm(Brange[0],Brange[1])
    cb3 = mpl.colorbar.ColorbarBase(cax3, cmap=cmap3, norm=norm3, orientation='horizontal',)
    cb3.ax.set_title(r'Magnetic Field Strenghth [$\mu $G]',fontsize=18)
    cb3.ax.tick_params(labelsize=20,direction='out') 
    
def run_R8(base,pid):
    for i in range(700):
        print("{} of {}...".format(i,700))
        ds=pa.AthenaDataSet('{}/{}/id0/{}.{:04d}.vtk'.format(base,pid,pid,i))
        dproj,Tproj,Bproj=get_proj_ims(ds,slab_range=[9,21])
        Nz,Nx=dproj.shape
        ratio=Nz/float(Nx)
        plt.style.use('default')
        fig=plt.figure(0,figsize=(15,15*ratio))
        draw_merged_proj(fig,dproj,Tproj,Bproj)
        fig.savefig('{}{}/proj_figures/{}.merged_proj.{:04d}.png'.format(base,pid,pid,i),dpi=100)
        fig.clf()
