import yt
import pyathena.yt_analysis.ytathena as ya
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import glob,os

class phase_parameters(object):
    def __init__(self):
        ''' 
        set global parameters

        (used in phase_part function)
        
        bin_fields: list 
            x and y field names.
            yt field names or ytathena added field names
            [[x_field1, y_field1],[x_field2, y_field2],...]
                    
        extrema: dict (min, max)
            extrema for the field.
        
        logs: dict (bool)
            if true or missing field name, log-scaled
        
        (used here)
        
        tigress_unit_system: dict (yt.UnitSystem)
        
        '''
        tigress_unit_system=yt.UnitSystem('tigress','pc','Msun','Myr',)
        tigress_unit_system['velocity']='km/s'
        tigress_unit_system['magnetic_field']='uG'
 
        bin_fields=[]
 
        bin_fields.append(['nH','pok'])
        #bin_fields.append(['nH','temperature'])
        bin_fields.append(['velocity_magnitude','sound_speed'])
        bin_fields.append(['nH','velocity_divergence'])
        bin_fields.append(['temperature','metallicity'])
        #bin_fields.append(['nH','velocity_magnitude'])
        #bin_fields.append(['radius','velocity_magnitude'])
        #bin_fields.append(['radius','temperature'])
        #bin_fields.append(['radius','nH'])
        #bin_fields.append(['nH','magnetic_field_magnitude'])
        #bin_fields.append(['temperature','magnetic_field_magnitude'])
        #bin_fields.append(['sound_speed','alfven_speed'])
        #bin_fields.append(['pok','ram_pok_z'])
        #bin_fields.append(['pok','mag_pok'])
        #bin_fields.append(['total_kinetic_energy','total_magnetic_energy'])
 
        extrema={}
        extrema['nH']=(1.e-5,1.e4)
        extrema['pok']=(1,1.e8)
        extrema['temperature']=(1,1.e9)
        extrema['velocity_magnitude']=(1.e-1,1.e4)
        extrema['sound_speed']=(0.1,1.e4)
        extrema['alfven_speed']=(1.e-1,1.e4)
        extrema['mag_pok']=(1,1.e8)
        extrema['ram_pok_z']=(1,1.e8)
        extrema['magnetic_field_magnitude']=(1.e-4,1.e2)
        extrema['total_kinetic_energy']=(1.e38,1.e48)
        extrema['total_magnetic_energy']=(1.e38,1.e48)
        extrema['velocity_divergence']=(-30,30)
        extrema['metallicity']=(0.02,0.2)
 
        logs={}
        logs['radius']=False
        logs['velocity_divergence']=False

        self.bin_fields = bin_fields
        self.extrema = extrema
        self.logs = logs
        self.unit_system = tigress_unit_system

def phase_by_region(filename,params=None,out_dir='',region='midplane',overwrite=False):
    '''
        phase_by_region(filename,out_dir='',region='midplane)
        
        calculate 2D joint PDFs from TIGRESS vtk dump.

        Inputs
        ------
        filename: string
            file name of vtk dump
        
        Parameters
        ----------
        params: class (phase_parameters)
        
        out_dir: string
            directory name for the output hdf5 file
            
        region: string
            setting the region for joint PDFs to be calculated
            'midplane' -- cubic region near the midplane
            'upper' or 'lower' -- region upper or lower than 512 pc
            'full' -- entire domain
    '''
    if params is None: params = phase_parameters()
    tigress_unit_system = params.unit_system

    ds=yt.load(filename,units_override=ya.unit_base,unit_system=tigress_unit_system)
    ya.add_yt_fields(ds,cooling=True,mhd=True,rotation=False)

    le=ds.domain_left_edge
    re=ds.domain_right_edge

    if region == 'midplane':
        le[2]=le[1]
        re[2]=re[1]
    elif region == 'upper':
        le[2]=yt.YTQuantity(512,'pc')
    elif region == 'lower':
        re[2]=yt.YTQuantity(-512,'pc')
    elif region == 'full':
        pass
    else:
        print('{} is unsupported region type.')
        print('region = [midplane, upper, lower, full]'.format(region))
        raise ValueError

    print('Joint PDFs to be calculated for quantities in')
    print('x = {} ~ {}'.format(le[0], re[0]))
    print('y = {} ~ {}'.format(le[1], re[1]))
    print('z = {} ~ {}'.format(le[2], re[2]))
    region=ds.box(le,re)

    phase_part(ds,region,out_dir,params=params,overwrite=overwrite)

def phase_by_slab(filename,params=None,dzslab=128,zmax=0.,verbose=50,overwrite=False):
    '''       
        calculate 2D joint PDFs from TIGRESS vtk dump, slab-by-slab
        
        resulting hdf5 files will be saved to 
            directory-name + 'phase/' + 'zl/zu' + slab-number
            zl, zu for lower and upper halfs
            slab-numbers = zmax_slab/dzslab 
        
        Inputs
        ------
        filename: string
            file name of vtk dump
        
        Parameters
        ----------
        params: class (phase_parameters)

        dzslab: float
            thickness of each slab
            
        zmax: float
            maximum altitude (|z|) to be claculated
            defalut is zmax=Lz/2
        
        verbose: int
            verbose level
    '''
    if params is None: params = phase_parameters()
    tigress_unit_system = params.unit_system

    yt.funcs.mylog.setLevel(verbose)
    ds=yt.load(filename,units_override=ya.unit_base,unit_system=tigress_unit_system)
    ya.add_yt_fields(ds,cooling=True,mhd=True,rotation=False)

    le=ds.domain_left_edge
    re=ds.domain_right_edge
    if zmax == 0. or zmax > re[2].value: zmax=re[2]

    for i,zmin in enumerate(np.arange(0,zmax.value,dzslab)):
        le[2] = yt.YTQuantity(zmin,'pc')
        re[2] = yt.YTQuantity(zmin+dzslab,'pc')
        out_dir=ds.directory.replace('id0','phase/zu{:02d}/'.format(i+1))
        region=ds.box(le,re)

        if verbose <=20:
            print('Joint PDFs to be calculated for quantities in')
            print('x = {} ~ {}'.format(le[0], re[0]))
            print('y = {} ~ {}'.format(le[1], re[1]))
            print('z = {} ~ {}'.format(le[2], re[2]))
            print('Total mass within the region is ',
                  region.sum('cell_mass').in_units('Msun'))
        if verbose <=50:
            print(out_dir)
            
        phase_part(ds,region,out_dir,params=params,overwrite=overwrite)
        
        le[2] = yt.YTQuantity(-zmin-dzslab,'pc')
        re[2] = yt.YTQuantity(-zmin,'pc')
        out_dir=ds.directory.replace('id0','phase/zl{:02d}/'.format(i+1))
        region=ds.box(le,re)

        if verbose <=20:
            print('Joint PDFs to be calculated for quantities in')
            print('x = {} ~ {}'.format(le[0], re[0]))
            print('y = {} ~ {}'.format(le[1], re[1]))
            print('z = {} ~ {}'.format(le[2], re[2]))
            print('Total mass within the region is ',region.sum('cell_mass').in_units('Msun'))
        if verbose <=50: 
            print(out_dir)
            
        phase_part(ds,region,out_dir,params=params,overwrite=overwrite)

def phase_part(ds,region,out_dir,params,overwrite=False):
    '''       
        calculate 2D joint PDFs using yt
                   
        Inputs
        ------
        ds: YTDataSet
            
        region: YTRegion
        
        out_dir: string
            directory name for output hdf files to be saved
            will be created if doesn't exist

        params: class (phase_parameters)

        Parameters
        ----------
        overwrite: bool
            if false, skip it if file exists
    '''
    logs = params.logs
    extrema = params.extrema
    bin_fields = params.bin_fields

    if not os.path.isdir(out_dir): os.mkdir(out_dir)

    for bf in bin_fields:
        xbin,ybin=bf
        outhead1='{}{}-{}-{}'.format(out_dir,ds,xbin,ybin)
        if not os.path.isfile(outhead1+'.h5') or overwrite:
            pdf=yt.create_profile(region,bf,['cell_mass','cell_volume'],
                                  extrema=extrema,
                                  logs=logs,
                                  n_bins=128,weight_field=None)
            pdf.save_as_dataset(outhead1)

def draw_joint_PDFs(pdfs,pdf_ds):
    '''
        draw 1D histograms and 2D joint PDFs from output hdf5 file
         
        Inputs
        ------
        pdfs: dict
         
        pdf_ds: yt.DataSet created by phase_part()         
    '''
    fig, axes = plt.subplots(2,2,figsize=(12,8))
    h2d=[]
    
    for i,field in enumerate(['mass','volume']):
        pdf = pdfs[field]
        xbin = pdf_ds.data[('data','x_bins')]
        ybin = pdf_ds.data[('data','y_bins')]

        # draw histogram along the x-field
        ax = axes[0,0]
        #ax.bar(xbin[:-1],pdf.sum(axis=1)/pdf.sum(),np.diff(xbin),alpha=0.5)
        ax.step(xbin[:-1],pdf.sum(axis=1)/pdf.sum(),lw=2,label=field)
        ax.set_xlabel(pdf_ds.x_field[1])
        ax.set_yscale('log')
        if pdf_ds.x_log: ax.set_xscale('log')
        
        # draw histogram along the y-field
        ax = axes[1,0]
        #ax.bar(ybin[:-1],pdf.sum(axis=0)/pdf.sum(),np.diff(ybin),alpha=0.5)
        ax.step(ybin[:-1],pdf.sum(axis=0)/pdf.sum(),lw=2,label=field)
        ax.set_xlabel(pdf_ds.y_field[1])
        ax.set_yscale('log')
        if pdf_ds.y_log: ax.set_xscale('log')
            
        # draw joind PDF
        ax = axes[i,1]
        im=ax.pcolormesh(xbin,ybin,pdf.T/pdf.sum(),
                         norm=LogNorm(vmin=1.e-7,vmax=1.e-2),
                         cmap=plt.cm.cubehelix_r)

        ax.set_ylabel(pdf_ds.y_field[1]+'\n'+r'$\;[{{\rm {}}}]$'.format(ybin.units))
        ax.set_xlabel(pdf_ds.x_field[1]+'\n'+r'$\;[{{\rm {}}}]$'.format(xbin.units))
        if pdf_ds.x_log: ax.set_xscale('log')
        if pdf_ds.y_log: ax.set_yscale('log')
        h2d.append(im)

    axes[1,0].legend(loc=0)
    time_str='t={:3.1f} Myr'.format(float(pdf_ds.current_time.to('Myr').v))
    axes[0,0].text(0.1,0.9,time_str,transform=axes[0,0].transAxes)

    plt.setp(axes[:,0],'ylabel',r'fraction of gas')
    plt.setp(axes[:,:-1],'ylim',(1.e-5,1.))
    fig.tight_layout()
    for ax,im,field in zip(axes[:,1],h2d,['mass fraction','volume fraction']):
        pos=ax.get_position()
        cax = fig.add_axes([pos.x1,pos.y0,0.02,pos.height])
        fig.colorbar(im,cax=cax,label=field)

    return fig

def plot_joint_PDFs(filename):
    '''       
                                  
        Inputs
        ------
        filename: string
            file name of the hdf5 file created by phase_part() function
        
        Outputs
        -------
        matplotlib figure
    '''
    pdf_ds=yt.load(filename)
    pdfs={'mass':pdf_ds.data['cell_mass'],'volume':pdf_ds.data['cell_volume']}
    fig=draw_joint_PDFs(pdfs,pdf_ds)
    
    return fig

def plot_joint_PDFs_slab(base,pid,itime,fields,pdir=None,slab_index=None,both=True):
    '''       
                                   
        Inputs
        ------
        base: string
        
        pid: string
        
        itime: int
        
        fields: list of string
        
        Parameters
        ----------
        pdir: string
        
        slab_index: list or array
        
        both: bool
        
        Outputs
        -------
        matplotlib figure
    '''
    import glob,os
    
    if pdir is None: pdir=pid
    basedir='{}{}/phase/'.format(base,pdir)
    basename='{}.{:04d}-{}-{}.h5'.format(pid,itime,fields[0],fields[1])
    
    filenames = glob.glob('{}zu??/{}'.format(basedir,basename))
    Nslab=len(filenames)
    if Nslab == 0:
        print(filenames)
        print('cannot find files for {}'.format(basename))
    if slab_index is None: 
        slab_index=range(1,Nslab+1)
    else:
        islab_min = np.min(slab_index)
        islab_max = np.max(slab_index)
        if islab_min < 1: 
            raise ValueError('minimum slab index should be 1')
        if islab_max > Nslab: 
            raise ValueError('maximum slab index should be {}'.format(Nslab))
            
    pdf_ds=yt.load(filenames[0])
    pdf_mass = np.zeros_like(pdf_ds.data['cell_mass'])
    pdf_volume = np.zeros_like(pdf_ds.data['cell_volume'])
    
    for islab in slab_index:
        filename='{}zu{:02d}/{}'.format(basedir,islab,basename)
        pdf_ds=yt.load(filename)
        pdf_mass += pdf_ds.data['cell_mass']
        pdf_volume += pdf_ds.data['cell_volume']
            
        if both:
            filename='{}zl{:02d}/{}'.format(basedir,islab,basename)
            pdf_ds=yt.load(filename)    
            pdf_mass += pdf_ds.data['cell_mass']
            pdf_volume += pdf_ds.data['cell_volume']

    pdfs={'mass':pdf_mass,'volume':pdf_volume}

    fig = draw_joint_PDFs(pdfs,pdf_ds)
    
    return fig
