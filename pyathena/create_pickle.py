from __future__ import print_function

import numpy as np
import pickle
import glob,os

import astropy.constants as c
import pyathena as pa
from pyathena.utils import compare_files
from pyathena.set_units import *
from pyathena.plot_tools.plot_projection import plot_projection
from pyathena.plot_tools.plot_slices import slice2 as plot_slice
from pyathena.plot_tools.set_aux import set_aux

coolftn=pa.coolftn()
unit=pa.set_units(muH=1.4271)

to_Myr=unit['time'].to('Myr').value
to_Pok=(unit['pressure']/c.k_B).cgs.value
to_microG=unit['magnetic_field'].value
to_surf=(unit['density']*unit['length']).to('Msun/pc^2').value

data_axis={'x':2,'y':1,'z':0}
domain_axis={'x':0,'y':1,'z':2}
proj_axis={'z':('x','y'),'y':('x','z'),'x':('y','z')}

def get_scalars(ds):
    scal_fields=[]
    for f in ds.field_list:
        if f.startswith('specific_scalar'):
            scal_fields.append(f)

    scal_fields.sort()
    return scal_fields

def create_surface_density(ds,surf_fname):
    '''
        specific function to create pickle file containing surface density
    '''
    
    time = ds.domain['time']*to_Myr
    dx=ds.domain['dx']

    surf_data={}
    le=ds.domain['left_edge']
    re=ds.domain['right_edge']
    pdata=ds.read_all_data('density')

    proj = np.nanmean(pdata,axis=0)
    proj *= ds.domain['Lx'][2]*to_surf
    bounds = np.array([le[0],re[0],le[1],re[1]])
    surf_data={'time':time,'data':proj,'bounds':bounds}
    pickle.dump(surf_data,open(surf_fname,'wb'),pickle.HIGHEST_PROTOCOL)
    
def create_projection(ds,proj_fname,field='density',conversion=1.0,weight_field=None):
    '''
        generic function to create pickle file containing projections of field along all axes
    
        field: field name to be projected
        
        weighting_field: field name to be multiplied to the projected field
        
        conversion: either be a scalar to be applied to all projections 
            or dict with different factors for each axis projection
    '''
    
    time = ds.domain['time']*to_Myr
    dx=ds.domain['dx']

    surf_data={}
    surf_data['time']=time
    le=ds.domain['left_edge']
    re=ds.domain['right_edge']
    if field is 'temperature':
        pdata=ds.read_all_data('T1')
        pdata=coolftn.get_temp(pdata)
    else:
        pdata=ds.read_all_data(field)

    if weight_field != None:
        if weight_field is 'temperature':
            wdata=ds.read_all_data('T1')
            wdata=coolftn.get_temp(wdata)
        else:
            wdata=ds.read_all_data(weight_field)
        pdata *= wdata

    for i,axis in enumerate(['x','y','z']):
        proj = np.nanmean(pdata,axis=data_axis[axis])
        if weight_field != None:
            wproj = np.nanmean(wdata,axis=data_axis[axis])
            proj /= wproj
        if type(conversion) == dict:
            if axis in conversion:
                proj *= conversion[axis]
        else:
            proj *= conversion
        bounds = np.array([le[domain_axis[proj_axis[axis][0]]],re[domain_axis[proj_axis[axis][0]]],
                           le[domain_axis[proj_axis[axis][1]]],re[domain_axis[proj_axis[axis][1]]]])
        surf_data[axis]={'data':proj,'bounds':bounds}
    pickle.dump(surf_data,open(proj_fname,'wb'),pickle.HIGHEST_PROTOCOL)

def create_slices(ds,slcfname,slc_fields,force_recal=False,factors={}):
    '''
        generic function to create pickle file containing slices of fields
    
        slc_field: list of field names to be sliced
        
        factors: multiplication factors for unit conversion
    '''
 
    time = ds.domain['time']*to_Myr
    dx=ds.domain['dx']
    c=ds.domain['center']
    le=ds.domain['left_edge']
    re=ds.domain['right_edge']
    cidx=pa.cc_idx(ds.domain,ds.domain['center']).astype('int')
 

    import copy
    field_to_slice = copy.copy(slc_fields)
    if os.path.isfile(slcfname) and not force_recal:
        slc_data = pickle.load(open(slcfname,'rb'))
        existing_fields = slc_data['z'].keys()
        for f in existing_fields:
            if f in field_to_slice: 
                print('{} is already there'.format(f))
                field_to_slice.remove(f)
    else:
        slc_data={}
        slc_data['time']=time
       
        for i,axis in enumerate(['x','y','z']):
            bounds = np.array([le[domain_axis[proj_axis[axis][0]]],re[domain_axis[proj_axis[axis][0]]],
                               le[domain_axis[proj_axis[axis][1]]],re[domain_axis[proj_axis[axis][1]]]])
            slc_data[axis]={}
            slc_data[axis+'extent']=bounds/1.e3

    for f in field_to_slice:
        print('slicing {} ...'.format(f))
        if f is 'temperature':
            pdata=ds.read_all_data('T1')
        elif f is 'magnetic_field_strength':
            pdata=ds.read_all_data('magnetic_field')
        elif f is 'ram_pok_z':
            pdata=ds.read_all_data('kinetic_energy3')*2.0
        elif f is 'pok':
            pdata=ds.read_all_data('pressure')
        elif f is 'velocity_z':
            pdata=ds.read_all_data('velocity3')
        elif f is 'velocity_y':
            pdata=ds.read_all_data('velocity2')
        elif f is 'velocity_x':
            pdata=ds.read_all_data('velocity1')
        elif f is 'mag_pok':
            pdata=ds.read_all_data('magnetic_pressure')
        elif f is 'nH':
            pdata=ds.read_all_data('density')
        else:
            pdata=ds.read_all_data(f)

        for i,axis in enumerate(['x','y','z']):
            if f is 'temperature':
                slc=coolftn.get_temp(pdata.take(cidx[i],axis=2-i))
            elif f is 'magnetic_field_strength':
                slc=np.sqrt((pdata.take(cidx[i],axis=2-i)**2).sum(axis=-1))
            else:
                slc=pdata.take(cidx[i],axis=2-i)

            if f in factors:
                slc_data[axis][f] = slc * factors[f]
            else:
                slc_data[axis][f] = slc

    pickle.dump(slc_data,open(slcfname,'wb'),pickle.HIGHEST_PROTOCOL)

def create_all_pickles(do_drawing=False, force_recal=False, force_redraw=False, verbose=True, **kwargs):
    dir = kwargs['base_directory']+kwargs['directory']
    if 'vtk_directory' in kwargs: vtkdir=kwargs['vtk_directory']
    else: vtkdir=dir

    fname=glob.glob(vtkdir+'id0/'+kwargs['id']+'.????.vtk')
    fname.sort()
    if len(fname) == 0:
        fname=glob.glob(vtkdir+kwargs['id']+'.????.vtk')
        print(vtkdir, kwargs['id'], fname)

    if kwargs['range'] != '':
        sp=kwargs['range'].split(',')
        start = eval(sp[0])
        end = eval(sp[1])
        fskip = eval(sp[2])
    else:
        start = 0
        end = len(fname)
        fskip = 1
    fname=fname[start:end:fskip]

    ngrids=len(glob.glob(vtkdir+'id*/'+kwargs['id']+'*'+fname[0][-8:]))

    ds=pa.AthenaDataSet(fname[0])
    mhd='magnetic_field' in ds.field_list
    cooling='pressure' in ds.field_list


    Omega=kwargs['rotation']
    rotation=kwargs['rotation'] != 0.
    if verbose:
        print("MHD:", mhd)
        print("cooling:", cooling)
        print("rotation:", rotation, Omega)

    #slc_fields=['nH','pok','temperature','velocity_z','ram_pok_z']
    slc_fields=['nH','pok','temperature','velocity_x','velocity_y','velocity_z','ram_pok_z']
    fields_to_draw=['star_particles','nH','temperature','pok','velocity_z']
    if mhd:
        slc_fields.append('magnetic_field_strength')
        slc_fields.append('mag_pok')
        fields_to_draw.append('magnetic_field_strength')
    mul_factors={'pok':to_Pok,'magnetic_field_strength':to_microG,'mag_pok':to_Pok,'ram_pok_z':to_Pok}
    scal_fields=get_scalars(ds)
    slc_fields+=scal_fields

    if not os.path.isdir(dir+'slice/'): os.mkdir(dir+'slice/')
    if not os.path.isdir(dir+'surf/'): os.mkdir(dir+'surf/')
    if not os.path.isdir(dir+'proj/'): os.mkdir(dir+'proj/')

    for i,f in enumerate(fname):
        slcfname=dir+'slice/'+kwargs['id']+f[-9:-4]+'.slice.p'
        surfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.surf.p'
        scalfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.scal0.p'
        projfname=dir+'proj/'+kwargs['id']+f[-9:-4]+'.ddproj.p'

        tasks={'slice':(not compare_files(f,slcfname)) or force_recal,
               'surf':(not compare_files(f,surfname)) or force_recal,
               'scal':(not compare_files(f,scalfname)) or force_recal,
               'proj':(not compare_files(f,projfname)) or force_recal,
        }

        do_task=(tasks['slice'] or tasks['surf'] or tasks['scal'] or tasks['proj'])
         
        if verbose: 
            print('file number: {} -- Tasks to be done ['.format(i),end='')
            for k in tasks: print('{}:{} '.format(k,tasks[k]),end='')
            print(']')
        if do_task:
            ds = pa.AthenaDataSet(f)
            if tasks['surf']: create_projection(ds,surfname,conversion={'z':ds.domain['Lx'][2]*to_surf})
            if tasks['slice']: create_slices(ds,slcfname,slc_fields,factors=mul_factors,force_recal=force_recal)
            if tasks['scal']: 
                for nscal,sf in enumerate(scal_fields):
                    nscalfname=scalfname.replace('scal0','scal{}'.format(nscal))
                    create_projection(ds,nscalfname,field=sf,weight_field='density')
            if tasks['proj']: 
                create_projection(ds,projfname,
                                 field='density',weight_field='density')
                create_projection(ds,projfname.replace('ddproj','dTproj'),
                                 field='temperature',weight_field='density')

    aux=set_aux(kwargs['id'])

    if do_drawing:
        for i,f in enumerate(fname):
            slcfname=dir+'slice/'+kwargs['id']+f[-9:-4]+'.slice.p'
            surfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.surf.p'

            starpardir='id0/'
            if os.path.isdir(dir+'starpar/'): starpardir='starpar/'
            starfname=dir+starpardir+kwargs['id']+f[-9:-4]+'.starpar.vtk'

            tasks={'slice':(not compare_files(f,slcfname+'ng')) or force_redraw,
                   'surf':(not compare_files(f,surfname+'ng')) or force_redraw,
            }
            do_task=(tasks['slice'] and tasks['surf'])

            if verbose:
                print('file number: {} -- Tasks to be done ['.format(i),end='')
                for k in tasks: print('{}:{} '.format(k,tasks[k]),end='')
                print(']')
            if tasks['surf']:
                plot_projection(surfname,starfname,runaway=False,aux=aux['surface_density'])
            if tasks['slice']:
                plot_slice(slcfname,starfname,fields_to_draw,aux=aux)
