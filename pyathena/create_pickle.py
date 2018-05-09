from __future__ import print_function

import numpy as np
import cPickle as pickle
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

    proj = pdata.mean(axis=0)
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
    pdata=ds.read_all_data(field)
    if weight_field != None:
        wdata=ds.read_all_data(weight_field)
        pdata *= wdata
    for i,axis in enumerate(['x','y','z']):
        proj = pdata.mean(axis=data_axis[axis])
        if weight_field != None:
            wproj = wdata.mean(axis=data_axis[axis])
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

def create_slices(ds,slcfname,slc_fields,factors={}):
    '''
        generic function to create pickle file containing slices of fields
    
        slc_field: list of fields name to be sliced
        
        factors: multiplication factors for unit conversion
    '''
 
    time = ds.domain['time']*to_Myr
    dx=ds.domain['dx']

    slc_data={}
    slc_data['time']=time
    c=ds.domain['center']
    le=ds.domain['left_edge']
    re=ds.domain['right_edge']
    cidx=pa.cc_idx(ds.domain,ds.domain['center']).astype('int')
    
    for i,axis in enumerate(['x','y','z']):
        bounds = np.array([le[domain_axis[proj_axis[axis][0]]],re[domain_axis[proj_axis[axis][0]]],
                           le[domain_axis[proj_axis[axis][1]]],re[domain_axis[proj_axis[axis][1]]]])
        slc_data[axis]={}
        slc_data[axis+'extent']=bounds/1.e3

    for f in slc_fields:
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

def create_all_pickles(force_recal=False, force_redraw=False, verbose=True, **kwargs):
    dir = kwargs['base_directory']+kwargs['directory']
    fname=glob.glob(dir+'id0/'+kwargs['id']+'.????.vtk')
    fname.sort()

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

    ngrids=len(glob.glob(dir+'id*/'+kwargs['id']+'*'+fname[0][-8:]))

    ds=pa.AthenaDataSet(fname[0])
    mhd='magnetic_field' in ds.field_list
    cooling='pressure' in ds.field_list


    Omega=kwargs['rotation']
    rotation=kwargs['rotation'] != 0.
    if verbose:
        print("MHD:", mhd)
        print("cooling:", cooling)
        print("rotation:", rotation, Omega)

    slc_fields=['nH','pok','temperature','velocity_z','ram_pok_z']
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

    for i,f in enumerate(fname):
        slcfname=dir+'slice/'+kwargs['id']+f[-9:-4]+'.slice.p'
        surfname=dir+'surf/'+kwargs['id']+f[-9:-4]+'.surf.p'

        tasks={'slice':(not compare_files(f,slcfname)) or force_recal,
               'surf':(not compare_files(f,surfname)) or force_recal,
        }

        do_task=(tasks['slice'] or tasks['surf'])
         
        if verbose: 
            print('file number: {} -- Tasks to be done ['.format(i),end='')
            for k in tasks: print('{}:{} '.format(k,tasks[k]),end='')
            print(']')
        if do_task:
            ds = pa.AthenaDataSet(f)
            if tasks['surf']: create_projection(ds,surfname,conversion={'z':ds.domain['Lx'][2]*to_surf})
            if tasks['slice']: create_slices(ds,slcfname,slc_fields,factors=mul_factors)

    aux=set_aux(kwargs['id'])

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
            plot_projection(surfname,starfname,runaway=True,aux=aux['surface_density'])
        if tasks['slice']:
            plot_slice(slcfname,starfname,fields_to_draw,aux=aux)
