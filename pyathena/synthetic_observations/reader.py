import pyathena as pa
import glob
import numpy as np

def data_reader(fname,vel=True,mhd=True):
    coolftn=pa.coolftn()
    dir, id, step, ext, mpi = pa.parse_filename(fname)
    ds=pa.AthenaDataSet(fname)
    rstfnames=glob.glob('%s/%s.par' % (ds.dir,ds.id))
    par,blocks,fields=pa.parse_par(rstfnames[0])

    domain=ds.domain
    domain['qshear']=eval(par['problem']['qshear'][0])
    domain['Omega']=eval(par['problem']['Omega'][0])
    fields=['density']
    if vel:
        fields.append('velocity1')
        fields.append('velocity2')
        fields.append('velocity3')
    if mhd: 
        fields.append('magnetic_field1')
        fields.append('magnetic_field2')
        fields.append('magnetic_field3')

    data={}
    for f in fields:
        data[f] = ds.read_all_data(f)
    data['temperature']=coolftn.get_temp(ds.read_all_data('T1'))
    fields.append('temperature')
    
    if vel:
        r3d,x3d,y3d,z3d=pa.pos3d(domain)
        vy0=-domain['qshear']*domain['Omega']*x3d
        data['velocity2'] -= vy0

    domain['losdir']=dir+'los/'
    domain['step']=step

    return data,domain

def setup_domain(fname,vel=True,mhd=True,shear=True):
    dir, id, step, ext, mpi = pa.parse_filename(fname)
    ds=pa.AthenaDataSet(fname)
    rstfnames=glob.glob('%s/%s.par' % (ds.dir,ds.id))
    par,blocks,fields=pa.parse_par(rstfnames[0])

    domain=ds.domain
    domain['qshear']=eval(par['problem']['qshear'][0])
    domain['Omega']=eval(par['problem']['Omega'][0])
    fields=['density']
    fields.append('temperature')
    if vel:
        fields.append('velocity1')
        fields.append('velocity2')
        fields.append('velocity3')
    if mhd: 
        fields.append('magnetic_field1')
        fields.append('magnetic_field2')
        fields.append('magnetic_field3')

    domain['fields']=fields
    domain['shear']=shear
    if shear:
        domain['losdir']=dir+'los/'
    else:
        domain['losdir']=dir+'los-periodic/'
    domain['step']=step

    return ds,domain


def read_data(ds,field,domain,vy0_subtract=True):

    if field is 'temperature':
        coolftn=pa.coolftn()
        data=coolftn.get_temp(ds.read_all_data('T1'))
    elif field is 'velocity2':
        data = ds.read_all_data(field)
        if vy0_subtract:
            r3d,x3d,y3d,z3d=pa.pos3d(domain)
            vy0=-domain['qshear']*domain['Omega']*x3d
            data -= vy0
    else:
        data = ds.read_all_data(field)
    
    return data

def select_phase(temperature,density,phase='warm'):

    T1=5050.
    T2=2.e4
    dmax=50
    if phase is 'whole': 
        return density
    elif phase is 'warm':
        idx = (temperature < T1) | (temperature > T2)
    elif phase is 'cold': 
        idx = (temperature >= T1) 
    elif phase is '2p': 
        idx = (temperature > T2) 
    elif phase is 'lowd': 
        idx = (temperature > T2) | (density > dmax)
    else:
        print("{} is not supported".format(phase))
        return -1
    
    dnew = np.copy(density)
    dnew[idx] = 0.0

    return dnew
