import pyathena as pa
import glob

def reader(fname,vel=True,mhd=True):
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


