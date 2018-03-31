import numpy as np
import struct
import string
import glob
import os
import re
import pandas as pd
from itertools import *
import astropy.constants as c
import astropy.units as u
import pickle as p

from ath_hst import test_pickle
from utils import *
from set_units import set_units

def parse_filename(filename):
    """
    #   PARSE_FILENAME    Break up a full-path filename into its component
    #   parts to check the extension, make it more readable, and extract the step
    #   number.  
    #
    #   PATH,BASENAME,STEP,EXT = PARSE_FILENAME(FILENAME)
    #
    #   E.g. If FILENAME='/home/Blast.0000.bin', then PATH='/home',
    #   BASENAME='Blast', STEP='0000', and EXT='bin'.
    #
    """


    path=os.path.dirname(filename)
    if path[-3:] == 'id0': 
        path=path[:-3]
        mpi_mode=True
    else:
        path=path+os.path.sep
        mpi_mode=False

    base=os.path.basename(filename)
    base_split=base.split('.')
    if len(base_split) == 3:
        id=base_split[0]
        step=base_split[1]
        ext=base_split[2]
    else:
        id=string.join(base_split[:-2],'.')
        step=base_split[-2]
        ext=base_split[-1]

    return path,id,step,ext,mpi_mode

def parse_line(line, grid):
    sp = line.strip().split()

    if b"vtk" in sp:
        grid['vtk_version'] = sp[-1]
    elif b"time=" in sp:
        time_index = sp.index(b"time=")
        grid['time'] = float(sp[time_index+1].rstrip(b','))
        if b'level' in sp: grid['level'] = int(sp[time_index+3].rstrip(b','))
        if b'domain' in sp: grid['domain'] = int(sp[time_index+5].rstrip(b','))  
        if sp[0] == b"PRIMITIVE": 
            grid['prim_var_type']=True
    elif b"DIMENSIONS" in sp:
        grid['Nx'] = np.array(sp[-3:]).astype('int')
    elif b"ORIGIN" in sp:
        grid['left_edge'] = np.array(sp[-3:]).astype('float64')
    elif b"SPACING" in sp:
        grid['dx'] = np.array(sp[-3:]).astype('float64')
    elif b"CELL_DATA" in sp:
        grid['ncells'] = int(sp[-1])
    elif b"SCALARS" in sp:
        grid['read_field'] = sp[1]
        grid['read_type'] = 'scalar'
    elif b"VECTORS" in sp:
        grid['read_field'] = sp[1]
        grid['read_type'] = 'vector'
    elif b"NSTARS" in sp:
        grid['nstar'] = eval(sp[1])
    elif b"POINTS" in sp:
        grid['nstar'] = eval(sp[1])
        grid['ncells'] = eval(sp[1])


class AthenaDomain(object):
    def __init__(self,filename,ds=None,setgrid=True,serial=False):
        self.flist = glob.glob(filename)
        if len(self.flist) == 0: 
            print(('no such file: %s' % filename))
        dir, id, step, ext, mpi = parse_filename(filename)
        self.dir = dir
        self.id = id
        self.step = step
        self.ext = ext
        self.starfile = os.path.join(dir+'id0/','%s.%s.%s.%s' % (id,step,'starpar',ext))
        if serial: mpi = False
        self.mpi = mpi
        self.ngrids = 1
        if mpi:
            self.ngrids += len(glob.glob(os.path.join(dir,'id*/%s-id*.%s.%s' % (id, step, ext))))
            for n in range(1,self.ngrids):
                self.flist.append(os.path.join(dir,'id%d/%s-id%d.%s.%s' % (n,id,n,step, ext)))
        if setgrid:
            if ds==None: 
                self.grids=self._setup_grid()
            else: 
                if ds.grids[0]['filename'] != self.flist[0]:
                    for g,f in zip(ds.grids,self.flist): g['filename']=f
                self.grids=ds.grids
            self.domain=self._setup_domain(self.grids)
            if ds==None: 
                self.domain['field_map']=None
            else:
                self.domain['field_map']=ds.domain['field_map']
            self._setup_mpi_grid()
            self._setup()
    
    def _setup(self):
        self.domain['data']={}

    def _setup_domain(self,grids):
        domain = {}
        ngrids=len(grids)
        left_edges = np.empty((ngrids,3), dtype='float32')
        dxs = np.empty((ngrids,3), dtype='float32')
        Nxs = np.ones_like(dxs)
        for nproc,g in enumerate(grids):
            left_edges[nproc,:] = g['left_edge']
            Nxs[nproc,:] = g['Nx']
            dxs[nproc,:] = g['dx']

        right_edges = left_edges + Nxs*dxs

        left_edge = left_edges.min(0)
        right_edge = right_edges.max(0)

        domain['left_edge'] = left_edge
        domain['right_edge'] = right_edge
        domain['dx'] = dxs[0,:]
        domain['Lx'] = right_edge - left_edge
        domain['center'] = 0.5*(right_edge + left_edge)
        domain['Nx'] = np.round(domain['Lx']/domain['dx']).astype('int')
        domain['ndim'] = 3 # should be revised
        file = open(self.flist[0],'rb')
        tmpgrid = {}
        tmpgrid['time']=None
        while tmpgrid['time'] is None:
            line = file.readline()
            parse_line(line,tmpgrid)
        file.close()
        domain['time'] = tmpgrid['time']

        return domain

    def _setup_mpi_grid(self):
        gnx = self.grids[0]['Nx']
        self.NGrids = (self.domain['Nx']/self.grids[0]['Nx']).astype(np.int)
        self.gid = np.arange(self.ngrids)
        i = 0
        for n in range(self.NGrids[2]):
            for m in range(self.NGrids[1]):
                for l in range(self.NGrids[0]):
                    self.grids[i]['is']=np.array([l*gnx[0],m*gnx[1],n*gnx[2]])
                    i += 1 

    def _setup_grid(self):
        grids=[]
        for nproc in range(self.ngrids):
            file = open(self.flist[nproc],'rb')
            grid = {}
            grid['filename']=self.flist[nproc]
            grid['read_field'] = None
            grid['read_type'] = None
            while grid['read_field'] is None:
                grid['data_offset']=file.tell()
                line = file.readline()
                parse_line(line, grid)
            file.close()
            grid['Nx'] -= 1
            grid['Nx'][grid['Nx'] == 0] = 1
            grid['dx'][grid['Nx'] == 1] = 1.
            grid['right_edge'] = grid['left_edge'] + grid['Nx']*grid['dx']
            #grid['field_map']=None

            grids.append(grid)

        return grids

class AthenaZprof(object):
    def __init__(self,filename,stitch=True,clean=False,forced_clean=False):
        path=os.path.dirname(filename)
        if path[-3:] == 'id0': 
            path=path[:-3]
            mpi=True
        else:
            path=path+os.path.sep
            mpi=False
        if filename.endswith('.p'):
            base=os.path.basename(filename[:-2])
        else:
            base=os.path.basename(filename)
        base_sp=base.split('.')
        ext = base_sp[-1]
        phase = base_sp[-2]
        step = base_sp[-3]
        id = string.join(base_sp[:-3],'.')
        self.path = path
        self.id = id
        self.step = step
        self.ext = ext
        self.mpi = mpi
        self.filename = filename
        #print id,step,phase,ext
        self._set_plist()
        if not filename.endswith('.p'): 
            self._set_time()
        if stitch:
            for p in self.plist: self._stitch(p)
        if clean:
            for p in self.plist: self._clean(p,forced=forced_clean)

    def _set_time(self):
        fp=open(self.filename,'r')
        line=fp.readline()
        fp.close()
        self.time=eval(re.split('=|\s',line)[-2])


    def _set_plist(self):
        zpfile = glob.glob(os.path.join(self.path,'%s.%s.*.%s' % (self.id, self.step, self.ext)))
        zpfile.sort()
        self.plist=[]
        for f in zpfile: 
            phase=f.split('.')[-2]
            self.plist.append(phase)

    def _stitch(self,phase):
        fname=self.zpfile(phase)
        flist = [fname]
        #nproc = len(glob.glob(os.path.join(self.path,'id*')))
        if self.mpi: flist += glob.glob(os.path.join(self.path,'id*/%s-id*.%s.%s.%s' % (self.id, self.step, phase, self.ext)))
        nf=len(flist)    
        if not test_pickle(fname): 
            if nf > 1:
#                for i in range(1,nproc):
#                    flist.append(os.path.join(self.path,'id%d/%s-id%d.%s.%s.%s' % (i, self.id, i, self.step, phase, self.ext)))
                flist.sort()
                print(("Read-Merge-Write for ",phase))
                data=pd.DataFrame()
                for f in flist:
                    df=pd.read_csv(f,comment='#',index_col=0)
                    data=data.add(df,fill_value=0)
                data.to_pickle(fname+'.p')
            else:
                print(("Read-Write for ",phase))
                df=pd.read_csv(fname,comment='#',index_col=0)
                df.to_pickle(fname+'.p')

    def _clean(self,phase,forced=False):
        fname=self.zpfile(phase)
        if test_pickle(fname) or forced: 
            if self.mpi: flist = glob.glob(os.path.join(self.path,'id*/%s-id*.%s.%s.%s' % (self.id, self.step, phase, self.ext)))
            if len(flist) != 0: 
                print(("Clean for ",phase))
                for f in flist: os.remove(f)

    def zpfile(self,phase):
        f_sp=self.filename.split('.')
        f_sp[-2]=phase
        return string.join(f_sp,'.')

    def read(self,phase='all'):
        if phase=='all': 
            zprof={}
            for p in self.plist:
                zprof[p]=pd.read_pickle(self.zpfile(p)+'.p')
        elif phase in self.plist: 
            zprof=pd.read_pickle(self.zpfile(phase)+'.p')
        else: 
            print((phase," is not defined"))

        return zprof
            
class AthenaDataSet(AthenaDomain):
    def _setup(self):
        for i,g in enumerate(self.grids):
#            if g['field_map']==None: 
#                print "Setting %d-th grid" % (i)
#                self._set_field_map(g)
#            else: 
            g['data']={}
        self.units=set_units()
        if self.domain['field_map']==None:
            self.domain['field_map'] = self._set_field_map(self.grids[0])
        fm = self.domain['field_map']
        if 'cell_centered_B' in list(fm.keys()):
            fm['magnetic_field']=fm['cell_centered_B']
            fm.pop('cell_centered_B')
        elif 'face_centered_B' in list(fm.keys()):
            fm['magnetic_field']=fm['face_centered_B']
            fm.pop('face_centered_B')
        nscal=0
        if 'specific_scalar[0]' in list(fm.keys()):
            keys=list(fm.keys())
            for k in keys:
                if k.startswith('specific_scalar'):
                    newkey=re.sub("\s|\W","",k)
                    fm[newkey] = fm.pop(k)
                    nscal += 1
        self.field_list=list(fm.keys())

        
        derived_field_list=[]
        derived_field_list_hd=[]
        derived_field_list_mhd=[]
        if 'magnetic_field' in self.field_list:
            derived_field_list_mhd.append('magnetic_field1')
            derived_field_list_mhd.append('magnetic_field2')
            derived_field_list_mhd.append('magnetic_field3')
            derived_field_list_mhd.append('magnetic_energy1')
            derived_field_list_mhd.append('magnetic_energy2')
            derived_field_list_mhd.append('magnetic_energy3')
            derived_field_list_mhd.append('magnetic_pressure')
            derived_field_list_mhd.append('plasma_beta')
            derived_field_list_mhd.append('alfven_velocity1')
            derived_field_list_mhd.append('alfven_velocity2')
            derived_field_list_mhd.append('alfven_velocity3')
            derived_field_list_mhd.append('magnetic_stress')
            derived_field_list_mhd.append('magnetic_stress1')
            derived_field_list_mhd.append('magnetic_stress2')
            derived_field_list_mhd.append('magnetic_stress3')
        if 'velocity' in self.field_list:
            derived_field_list_hd.append('velocity1')
            derived_field_list_hd.append('velocity2')
            derived_field_list_hd.append('velocity3')
            derived_field_list_hd.append('velocity_magnitude')
            derived_field_list_hd.append('kinetic_energy1')
            derived_field_list_hd.append('kinetic_energy2')
            derived_field_list_hd.append('kinetic_energy3')
            derived_field_list_hd.append('momentum1')
            derived_field_list_hd.append('momentum2')
            derived_field_list_hd.append('momentum3')
            derived_field_list_hd.append('reynold_stress')
            derived_field_list_hd.append('reynold_stress1')
            derived_field_list_hd.append('reynold_stress2')
            derived_field_list_hd.append('reynold_stress3')
        if 'pressure' in self.field_list:
            derived_field_list_hd.append('sound_speed')
            derived_field_list_hd.append('temperature')
            derived_field_list_hd.append('T1')
        if 'gravitational_potential' in self.field_list:
            derived_field_list_hd.append('potential_energy')
            derived_field_list_hd.append('gravity_stress')
            derived_field_list_hd.append('gravity_stress1')
            derived_field_list_hd.append('gravity_stress2')
            derived_field_list_hd.append('gravity_stress3')
        derived_field_list_hd.append('number_density')
        if nscal > 0:
            for n in range(nscal):
                derived_field_list_hd.append('scalar%d' % n)
        self.domain['nscal']=nscal
        self.derived_field_list=derived_field_list_hd+derived_field_list_mhd
        self.derived_field_list_hd=derived_field_list_hd
        self.derived_field_list_mhd=derived_field_list_mhd
        for f in self.derived_field_list+self.field_list:
            if f not in self.units: self.units[f]=1.0*u.dimensionless_unscaled

    def _set_field_map(self,grid):
        return set_field_map(grid)

    def _read_field(self,file_pointer,field_map):
        return read_field(file_pointer,field_map)

    def _read_grid_data(self,grid,field):
        gd=grid['data']
        if field in gd:
            return

        file=open(grid['filename'],'rb')
        #fm=grid['field_map']
        fm=self.domain['field_map']
        nx1=grid['Nx'][0]
        nx2=grid['Nx'][1]
        nx3=grid['Nx'][2]

        if field == 'face_centered_B1': nx1=nx1+1
        if field == 'face_centered_B2': nx2=nx2+1
        if field == 'face_centered_B3': nx3=nx3+1

        nvar=fm[field]['nvar']
        var = self._read_field(file,fm[field])
        if nvar == 1: 
            var.shape = (nx3, nx2, nx1)
        else: 
            var.shape = (nx3, nx2, nx1, nvar)
        file.close()
        grid['data'][field]=var
        if nvar == 3: self._set_vector_field(grid,self.units[field],field)


    def _get_grid_data(self,grid,field):
        gd=grid['data']
        if field in gd:
            return gd[field]
        elif field in self.field_list:
            self._read_grid_data(grid,field)
            return gd[field]
        elif field in self.derived_field_list:
            data=self._get_derived_field(grid,field)
            return data
        else:
            print((field,' is not supported'))

    def _set_vector_field(self,grid,unit,vfield):
        gd=grid['data']
        gd[vfield+'1'] = gd[vfield][:,:,:,0]
        gd[vfield+'2'] = gd[vfield][:,:,:,1]
        gd[vfield+'3'] = gd[vfield][:,:,:,2]
        self.units[vfield+'1'] = unit
        self.units[vfield+'2'] = unit
        self.units[vfield+'3'] = unit

    def _get_derived_field(self,grid,field):
        import astropy.constants as c
        u=self.units
        gd=grid['data']
        if field in gd:
            return gd[field]
        elif field.startswith('velocity'):
            self._read_grid_data(grid,'velocity')
            if field is 'velocity_magnitude': 
                v1=gd['velocity1']
                v2=gd['velocity2']
                v3=gd['velocity3']
                vmag=np.sqrt(v1**2+v2**2+v3**2)
                return vmag
            else: return gd[field]
        elif field.startswith('magnetic_field'):
            self._read_grid_data(grid,'magnetic_field')
            return gd[field]
        elif field.startswith('number_density'):
            self._read_grid_data(grid,'density')
            return gd['density']
        elif field.startswith('kinetic_energy'):
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,'velocity')
            den=gd['density']
            v1=gd['velocity1']
            v2=gd['velocity2']
            v3=gd['velocity3']
            u['kinetic_energy1']=u['density']*u['velocity']**2
            u['kinetic_energy2']=u['density']*u['velocity']**2
            u['kinetic_energy3']=u['density']*u['velocity']**2
            if field is 'kinetic_energy1': return 0.5*den*v1**2
            if field is 'kinetic_energy2': return 0.5*den*v2**2
            if field is 'kinetic_energy3': return 0.5*den*v3**2
        elif field.startswith('momentum'):
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,'velocity')
            den=gd['density']
            v1=gd['velocity1']
            v2=gd['velocity2']
            v3=gd['velocity3']
            u['momentum1']=u['density']*u['velocity']
            u['momentum2']=u['density']*u['velocity']
            u['momentum3']=u['density']*u['velocity']
            if field is 'momentum1': return den*v1
            if field is 'momentum2': return den*v2
            if field is 'momentum3': return den*v3
        elif field.startswith('magnetic_energy'):
            self._read_grid_data(grid,'magnetic_field')
            B1=gd['magnetic_field1']
            B2=gd['magnetic_field2']
            B3=gd['magnetic_field3']
            u['magnetic_energy1']=u['magnetic_field']**2
            u['magnetic_energy2']=u['magnetic_field']**2
            u['magnetic_energy3']=u['magnetic_field']**2
            if field is 'magnetic_energy1': return B1**2/(8*np.pi)
            if field is 'magnetic_energy2': return B2**2/(8*np.pi)
            if field is 'magnetic_energy3': return B3**2/(8*np.pi)
        elif field.startswith('magnetic_pressure'):
            self._read_grid_data(grid,'magnetic_field')
            B1=gd['magnetic_field1']
            B2=gd['magnetic_field2']
            B3=gd['magnetic_field3']
            u['magnetic_pressure']=u['magnetic_field']**2
            if field is 'magnetic_pressure': return (B1**2+B2**2+B3**2)/(8*np.pi)
        elif field.startswith('plasma_beta'):
            vfield='magnetic_field'
            unit=u[vfield]
            self._read_grid_data(grid,'pressure')
            self._read_grid_data(grid,vfield)
            B1=gd[vfield+'1']
            B2=gd[vfield+'2']
            B3=gd[vfield+'3']
            press=gd['pressure']
            u['plasma_beta']=u['pressure']/unit**2
            if field is 'plasma_beta': return press*(8.0*np.pi)/(B1**2+B2**2+B3**2)
        elif field.startswith('alfven_velocity'):
            vfield='magnetic_field'
            unit=u[vfield]
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,vfield)
            den=gd['density']
            B1=gd[vfield+'1']
            B2=gd[vfield+'2']
            B3=gd[vfield+'3']
            u['alfven_velocity1']=np.sqrt(unit**2/u['density'])
            u['alfven_velocity2']=np.sqrt(unit**2/u['density'])
            u['alfven_velocity3']=np.sqrt(unit**2/u['density'])
            if field is 'alfven_velocity1': return B1/np.sqrt(4*np.pi*den)
            if field is 'alfven_velocity2': return B2/np.sqrt(4*np.pi*den)
            if field is 'alfven_velocity3': return B3/np.sqrt(4*np.pi*den)
        elif field.startswith('sound_speed'):
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,'pressure')
            den=gd['density']
            press=gd['pressure']
            u['sound_speed']=np.sqrt(u['pressure']/u['density'])
            return np.sqrt(press/den)
        elif field.startswith('temperature'):
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,'pressure')
            den=gd['density']*u['density']
            press=gd['pressure']*u['pressure']
            T1=(press/den*c.m_p/c.k_B).cgs
            return T1*1.4/1.1
        elif field.startswith('T1'):
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,'pressure')
            den=gd['density']*u['density']
            press=gd['pressure']*u['pressure']
            u['T1']=u['temperature']
            return (press/den*c.m_p/c.k_B).cgs
        elif field.startswith('potential'):
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,'gravitational_potential')
            den=gd['density']
            pot=gd['gravitational_potential']
            u['potential_energy']=u['density']*u['gravitational_potential']
            return -den*pot
        elif field.startswith('magnetic_stress'):
            vfield='magnetic_field'
            unit=u[vfield]
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,vfield)
            B1=gd[vfield+'1']
            B2=gd[vfield+'2']
            B3=gd[vfield+'3']
            u['magnetic_stress']=unit**2
            u['magnetic_stress1']=unit**2
            u['magnetic_stress2']=unit**2
            u['magnetic_stress3']=unit**2
            if field is 'magnetic_stress1': return B2*B3/(4*np.pi) 
            if field is 'magnetic_stress2': return B1*B3/(4*np.pi) 
            if field is 'magnetic_stress3': return B1*B2/(4*np.pi) 
            return B1*B2/(4*np.pi)
        elif field.startswith('reynold_stress'):
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,'velocity')
            den=gd['density']
            v1=gd['velocity1']
            v2=gd['velocity2']
            v3=gd['velocity3']
            u['reynold_stress']=u['density']*u['velocity']**2
            u['reynold_stress1']=u['density']*u['velocity']**2
            u['reynold_stress2']=u['density']*u['velocity']**2
            u['reynold_stress3']=u['density']*u['velocity']**2
            if field is 'reynold_stress1': return den*v2*v3
            if field is 'reynold_stress2': return den*v1*v3
            if field is 'reynold_stress3': return den*v1*v2
            return den*v1*v2
        elif field.startswith('gravity_stress'):
            self._read_grid_data(grid,'gravitational_potential')
            phi=gd['gravitational_potential']
            dx=grid['dx']
            u['gravity_stress']=u['gravitational_potential']**2/u['length']**2/c.G.cgs
            u['gravity_stress1']=u['gravitational_potential']**2/u['length']**2/c.G.cgs
            u['gravity_stress2']=u['gravitational_potential']**2/u['length']**2/c.G.cgs
            u['gravity_stress3']=u['gravitational_potential']**2/u['length']**2/c.G.cgs
            g1,g2,g3=gradient(phi,dx)
            if field is 'gravity_stress1': return g2*g3/4/np.pi
            if field is 'gravity_stress2': return g1*g3/4/np.pi
            if field is 'gravity_stress3': return g1*g2/4/np.pi
            return  g1*g2/4/np.pi
        elif field.startswith('scalar'):
            scal = field[6:]
            self._read_grid_data(grid,'density')
            self._read_grid_data(grid,'specific_scalar'+scal)
            den=gd['density']
            sscal=gd['specific_scalar'+scal]
            return sscal*den
        
    def _set_data_array(self,field,dnx):
        fm=self.domain['field_map']
        
        if field in self.field_list:
            if fm[field]['nvar']==3:
                data=np.empty((dnx[2],dnx[1],dnx[0],3),dtype=fm[field]['dtype'])
            else:
                data=np.empty((dnx[2],dnx[1],dnx[0]),dtype=fm[field]['dtype'])
            if field is 'face_centered_B1':
                data=np.empty((dnx[2],dnx[1],dnx[0]+1),dtype=fm[field]['dtype'])
            if field is 'face_centered_B2':
                data=np.empty((dnx[2],dnx[1]+1,dnx[0]),dtype=fm[field]['dtype'])
            if field is 'face_centered_B3':
                data=np.empty((dnx[2]+1,dnx[1],dnx[0]),dtype=fm[field]['dtype'])
        elif field in self.derived_field_list:
            data=np.empty((dnx[2],dnx[1],dnx[0]),dtype=fm['density']['dtype'])
        return data

    def _get_slab_grid(self,slab=1,verbose=False):
        if slab > self.NGrids[2]: 
            print(("%d is lareger than %d" % (slab,self,NGrids[2])))
        NxNy=self.NGrids[0]*self.NGrids[1]
        gidx, = np.where(slab == self.gid/NxNy+1)
        grids = []
        for i in gidx:
            grids.append(self.grids[i])
        if verbose: print(("XY slab from z=%g to z=%g" % (grids[0]['left_edge'][2],grids[0]['right_edge'][2])))
        return grids

    def read_all_data(self,field,slab=False,verbose=False):
        #fm=self.grids[0]['field_map']
        fm=self.domain['field_map']
        dnx=np.copy(self.domain['Nx'])
        if slab: 
            dnx[2]=self.grids[0]['Nx'][2]
            grids = self._get_slab_grid(slab=slab,verbose=verbose)
        else:
            grids = self.grids
        data = self._set_data_array(field,dnx)
        for g in grids:
            gis=np.copy(g['is'])
            if slab: gis[2]=0
            gnx=np.copy(g['Nx'])
            gie=gis+gnx
            gd=self._get_grid_data(g,field)
            if field in self.field_list and fm[field]['nvar']==3:
                data[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0],:]=gd
            else:
                if gie[0] == dnx[0] and field is 'face_centered_B1':
                    data[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]+1]=gd
                elif gie[1] == dnx[1] and field is 'face_centered_B2':
                    data[gis[2]:gie[2],gis[1]:gie[1]+1,gis[0]:gie[0]]=gd
                elif gie[2] == dnx[2] and field is 'face_centered_B3':
                    data[gis[2]:gie[2]+1,gis[1]:gie[1],gis[0]:gie[0]]=gd
                else:
                    gd=gd[0:gnx[2],0:gnx[1],0:gnx[0]]
                    data[gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]=gd
            
        return data


class AthenaRegion(object):
    def __init__(self,ds,*args,**kwargs):
        self._region_selector(ds,*args,**kwargs)

    def _region_selector(self,ds,*args,**kwargs):
        return

    def _plane_grid_selector(self,ds,region):

        grid_list=[]
        for g in ds.grids:
            axis=region['axis']
            c=region['center'][axis]
            le=region['left_edge']
            re=region['right_edge']
            aidx=region['axis_idx']

            zmin=g['left_edge'][axis]
            zmax=g['right_edge'][axis]
            gle=g['left_edge'][aidx]
            gre=g['right_edge'][aidx]
            
            if (zmin-c)*(zmax-c)<=0:
                dle = gle-le
                dre = gre-re
                if dle[0]*dre[0]<=0 and dle[1]*dre[1]<=0:
                    #print g['left_edge'],g['right_edge'],c
                    grid_list.append(g)

        #print len(grid_list)
        return grid_list

    def _cubic_region_selector(self,ds,le,re):

        grid_list=[]
        le=np.array(le)
        re=np.array(re)
        for g in ds.grids:
            gle=g['left_edge']
            gre=g['right_edge']
            
            # get all grids with gre > le
            if (gre > le).all() and (gle <re).all():
                print((g['left_edge'],g['right_edge']))
                grid_list.append(g)
        self.grid_list=grid_list
        self.le=le
        self.re=re
        self.ngrid=len(grid_list)
        return grid_list

    
class AthenaSlice(AthenaRegion):
    def _region_selector(self,ds,*args,**kwargs):
        self.slice(ds,*args,**kwargs)

    def slice(self,ds,*args,**kwargs):
        if 'axis' in kwargs: axis=kwargs['axis']
        else: print('need to specify symmetric axis')

        if 'center' in kwargs: c=kwargs['center']
        else: print('need to specify center of the plane')

        if 'field' in kwargs: field=kwargs['field']
        else: print('need to specify field')

        if axis == 'x': axis = 2
        elif axis == 'y': axis = 1
        elif axis == 'z': axis = 0

        if c=='center': c=ds.domain['center']

        aidx=[0,1,2]
        aidx.remove(2-axis)
        self.axis_labels=['x','y','z']
        self.axis_labels.remove(self.axis_labels[2-axis])
        le=ds.domain['left_edge'][aidx]
        re=ds.domain['right_edge'][aidx]
        nx=ds.domain['Nx'][aidx]
        self.time=ds.domain['time']
        self.field=field
        self.wfield=None

        region={'axis':2-axis,'center':c,'left_edge':le,'right_edge':re,'axis_idx':aidx}
        self.data=np.empty((nx[1],nx[0]),dtype='float32')
        self.bound=(le[0],re[0],le[1],re[1])

        grid_list=self._plane_grid_selector(ds,region)

        #ng=0
        for g in grid_list:
            gd=g['data']
            cidx=np.array(cc_idx(ds.domain,c)-g['is'],dtype='int')
            #print "Reading:",g['filename']
            #print cidx

            gis=g['is'][aidx]
            gnx=g['Nx'][aidx]
            gie=gis+gnx
            if field not in gd: 
                #ng = ng+1
                #print "Reading:",g['filename'],ng
                gdata=ds._get_grid_data(g,field)
            else:
                gdata=gd[field]
            if axis==2:
                data=gdata[:,:,cidx[0]]
            if axis==1:
                data=gdata[:,cidx[1],:]
            if axis==0:
                data=gdata[cidx[2],:,:]
            self.data[gis[1]:gie[1],gis[0]:gie[0]]=data

        self.units=ds.units[field]

class AthenaSurf(object):
    def __init__(self,ds,*args,**kwargs):
        if 'axis' in kwargs: axis=kwargs['axis']
        else: axis=1
        if 'field' in kwargs: field=kwargs['field']
        else: field='density'
        if 'weighted_field' in kwargs: wfield=kwargs['weighted_field']
        else: wfield=None

        self.axis=axis
        self.field=field
        self.wfield=wfield
        self.time=ds.domain['time']

        self.data = self._get_mean(ds,axis,field,wfield)
        self._setup(ds)
        self.units=ds.units[field]

    def _setup(self,ds):
        aidx=[0,1,2]
        aidx.remove(2-self.axis)
        self.axis_labels=['x','y','z']
        self.axis_labels.remove(self.axis_labels[2-self.axis])
        le=ds.domain['left_edge'][aidx]
        re=ds.domain['right_edge'][aidx]
        nx=ds.domain['Nx'][aidx]

        self.bound=(le[0],re[0],le[1],re[1])

    def _get_mean(self,ds,axis,field,wfield):
        data=ds.read_all_data(field)
        if wfield != None:
            wdata=ds.read_all_data(wfield)
            data *= wdata
            wdata *= ds.domain['dx'][2-axis]
            wdata = wdata.sum(axis=axis)
        else:
            wdata = ds.domain['Lx'][2-axis]
        data *= ds.domain['dx'][2-axis]
        data = data.sum(axis=axis)

        data = data/wdata

        return data

def set_field_map(grid):
    file=open(grid['filename'],'rb')
    file.seek(0,2)
    eof = file.tell()
    offset = grid['data_offset']
    file.seek(offset)

    field_map={}

    if 'Nx' in grid: Nx=grid['Nx']

    while offset < eof:

        line=file.readline()
        sp = line.strip().split()
        if len(sp) == 0:
            line=file.readline()
            sp = line.strip().split()
        #print line,sp,len(sp)
        field=sp[1].decode('utf-8')    
        field_map[field] = {}
        field_map[field]['read_table']=False

        if b"SCALARS" in line:
            tmp=file.readline()
            field_map[field]['read_table']=True
            field_map[field]['nvar'] = 1
        elif b"VECTORS" in line:
            field_map[field]['nvar'] = 3
        else:
            print(('Error: '+sp[0] + ' is unknown type'))
            raise TypeError

        field_map[field]['offset']=offset
        field_map[field]['ndata']=field_map[field]['nvar']*grid['ncells']
        if field == 'face_centered_B1':
            field_map[field]['ndata']=(Nx[0]+1)*Nx[1]*Nx[2]
        elif field == 'face_centered_B2':
            field_map[field]['ndata']=Nx[0]*(Nx[1]+1)*Nx[2]
        elif field == 'face_centered_B3':
            field_map[field]['ndata']=Nx[0]*Nx[1]*(Nx[2]+1)
                
        if sp[2]==b'int': dtype='i'
        elif sp[2]==b'float': dtype='f'
        elif sp[2]==b'double': dtype='d'
        field_map[field]['dtype']=dtype
        field_map[field]['dsize']=field_map[field]['ndata']*struct.calcsize(dtype)
        file.seek(field_map[field]['dsize'],1)
        offset = file.tell()
        tmp=file.readline()
        if len(tmp)>1: file.seek(offset)
        else: offset = file.tell()

    #grid['field_map'] = field_map
    #grid['data']={}
    return field_map

def read_field(file_pointer,field_map):
    ndata=field_map['ndata']
    dtype=field_map['dtype']
    file_pointer.seek(field_map['offset'])
    file_pointer.readline() # HEADER
    if field_map['read_table']: file_pointer.readline()
    data = file_pointer.read(field_map['dsize'])
    var = np.asarray(struct.unpack('>'+ndata*dtype,data))

    return var

def read_starvtk(starfile,time_out=False):
    file=open(starfile,'rb')
    star = {}
    star['filename']=starfile
    star['read_field'] = None
    star['read_type'] = None
    while star['read_field'] is None:
        star['data_offset']=file.tell()
        line = file.readline()
        parse_line(line, star)

    time=star['time']
    nstar=star['nstar']
    #print nstar
    fm=set_field_map(star)
    #print fm.keys()
    id=read_field(file,fm['star_particle_id'])
    mass=read_field(file,fm['star_particle_mass'])
    age=read_field(file,fm['star_particle_age'])
    pos=read_field(file,fm['star_particle_position']).reshape(nstar,3)
    vel=read_field(file,fm['star_particle_velocity']).reshape(nstar,3)
    file.close()
    star=[]
    for i in range(nstar):
        star.append({})

    for i in range(nstar):
        star_dict = star[i]
        star_dict['id']=id[i]
        star_dict['mass']=mass[i]
        star_dict['age']=age[i]
        star_dict['v1']=vel[i][0]
        star_dict['v2']=vel[i][1]
        star_dict['v3']=vel[i][2]
        star_dict['x1']=pos[i][0]
        star_dict['x2']=pos[i][1]
        star_dict['x3']=pos[i][2]
        star_dict['time']=time

    if time_out:
        return time,pd.DataFrame(star)
    else:
        return pd.DataFrame(star)

def starparvtk_to_starparhist(spfiles):
    import pandas as pd
    spfiles.sort()
    sphst={}
    for f in spfiles:
        sp=read_starvtk(f)
        if len(sp) != 0:
            idx=sp['mass'] != 0
        sp=sp[idx]
        sp.index=sp.time
        if len(sp) != 0:
            for id in sp.id:
                if id in sphst:
                    sphst[id]=pd.DataFrame.append(sphst[id],sp[sp.id == id])
                else:
                    sphst[id]=pd.DataFrame(sp[sp.id == id])
    sppd=pd.Panel(sphst)
    fsp=f.split('.')
    sppd.to_pickle(fsp[0]+'.starpar.p')
    return sppd
