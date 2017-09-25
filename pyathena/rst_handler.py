import struct
import numpy as np
import glob
import os

#writer 

def parse_misc_info(rstfile):
    fp=open(rstfile,'rb')
    search_block=['par','time','data','star','user']
    start={}
    size={}
    start['par']=0
    iblock=0

    while 1:
        block=search_block[iblock]
        size[block]=fp.tell()-start[block]
        
        l=fp.readline()
        if not l: break
    
        if l.startswith('N_STEP') or l.startswith('DENSITY') or \
           l.startswith('STAR') or l.startswith('USER'): 
            iblock+=1
            start[search_block[iblock]]=start[block]+size[block]

    data={}
    search_block=['par','time','star','user']
    for block in search_block:
        #fp.seek(start[block])
        #print fp.readline()
        fp.seek(start[block])
        data[block]=fp.read(size[block])

    fp.close()
    
    return data

def write_onefile(newfile,data_part,data_par):

    fp=open(newfile,'wb')
    fields=['DENSITY', '1-MOMENTUM', '2-MOMENTUM', '3-MOMENTUM', 'ENERGY','POTENTIAL',
        '1-FIELD', '2-FIELD', '3-FIELD', 'SCALAR 0','SCALAR 1','SCALAR 2']
    for block in ['par','time']: fp.write(data_par[block])

    fp.write('DENSITY\n')
    fp.write(data_part['DENSITY'].flatten().tobytes('C'))
    for f in fields[1:]:
        if f in data_part.keys():
        #print f,data_part[f].shape
            fp.write('\n%s\n' % f)
            fp.write(data_part[f].flatten().tobytes('C'))
    fp.write('\n')
    for block in ['star','user']: fp.write(data_par[block])
    fp.close()

    return

def write_allfile(pardata,rstdata,grids,id='newrst',dname='/tigress/changgoo/rst/',itime=0,verbose=False,scalar=0):
    ngrids=len(grids)
#    if not (ds.domain['Nx'][::-1] == rstdata['DENSITY'].shape).all():
#       print 'mismatch in DIMENSIONS!!'
#       print 'restart data dimension:', rstdata['DENSITY'].shape
#       print 'new grid data dimension:', ds.domain['Nx'][::-1] 
#
#       return -1

    fields = rstdata.keys()

    cc_varnames=['DENSITY','1-MOMENTUM','2-MOMENTUM','3-MOMENTUM',\
                 'ENERGY','POTENTIAL']
    fc_varnames=['1-FIELD','2-FIELD','3-FIELD']

    for i in range(ngrids):
        if i == 0:
          fname=id+'.%4.4d.rst' % itime
        else:
          fname=id+'-id%d.%4.4d.rst' % (i,itime)

        g=grids[i]
        gis=g['is']
        gnx=g['Nx']
        gie=gis+gnx

        #print g['filename'].split('.'),gis,gie,fname
        data={}
        for f in cc_varnames:
            if f in fields:
                data[f]=rstdata[f][gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]
 
        for f in fc_varnames:
            ib,jb,kb=(0,0,0)
            if f in fields:
                if f.startswith('1'): ib=1
                if f.startswith('2'): jb=1
                if f.startswith('3'): kb=1
                data[f]=rstdata[f][gis[2]:gie[2]+kb,gis[1]:gie[1]+jb,gis[0]:gie[0]+ib]
 
        for ns in range(scalar):
            f='SCALAR %d' % ns
            if f in fields:
                data[f]=rstdata[f][gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]
        write_onefile(dname+fname,data,pardata)

    return

def get_eint(rstdata,neg_correct=True):
    eint=rstdata['ENERGY'].copy()
    eint -= 0.5*rstdata['1-MOMENTUM']**2/rstdata['DENSITY']
    eint -= 0.5*rstdata['2-MOMENTUM']**2/rstdata['DENSITY']
    eint -= 0.5*rstdata['3-MOMENTUM']**2/rstdata['DENSITY']
    
    for i,f in enumerate(['1-FIELD','2-FIELD','3-FIELD']):
        if f is '1-FIELD': Bc=0.5*(rstdata[f][:,:,:-1]+rstdata[f][:,:,1:])
        elif f is '2-FIELD': Bc=0.5*(rstdata[f][:,:-1,:]+rstdata[f][:,1:,:])
        elif f is '3-FIELD': Bc=0.5*(rstdata[f][:-1,:,:]+rstdata[f][1:,:,:])
        eint -= 0.5*Bc**2
    
    if neg_correct:
        k,j,i=np.where(eint<0)
        eavg=[]
        for kk,jj,ii in zip(k,j,i):
            epart=eint[kk-1:kk+2,jj-1:jj+2,ii-1:ii+2]
            e_neg=epart[epart<0]
            Nneg=len(e_neg)
            eavg.append((epart.sum()-e_neg.sum())/(epart.size-e_neg.size))
            print eint[kk,jj,ii],eavg[-1],epart.sum(),e_neg.sum()
        eint[k,j,i]=np.array(eavg)
 
    return eint

def to_etot(rstdata):
    eint=rstdata['ENERGY'].copy()
   
    eint += 0.5*rstdata['1-MOMENTUM']**2/rstdata['DENSITY']
    eint += 0.5*rstdata['2-MOMENTUM']**2/rstdata['DENSITY']
    eint += 0.5*rstdata['3-MOMENTUM']**2/rstdata['DENSITY']
    
    for i,f in enumerate(['1-FIELD','2-FIELD','3-FIELD']):
        if f is '1-FIELD': Bc=0.5*(rstdata[f][:,:,:-1]+rstdata[f][:,:,1:])
        elif f is '2-FIELD': Bc=0.5*(rstdata[f][:,:-1,:]+rstdata[f][:,1:,:])
        elif f is '3-FIELD': Bc=0.5*(rstdata[f][:-1,:,:]+rstdata[f][1:,:,:])
        eint += 0.5*Bc**2
    return eint

def degrade(rstdata,scalar=0):
    
    cc_varnames=['DENSITY','1-MOMENTUM','2-MOMENTUM','3-MOMENTUM',\
                 'ENERGY','POTENTIAL']
    fc_varnames=['1-FIELD','2-FIELD','3-FIELD']
    scalar_varnames=[]
    for ns in range(scalar):
        scalar_varnames.append('SCALAR %d' % ns)
    rstdata_new={}
    for f in cc_varnames:
        if f is 'ENERGY':
            data=get_eint(rstdata)
        else:
            data=rstdata[f].copy()
        shape=np.array(data.shape)/2
        newdata=np.zeros(shape,dtype='d')
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    newdata += data[k::2,j::2,i::2]
        rstdata_new[f]=newdata*0.125
        
    for f in fc_varnames:
        data=rstdata[f].copy()
        shape=np.array(data.shape)/2
        if f is '1-FIELD':
            newdata=np.zeros(shape+np.array([0,0,1]),dtype='d')
            for j in range(2):
                for k in range(2):
                    newdata += data[k::2,j::2,::2]
        if f is '2-FIELD':
            newdata=np.zeros(shape+np.array([0,1,0]),dtype='d')
            for i in range(2):
                for k in range(2):
                    newdata += data[k::2,::2,i::2]
        if f is '3-FIELD':
            newdata=np.zeros(shape+np.array([1,0,0]),dtype='d')
            for j in range(2):
                for i in range(2):
                    newdata += data[::2,j::2,i::2]
        rstdata_new[f]=newdata*0.25
        
    rstdata_new['ENERGY']=to_etot(rstdata_new)
    return rstdata_new

def refine(rstdata,scalar=0):
    
    cc_varnames=['DENSITY','1-MOMENTUM','2-MOMENTUM','3-MOMENTUM',\
                 'ENERGY','POTENTIAL']
    fc_varnames=['1-FIELD','2-FIELD','3-FIELD']
    scalar_varnames=[]
    for ns in range(scalar):
        scalar_varnames.append('SCALAR %d' % ns)
    
    if scalar: cc_varnames += scalar_varnames
    rstdata_new={}
    for f in cc_varnames:
        if f is 'ENERGY':
            data=get_eint(rstdata)
        else:
            data=rstdata[f]
        shape=np.array(data.shape)*2
        newdata=np.zeros(shape,dtype='d')
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    newdata[k::2,j::2,i::2] = data.copy()
        rstdata_new[f]=newdata
        
    for f in fc_varnames:
        data=rstdata[f]
        shape=np.array(data.shape)*2
        if f is '1-FIELD':
            newdata=np.zeros(shape-np.array([0,0,1]),dtype='d')
            idata = 0.5*(data[:,:,:-1]+data[:,:,1:])

            for j in range(2):
                for k in range(2):
                    newdata[k::2,j::2,::2] = data.copy()
                    newdata[k::2,j::2,1::2] = idata.copy()

        if f is '2-FIELD':
            newdata=np.zeros(shape-np.array([0,1,0]),dtype='d')
            idata = 0.5*(data[:,:-1,:]+data[:,1:,:])
            for i in range(2):
                for k in range(2):
                    newdata[k::2,::2,i::2] = data.copy()
                    newdata[k::2,1::2,i::2] = idata.copy()
                    
        if f is '3-FIELD':
            newdata=np.zeros(shape-np.array([1,0,0]),dtype='d')
            idata = 0.5*(data[:-1,:,:]+data[1:,:,:])
            for j in range(2):
                for i in range(2):
                    newdata[::2,j::2,i::2] = data.copy()
                    newdata[1::2,j::2,i::2] = idata.copy()
        rstdata_new[f]=newdata
        
    rstdata_new['ENERGY']=to_etot(rstdata_new)
    return rstdata_new

def calculate_grid(Nx,NBx):
    NGrids=np.array(Nx)/np.array(NBx)
    NProcs=NGrids[0]*NGrids[1]*NGrids[2]
    grids=[]
    i=0
    print Nx, NBx, NGrids, NProcs
    for n in range(NGrids[2]):
       for m in range(NGrids[1]):
           for l in range(NGrids[0]):
               grid={}
               grid['id']=i
               grid['is']=np.array([l*NBx[0],m*NBx[1],n*NBx[2]])
               grid['Nx']=np.array(NBx)
               grids.append(grid)
               i += 1 

    return grids,NGrids

# reader

def parse_par(rstfile):

    fp=open(rstfile,'rb')
    par={}
    line=fp.readline()
    
    while 1:

        if line.startswith('<'):
            block=line[1:line.rfind('>')]
            if block == 'par_end': break
            par[block]={}
        line=fp.readline()

        if block in ['problem','domain1','time']:
            sp = line.strip().split()
            if len(sp) >= 3: par[block][sp[0]]=eval(sp[2])
        else:
            sp=line.split('=')
            if len(sp) == 2: par[block][sp[0].strip()]=sp[1].split('#')[0].strip()

    par[block]=fp.tell()

    fp.close()
    
    return par

def parse_rst(var,par,fm):
    
    starpar=False
    if par['configure'].has_key('star particles'): 
        if par['configure']['star particles'] == 'none':
            starpar=False
        else:
            starpar=True
    vtype='param'
    cc_varnames=['DENSITY','1-MOMENTUM','2-MOMENTUM','3-MOMENTUM','ENERGY','POTENTIAL']
    fc_varnames=['1-FIELD','2-FIELD','3-FIELD']
    dm=par['domain1']
    nx1=dm['Nx1']/dm['NGrid_x1']
    nx2=dm['Nx2']/dm['NGrid_x2']
    nx3=dm['Nx3']/dm['NGrid_x3']

    if var=='N_STEP':
        ndata=1
        dtype='i'
    elif var=='TIME':
        ndata=1
        dtype='d'
    elif var=='TIME_STEP':
        ndata=1
        if starpar: ndata+=1
        dtype='d'
    elif var in cc_varnames:
        ndata=nx1*nx2*nx3
        dtype='d'
        vtype='ccvar'
    elif var in fc_varnames:
        if var.startswith('1'): nx1 += 1
        if var.startswith('2'): nx2 += 1
        if var.startswith('3'): nx3 += 1
            
        ndata=nx1*nx2*nx3
        dtype='d'
        vtype='fcvar'
    elif var.startswith('SCALAR'):
        ndata=nx1*nx2*nx3
        dtype='d'
        vtype='ccvar'
    elif var.startswith('STAR PARTICLE LIST'):
        ndata=1
        dtype='i'
        vtype='star'
    else:
        return 0

    fm[var]={}
    
    fm[var]['ndata']=ndata
    fm[var]['dtype']=dtype
    fm[var]['vtype']=vtype
    
    if vtype == 'ccvar' or vtype == 'fcvar':
        fm[var]['nx']=(nx3,nx2,nx1)
        
    return 1

def read_star(fp):
# This works for MST_4pc
#    ivars=['id','merge_history','isnew','active']
#    dvars=['m','x1','x2','x3','v1','v2','v3','age','mage','mdot',\
#           'x1_old','x2_old','x3_old',\
#           'm_old','M1_old','M2_old','M3_old',\
#           'navg','n2avg','v1avg','v2avg','v3avg',\
#           'eavg','Vol','radius','SFUV','SNRate',\
#'SNprob',\
#           'x1sn','x2sn','x3sn',\
#          ]
# Latest restart file
    ivars=['id','merge_history','isnew','active']
    dvars=['m','x1','x2','x3','v1','v2','v3','age','mage','mdot',\
           'x1_old','x2_old','x3_old',\
          ]

    star_dict={}
    dtype='i'
    for var in ivars:
       data=fp.read(struct.calcsize(dtype))
       star_dict[var],=struct.unpack('<'+dtype,data)

    dtype='d'
    for var in dvars:
       data=fp.read(struct.calcsize(dtype))
       star_dict[var],=struct.unpack('<'+dtype,data)

    return star_dict

def read_rst_grid(rstfile,verbose=False):
    
    par=parse_par(rstfile)

    fp=open(rstfile,'rb')
    fp.seek(par['par_end'])
    rst={}
    data_array={}
    while 1:
        l=fp.readline()
        var=l.strip()

        if parse_rst(var,par,rst):
            dtype=rst[var]['dtype']
            ndata=rst[var]['ndata']
            vtype=rst[var]['vtype']
            dsize=ndata*struct.calcsize(dtype)
            data=fp.read(dsize)
            if vtype == 'param': 
                if verbose: print var,struct.unpack('<'+ndata*dtype,data)
            elif vtype == 'star':
                nstar,=struct.unpack('<'+ndata*dtype,data)
                data=fp.read(dsize)
                star_list=[] 
                if nstar > 0:
                  for i in range(nstar):
                      star_list.append(read_star(fp))
                  if verbose: 
                      print var, nstar
                      print star_list[0]
                      print star_list[nstar-1]
                data_array[var]=star_list
            else: 
                arr=np.asarray(struct.unpack('<'+ndata*dtype,data))
                arr.shape = rst[var]['nx']
                data_array[var]=arr
                if verbose: print var, arr.mean(), arr.shape
            fp.readline()
        else: 
            break
    if verbose: print l, fp.tell()
    fp.close()

    return rst,data_array

def read(rstfile,grids,NGrids,parfile=None,verbose=False):
    if parfile==None: par=parse_par(rstfile)
    else: par=parse_par(parfile)
    nprocs=len(grids)#par['domain1']['AutoWithNProc']
    field_maps=[]
    rstdata={}
    nx=NGrids*grids[0]['Nx']
    nx=nx[::-1]
    #nx=ds.domain['Nx'][::-1]
    print nx,nprocs
    dirname=os.path.dirname(rstfile)
    basename=os.path.basename(rstfile)

    fm,data=read_rst_grid(rstfile,verbose=verbose)

    g=grids[0]
    gis=g['is']
    gnx=g['Nx']
    gie=gis+gnx

    print fm['DENSITY']['nx'],gnx


    for k in fm:
        ib,jb,kb=(0,0,0)
        if fm[k]['vtype'] == 'ccvar':
            rstdata[k]=np.empty(nx,dtype=fm[k]['dtype'])
            rstdata[k][gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]=data[k]
        elif fm[k]['vtype'] == 'fcvar':
            if k.startswith('1'): ib=1
            if k.startswith('2'): jb=1
            if k.startswith('3'): kb=1
            rstdata[k]=np.empty((nx[0]+kb,nx[1]+jb,nx[2]+ib),dtype=fm[k]['dtype'])
            rstdata[k][gis[2]:gie[2]+kb,gis[1]:gie[1]+jb,gis[0]:gie[0]+ib]=data[k]
#for i in range(nprocs):
    for i in range(1,nprocs):
        g=grids[i]
        gis=g['is']
        gnx=g['Nx']
        gie=gis+gnx
#        if i % 50 == 0: 
#            print i,gis,gie
#            print rstfile,g['filename']
        rstfname = '%s/%s-id%d%s' % (dirname,basename[:-9],i,basename[-9:])
        if not os.path.isfile(rstfname):
            rstfname = '%s/../id%d/%s-id%d%s' % (dirname,i,basename[:-9],i,basename[-9:])
        fm,data=read_rst_grid(rstfname)

        if verbose > 1: print i,fm['DENSITY']['nx'],gnx

        for k in fm:
            ib,jb,kb=(0,0,0)
            if fm[k]['vtype'] == 'ccvar':
                rstdata[k][gis[2]:gie[2],gis[1]:gie[1],gis[0]:gie[0]]=data[k]
            elif fm[k]['vtype'] == 'fcvar':
                if k.startswith('1'): ib=1
                if k.startswith('2'): jb=1
                if k.startswith('3'): kb=1
                rstdata[k][gis[2]:gie[2]+kb,gis[1]:gie[1]+jb,gis[0]:gie[0]+ib]=data[k]

    return rstdata

def set_xpos_with_dm(dm):
    le=np.array([dm['x1min'],dm['x2min'],dm['x3min']])
    re=np.array([dm['x1max'],dm['x2max'],dm['x3max']])
    Lx=re-le
    Nx=np.array([dm['Nx1'],dm['Nx2'],dm['Nx3']])
    dx=Lx/Nx
    xc={}
    xf={}
    for i,ax in zip(range(3),['x','y','z']):
        xf[ax]=np.arange(le[i],re[i]+dx[i],dx[i])
        xc[ax]=np.arange(le[i],re[i],dx[i])+0.5*dx[i]
    return xf,xc


def set_xpos(ds):
    le=ds.domain['left_edge']
    re=ds.domain['right_edge']
    dx=ds.domain['dx']
    xc={}
    xf={}
    for i,ax in zip(range(3),['x','y','z']):
        xf[ax]=np.arange(le[i],re[i]+dx[i],dx[i])
        xc[ax]=np.arange(le[i],re[i],dx[i])+0.5*dx[i]
    return xf,xc

def to_hdf5(h5file,rstdata,ds):
    import h5py

    Bx=rstdata['1-FIELD']
    By=rstdata['2-FIELD']
    Bz=rstdata['3-FIELD']
    xf,xc=set_xpos(ds)

    f=h5py.File(h5file,'a')
    for name in ['Bfields','cell_centered_coord','face_centered_coord']:
        if name in f.keys():
            grp=f[name]
        else:
            grp=f.create_group(name)
        print name

    grp=f['Bfields']
    for name,B in zip(['Bx','By','Bz'],[Bx,By,Bz]):
        if name in grp.keys():
            dset=grp[name]
        else:
            dset=grp.create_dataset(name,B.shape,data=B,dtype=B.dtype)

    for k in grp.keys():
        for i,ax in enumerate(['z','y','x']):
            grp[k].dims[i].label=ax

    bfield=f['Bfields']
    ccoord=f['cell_centered_coord']
    fcoord=f['face_centered_coord']
    for ax in ['x','y','z']:
        if ax in ccoord.keys():
            print ax
        else:
            ccoord[ax] = xc[ax]
        
        if ax in fcoord.keys():
            print ax
        else:
            fcoord[ax] = xf[ax]

    for b in bfield.keys():
        bax=b[-1]

        for i,ax in enumerate(['z','y','x']):
            if ax == bax:
                bfield[b].dims[i].attach_scale(fcoord[ax])
            else:
                bfield[b].dims[i].attach_scale(ccoord[ax])

    f.close()

def divB(rstdata):
    Bx=rstdata['1-FIELD']
    By=rstdata['2-FIELD']
    Bz=rstdata['3-FIELD']
    dBx=np.diff(Bx,axis=2)
    dBy=np.diff(By,axis=1)
    dBz=np.diff(Bz,axis=0)
    dB = dBx+dBy+dBz
    return dB
