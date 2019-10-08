import struct
import numpy as np
import glob
import os
import sys

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
    
        if l.startswith(b'N_STEP') or l.startswith(b'DENSITY') or \
           l.startswith(b'STAR') or l.startswith(b'USER'): 
            iblock+=1
            start[search_block[iblock]]=start[block]+size[block]

    data={}
    search_block=['par','time','star','user']
    for block in search_block:
        if block in start:
            fp.seek(start[block])
            data[block]=fp.read(size[block])

    fp.close()
    
    return data

def write_onefile(newfile,data_part,data_par):

    fp=open(newfile,'wb')
    fields=['DENSITY', '1-MOMENTUM', '2-MOMENTUM', '3-MOMENTUM', 'ENERGY','POTENTIAL',
            '1-FIELD', '2-FIELD', '3-FIELD',
            'SCALAR 0','SCALAR 1','SCALAR 2','SCALAR 3','SCALAR 4',
            'SCALAR 5','SCALAR 6','SCALAR 7','SCALAR 8','SCALAR 9']
    for block in ['par','time']: fp.write(data_par[block])

    fp.write(b'DENSITY\n')
    fp.write(data_part['DENSITY'].flatten().tobytes('C'))
    for f in fields[1:]:
        if f in list(data_part.keys()):
        #print f,data_part[f].shape
            fp.write('\n{}\n'.format(f).encode())
            fp.write(data_part[f].flatten().tobytes('C'))
    fp.write(b'\n')
    for block in ['star','user']: 
      if block in data_par: fp.write(data_par[block])
    fp.close()

    return

def write_allfile(pardata,rstdata,grids,grid_disp=np.array([0,0,0]),
  id='newrst',dname='/tigress/changgoo/rst/',itime=0,verbose=False,scalar=0):
    ngrids=len(grids)
#    if not (ds.domain['Nx'][::-1] == rstdata['DENSITY'].shape).all():
#       print 'mismatch in DIMENSIONS!!'
#       print 'restart data dimension:', rstdata['DENSITY'].shape
#       print 'new grid data dimension:', ds.domain['Nx'][::-1] 
#
#       return -1

    fields = list(rstdata.keys())

    cc_varnames=['DENSITY','1-MOMENTUM','2-MOMENTUM','3-MOMENTUM',\
                 'ENERGY','POTENTIAL']
    fc_varnames=['1-FIELD','2-FIELD','3-FIELD']

    for g in grids:
        i=g['id']
        if i == 0:
          fname=id+'.%4.4d.rst' % itime
        else:
          fname=id+'-id%d.%4.4d.rst' % (i,itime)

        gis=g['is']-grid_disp
        gnx=g['Nx']
        gie=gis+gnx

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
        if verbose: dname+fname
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
        k_end,j_end,i_end = eint.shape
        k_str=j_str=i_str = 0
        k,j,i=np.where(eint<0)
        eavg=[]
        for kk,jj,ii in zip(k,j,i):
            kl=kk if kk==k_str else kk-1
            kh=kk+1 if kk==(k_end-1) else kk+2
            jl=jj if jj==j_str else jj-1
            jh=jj+1 if jj==(j_end-1) else jj+2
            il=ii if ii==i_str else ii-1
            ih=ii+1 if ii==(i_end-1) else ii+2
            epart=eint[kl:kh,jl:jh,il:ih]
            e_neg=epart[epart<0]
            Nneg=len(e_neg)
            eavg.append((epart.sum()-e_neg.sum())/(epart.size-e_neg.size))
            print(kk,jj,ii,eint[kk,jj,ii],eavg[-1],epart.sum(),e_neg.sum())
        eint[k,j,i]=np.array(eavg)
        if len(eint[eint<0]) > 0: sys.exit("negative energy persist!")
 
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
    if scalar: cc_varnames += scalar_varnames

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
                 'ENERGY']
    if 'POTENTIAL' in rstdata: cc_varnames += ['POTENTIAL']
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
    NGrids=(np.array(Nx)/np.array(NBx)).astype('int')
    NProcs=NGrids[0]*NGrids[1]*NGrids[2]
    grids=[]
    i=0
    print('Domain Size:',Nx)
    print('Grid Size:', NBx)
    print('Processor configuration:', NGrids)
    print('Number of Processors:', NProcs)
    for n in range(NGrids[2]):
       for m in range(NGrids[1]):
           for l in range(NGrids[0]):
               grid={}
               grid['id']=i
               grid['is']=np.array([l*NBx[0],m*NBx[1],n*NBx[2]]).astype('int')
               grid['Nx']=np.array(NBx).astype('int')
               grids.append(grid)
               i += 1 

    return grids,NGrids

# reader

def parse_par(rstfile):

    fp=open(rstfile,'rb')
    par={}
    line=fp.readline().decode('utf-8')
    
    while 1:

        if line.startswith('<'):
            block=line[1:line.rfind('>')]
            if block == 'par_end': break
            par[block]={}
        line=fp.readline().decode('utf-8')

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
    if 'star particles' in par['configure']:
        if par['configure']['star particles'] == 'none':
            starpar=False
        else:
            starpar=True
    vtype='param'
    cc_varnames=['DENSITY','1-MOMENTUM','2-MOMENTUM','3-MOMENTUM','ENERGY','POTENTIAL']
    fc_varnames=['1-FIELD','2-FIELD','3-FIELD']
    dm=par['domain1']
    nx1=int(dm['Nx1']/dm['NGrid_x1'])
    nx2=int(dm['Nx2']/dm['NGrid_x2'])
    nx3=int(dm['Nx3']/dm['NGrid_x3'])

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

def read_star(fp,nscal=0,ghost=True):
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
# additional fields depending on the version
    for i in range(nscal):
        dvars += ['metal{}'.format(i)]

    if ghost:
        dvars += ['mghost','M1ghost','M2ghost','M3ghost']
        for i in range(nscal):
            dvars += ['Sghost{}'.format(i)]

    star_dict={}
    dtype='i'
    for var in ivars:
       data=fp.read(struct.calcsize(dtype))
       tmp=struct.unpack('<'+dtype,data)
       star_dict[var]=tmp

    dtype='d'
    for var in dvars:
       data=fp.read(struct.calcsize(dtype))
       tmp=struct.unpack('<'+dtype,data)
       #if var is 'm': print(var,tmp)
       star_dict[var]=tmp

    return star_dict

def read_rst_grid(rstfile,verbose=False,starghost=True):
    
    par=parse_par(rstfile)

    fp=open(rstfile,'rb')
    fp.seek(par['par_end'])
    rst={}
    data_array={}
    nscal=0
    while 1:
        l=fp.readline().decode('utf-8')
        var=l.strip()

        if parse_rst(var,par,rst):
            dtype=rst[var]['dtype']
            ndata=rst[var]['ndata']
            vtype=rst[var]['vtype']
            dsize=ndata*struct.calcsize(dtype)
            data=fp.read(dsize)
            if vtype == 'param': 
                if verbose: print(var,struct.unpack('<'+ndata*dtype,data))
            elif vtype == 'star':
                nstar,=struct.unpack('<'+ndata*dtype,data)
                data=fp.read(dsize)
                star_list=[] 
                if nstar > 0:
                  for i in range(nstar):
                      star_list.append(read_star(fp,nscal=nscal,ghost=starghost))
                  if verbose: 
                      print(var, nstar)
                      print(star_list[0])
                      print(star_list[nstar-1])
                data_array[var]=star_list
            else: 
                arr=np.asarray(struct.unpack('<'+ndata*dtype,data))
                arr.shape = rst[var]['nx']
                data_array[var]=arr
                if verbose: print(var, arr.mean(), arr.shape)
                if var.startswith('SCALAR'): nscal += 1
            fp.readline()
        else: 
            break
    if verbose: print(l, fp.tell())
    fp.close()

    return rst,data_array

def read(rstfile,grids,NGrids,parfile=None,verbose=False,starghost=True):
    if parfile==None: par=parse_par(rstfile)
    else: par=parse_par(parfile)
    nprocs=len(grids)#par['domain1']['AutoWithNProc']
    field_maps=[]
    rstdata={}
    nx=NGrids*grids[0]['Nx']
    nx=nx[::-1]
    #nx=ds.domain['Nx'][::-1]
    print(nx,nprocs)
    dirname=os.path.dirname(rstfile)
    basename=os.path.basename(rstfile)

    fm,data=read_rst_grid(rstfile,verbose=verbose,starghost=starghost)

    g=grids[0]
    gis=g['is']
    gnx=g['Nx']
    gie=gis+gnx

    print(fm['DENSITY']['nx'],gnx)


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
        fm,data=read_rst_grid(rstfname,starghost=starghost)

        if verbose > 1: print(i,fm['DENSITY']['nx'],gnx)

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

def read_part(rstfile,grids,nx,verbose=False):
    nprocs=len(grids)
    field_maps=[]
    rstdata={}
    print(nx,nprocs)

    basename=os.path.basename(rstfile)
    pid=basename[:-9]
    fm,data=read_rst_grid(rstfile,verbose=verbose)

    g=grids[0]
    gis=g['is']
    gnx=g['Nx']
    gie=gis+gnx
    ks=gis[2]

    print(fm['DENSITY']['nx'],gnx)


    for k in fm:
        ib,jb,kb=(0,0,0)
        if fm[k]['vtype'] == 'ccvar':
            rstdata[k]=np.empty(nx,dtype=fm[k]['dtype'])
        elif fm[k]['vtype'] == 'fcvar':
            if k.startswith('1'): ib=1
            if k.startswith('2'): jb=1
            if k.startswith('3'): kb=1
            rstdata[k]=np.empty((nx[0]+kb,nx[1]+jb,nx[2]+ib),dtype=fm[k]['dtype'])

    for i in range(nprocs):
        g=grids[i]
        gis=g['is']
        gnx=g['Nx']
        gie=gis+gnx
        gid=g['id']
        if gid > 0:
            rstfname = rstfile.replace('{}.'.format(pid),'{}-id{}.'.format(pid,gid))
        else:
            rstfname = rstfile
        if not os.path.isfile(rstfname):
            rstfname = rstfile.replace('id{}/{}.'.format(gid,pid),
                                       'id{}/{}-id{}.'.format(gid,pid,gid))

        fm,data=read_rst_grid(rstfname)

        if verbose > 1: print(i,fm['DENSITY']['nx'],gnx)

        for k in fm:
            ib,jb,kb=(0,0,0)
            if fm[k]['vtype'] == 'ccvar':
                rstdata[k][gis[2]-ks:gie[2]-ks,gis[1]:gie[1],gis[0]:gie[0]]=data[k]
            elif fm[k]['vtype'] == 'fcvar':
                if k.startswith('1'): ib=1
                if k.startswith('2'): jb=1
                if k.startswith('3'): kb=1
                rstdata[k][gis[2]-ks:gie[2]-ks+kb,gis[1]:gie[1]+jb,gis[0]:gie[0]+ib]=data[k]

    return rstdata


def set_xpos_with_dm(dm):
    le=np.array([dm['x1min'],dm['x2min'],dm['x3min']])
    re=np.array([dm['x1max'],dm['x2max'],dm['x3max']])
    Lx=re-le
    Nx=np.array([dm['Nx1'],dm['Nx2'],dm['Nx3']])
    dx=Lx/Nx
    xc={}
    xf={}
    for i,ax in zip(list(range(3)),['x','y','z']):
        xf[ax]=np.arange(le[i],re[i]+dx[i],dx[i])
        xc[ax]=np.arange(le[i],re[i],dx[i])+0.5*dx[i]
    return xf,xc


def set_xpos(ds):
    le=ds.domain['left_edge']
    re=ds.domain['right_edge']
    dx=ds.domain['dx']
    xc={}
    xf={}
    for i,ax in zip(list(range(3)),['x','y','z']):
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
        if name in list(f.keys()):
            grp=f[name]
        else:
            grp=f.create_group(name)
        print(name)

    grp=f['Bfields']
    for name,B in zip(['Bx','By','Bz'],[Bx,By,Bz]):
        if name in list(grp.keys()):
            dset=grp[name]
        else:
            dset=grp.create_dataset(name,B.shape,data=B,dtype=B.dtype)

    for k in list(grp.keys()):
        for i,ax in enumerate(['z','y','x']):
            grp[k].dims[i].label=ax

    bfield=f['Bfields']
    ccoord=f['cell_centered_coord']
    fcoord=f['face_centered_coord']
    for ax in ['x','y','z']:
        if ax in list(ccoord.keys()):
            print(ax)
        else:
            ccoord[ax] = xc[ax]
        
        if ax in list(fcoord.keys()):
            print(ax)
        else:
            fcoord[ax] = xf[ax]

    for b in list(bfield.keys()):
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
