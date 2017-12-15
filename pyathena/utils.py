import numpy as np
import pandas as pd
import os

def cc_arr(domain):
    le=domain['left_edge']
    re=domain['right_edge']
    dx=domain['dx']
    x=[]
    for i in range(3):
        x.append(np.arange(le[i],re[i],dx[i])+0.5*dx[i])

    return x

def fc_arr(domain):
    le=domain['left_edge']
    re=domain['right_edge']
    dx=domain['dx']
    x=[]
    for i in range(3):
        x.append(np.arange(le[i],re[i]+dx[i],dx[i]))

    return x

def cc_pos(domain,idx):
    le=domain['left_edge']
    dx=domain['dx']
    return le+0.5*dx+dx*np.array(idx)

def cc_idx(domain,pos):
    le=domain['left_edge']
    dx=domain['dx']
    if np.array(pos).ndim == 2:
        le=le[:,np.newaxis]
        dx=dx[:,np.newaxis]
    elif np.array(pos).ndim == 3:
        le=le[:,np.newaxis,np.newaxis]
        dx=dx[:,np.newaxis,np.newaxis]

    idx=(pos-le-0.5*dx)/dx
    return idx

def vecpot(bx,by):
    nx=bx.data.shape[0]
    ny=bx.data.shape[1]
    dx=(bx.bound[1]-bx.bound[0])/nx
    dy=(bx.bound[3]-bx.bound[2])/ny
    vecp=np.empty(bx.data.shape, dtype='float32')
    vecp[0,0] = 0.0
    for j in np.arange(1,ny): vecp[j,0] = vecp[j-1,0] + dy*bx.data[j,0]
    for i in np.arange(1,nx): vecp[:,i] = vecp[:,i-1] - dx*by.data[:,i]

    return vecp

def fline(bx,by,starts,nitermax=100,ds=1.0):
    nx=bx.data.shape[0]
    ny=bx.data.shape[1]
    xmax=bx.bound[1]
    xmin=bx.bound[0]
    ymax=bx.bound[3]
    ymin=bx.bound[2]
    dx=(xmax-xmin)/nx
    dy=(ymax-ymin)/ny
    xarr=np.arange(xmin,xmax,dx)
    yarr=np.arange(ymin,ymax,dy)

    flinesx=[]
    flinesy=[]
    for s in starts:
        flinex=[s[0]]
        fliney=[s[1]]
        x=s[0]
        y=s[1]
        niter=0
        while ((x-xmin)*(x-xmax) < 0) and ((y-ymin)*(x-ymax) < 0):
            xind=int((x-xmin)/dx)
            yind=int((y-ymin)/dy)
            #print x,y,xind,yind
            if xind+2 > nx or yind+2 > ny: break
            bxint=bilinear(bx.data[yind:yind+2,xind:xind+2],x-xarr[xind],y-yarr[yind],dx,dy)
            byint=bilinear(by.data[yind:yind+2,xind:xind+2],x-xarr[xind],y-yarr[yind],dx,dy)
            if np.abs(bxint)>np.abs(byint):
                xnext = x+ds
                ynext = y+byint/bxint*ds
            else:
                ynext = y+ds
                xnext = x+bxint/byint*ds
            flinex.append(xnext)
            fliney.append(ynext)
            x=xnext
            y=ynext
            niter += 1
            if niter > nitermax: break
        flinesx.append(flinex)
        flinesy.append(fliney)
    return flinesx,flinesy

def bilinear(f,x,y,dx,dy):
    f1=f[0,0]*(1-x/dx)+f[0,1]*x/dx
    f2=f[1,0]*(1-x/dx)+f[1,1]*x/dx
    fint = f1*(1-y/dy)+f2*y/dy

    return fint

def trilinear(f,x,y,z,dx,dy,dz):
    f1=f[0,0]*(1-x/dx)+f[0,1]*x/dx
    f2=f[1,0]*(1-x/dx)+f[1,1]*x/dx
    fint = f1*(1-y/dy)+f2*y/dy

    return fint

def plane_array(le,re,nx,ny,z=0.0):
    xarr=np.linspace(le[0],re[0],num=nx)
    yarr=np.linspace(le[1],re[1],num=ny)
    arr=np.empty((nx*ny,3),dtype='float32')
    for i,x in enumerate(xarr):
        for j,y in enumerate(yarr):
            arr[j+i*ny,0]=x
            arr[j+i*ny,1]=y
            arr[j+i*ny,2]=z

    return arr

def draw_cube(ax,le,re,color='b'):
    xr=[le[0],re[0]]
    yr=[le[1],re[1]]
    zr=[le[2],re[2]]
    dx=np.array(re)-np.array(le)
    plist=list(product(xr,yr,zr))
    pset = combinations(np.array(plist),2)
    for p1,p2 in pset:
        if (np.sum(np.abs(p1-p2)) == dx).any():
            ax.plot3D(*list(zip(p1,p2)),color=color)

def dump_visit(region,ds,ni=0,ne=1,fileout=True):
    visit_fname=ds.dir+ds.id+'.region.visit'
    if fileout:
        fp=open(visit_fname,'w')
        print(visit_fname)
        print('!NBLOCKS',region.ngrid)
        fp.write('!NBLOCKS %i' % region.ngrid)
        fp.write('\n')
    for i in range(ni,ne):
        tstring='.%4.4i.vtk' % i
        for g in region.grid_list:
            fname=g['filename'][:-9]+tstring
            if fileout:
                fp.write(fname)
                fp.write('\n')
            else:
                print(fname)

def gradient(phi,dx):
    Nx=phi.shape

    g1=np.empty(Nx)
    g2=np.empty(Nx)
    g3=np.empty(Nx)

    g1[:,:,1:-1]=(phi[:,:,2:]-phi[:,:,:-2])/dx[0]/2.0
    g1[:,:,0 ]=(phi[:,:,1 ]-phi[:,:,0 ])/dx[0]
    g1[:,:,-1]=(phi[:,:,-1]-phi[:,:,-2])/dx[0]

    g2[:,1:-1,:]=(phi[:,2:,:]-phi[:,:-2,:])/dx[1]/2.0
    g2[:,0 ,:]=(phi[:,1 ,:]-phi[:,0 ,:])/dx[1]
    g2[:,-1,:]=(phi[:,-1,:]-phi[:,-2,:])/dx[1]

    g3[1:-1,:,:]=(phi[2:,:,:]-phi[:-2,:,:])/dx[2]/2.0
    g3[0 ,:,:]=(phi[1 ,:,:]-phi[0 ,:,:])/dx[2]
    g3[-1,:,:]=(phi[-1,:,:]-phi[-2,:,:])/dx[2]

    return g1,g2,g3

def texteffect(fontsize=12):
  try:
    from matplotlib.patheffects import withStroke
    myeffect = withStroke(foreground="w", linewidth=3)
    kwargs = dict(path_effects=[myeffect], fontsize=fontsize)
  except ImportError:
    kwargs = dict(fontsize=fontsize)
  return kwargs

def pos3d(domain,x0=0.,y0=0.,z0=0.):
    x=cc_arr(domain)
    x1=x[0]
    x2=x[1]
    x3=x[2]

    nx1=len(x1)
    nx2=len(x2)
    nx3=len(x3)

    x3d = np.tile(x1.reshape(1,1,nx1),(nx3,nx2,1))-x0
    y3d = np.tile(x2.reshape(1,nx2,1),(nx3,1,nx1))-y0
    z3d = np.tile(x3.reshape(nx3,1,1),(1,nx2,nx1))-z0

    r3d = np.sqrt((x3d)**2+(y3d)**2+(z3d)**2)

    return r3d,x3d,y3d,z3d

def compare_files(source, output):
    if os.path.isfile(source):
        smtime=os.path.getmtime(source)
    else:
        return True
    if os.path.isfile(output):
        omtime=os.path.getmtime(output)
        if omtime < smtime:
            return False
        else:
            return True
    else:
        return False
