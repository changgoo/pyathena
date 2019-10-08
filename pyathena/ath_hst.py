#!python
import numpy as np
import re
import os
"""
  handling history data files from Athena
"""

def test_pickle(filename):
    hstmtime=os.path.getmtime(filename)
    if os.path.isfile(filename+'.p'):
        pmtime=os.path.getmtime(filename+'.p')
    else:
        pmtime=0.
    if pmtime < hstmtime:
        return False
    else:
        return True

def get_volume(filename):
    file=open(filename,'r')
    line=file.readline()
    file.close()

    match=re.search('-?\d(\.\d+|)[Ee][+\-]\d\d',line)

    return eval(match.group(0))

def get_varlist(filename,snfile=False):

    base=os.path.basename(filename)
    split=re.split('\.',base)
    ext=split[-1]
    file=open(filename,'r')

    if not snfile: line=file.readline()
    header=file.readline()
    file.close()

    if not snfile:
        varlist=re.split("\[\d+]\=|\n",header)
        for i in range(len(varlist)): varlist[i]=re.sub("\s|\W","",varlist[i])
        varlist=varlist[1:-1]
    else:
        varlist=re.split("\,",header[1:])
        for i in range(len(varlist)): varlist[i]=re.sub("\s|\W","",varlist[i])

    return varlist

def read_w_pandas(filename,silent=False,snfile=False,write=True):
    varlist=get_varlist(filename,snfile=snfile)

    import pandas as pd
    if not snfile: nheader=3
    else: nheader=1
    if test_pickle(filename):
        hst = pd.read_pickle(filename+'.p')
        if not silent: print("Reading a history file:" + filename+'.p')
    else: 
        hst = pd.read_csv(filename,skiprows=nheader,names=varlist,
                            sep='\s+',engine='python',comment='#')
        if not silent: print("Reading a history file:" + filename)
        hst.to_pickle(filename+'.p') 

    return hst

def read(filename,sortable_key=False,silent=False):

    if not silent: print("Reading a history file:" + filename)
    base=os.path.basename(filename)
    split=re.split('\.',base)
    ext=split[-1]
    file=open(filename,'r')
    lines=file.readlines()
    file.close()

    data={}

    header = lines[1]

    if ext == 'hst':
        match=re.search('-?\d(\.\d+|)[Ee][+\-]\d\d',lines[0])
        data['vol']=eval(match.group(0))

    if header[0] != "#":
        split=re.findall('[+\-]?\d\.?\d*[Ee][+\-]\d\d?',header)
        nvar=len(split)
        varlist=[]
        for i in range(nvar):
            varlist.append('var%d' % i)
    else:
        varlist=re.split("\[\d+]\=|\n",header)
        for i in range(len(varlist)): varlist[i]=re.sub("\s|\W","",varlist[i])
        varlist=varlist[1:-1]
        nvar=len(varlist)

    if sortable_key:
        for i in range(nvar): 
            head= '%02d' % (i+1)
            varlist[i] = head+varlist[i] 

    for var in varlist:
        data[var] = []

    for line in lines:
        if ext == 'hst':
            split=re.findall('[+\-]?\d\.?\d*[Ee][+\-]\d\d?',line)
        elif ext == 'sn':
            split=line.split()
        if nvar == len(split):
            for var, value  in zip(varlist, split):
                data[var].append(eval(value))
    for var in varlist:
        data[var] = np.array(data[var])

    return data

def readsn(filename,sortable_key=False,silent=False):

    if not silent: print("Reading a sn file:" + filename)
    base=os.path.basename(filename)
    split=re.split('\.',base)
    ext=split[-1]
    file=open(filename,'r')
    lines=file.readlines()
    file.close()

    data={}

    header = lines[0]

    if header[0] != "#":
        split=re.findall('[+\-]?\d\.?\d*[Ee][+\-]\d\d?',header)
        nvar=len(split)
        varlist=[]
        for i in range(nvar):
            varlist.append('var')
    else:
        varlist=re.split("\,",header[1:])
        for i in range(len(varlist)): varlist[i]=re.sub("\s|\W","",varlist[i])
        nvar=len(varlist)

    if sortable_key:
        for i in range(nvar): 
            head= '%02d' % (i+1)
            varlist[i] = head+varlist[i] 

    for var in varlist:
        data[var] = []

    for line in lines[1:]:
        split=re.split('\s*',line[:-1])
        if nvar == len(split):
            for var, value  in zip(varlist, split):
                data[var].append(eval(value))

    for var in varlist:
        data[var] = np.array(data[var])

    return data
