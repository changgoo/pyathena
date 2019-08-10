def write_par_from_rst(rstfile,parfile):
    fp=open(rstfile,'rb')
    search_block='par'
    start=0
    first=True
    while 1:
        l=fp.readline()
        if not l: break
        if l.startswith(b'# --------------------- PAR_DUMP '):
            if first: 
                start=fp.tell()
                first=False
            else: 
                size=fp.tell()-start-len(l)
                break
        if l.startswith(b'<par_end>'):
            size=fp.tell()-start
            break

    fp.seek(start)
    data=fp.read(size)

    fp.close()

    fp=open(parfile,'wb')
    fp.write(data)
    fp.close()

def parse_par(rstfile):

    fp=open(rstfile,'r')
    par={}
    fields={}
    blocks=[]
    line=fp.readline()

    while 1:

        if line.startswith('<'):
            block=line[1:line.rfind('>')]
            if block == 'par_end': break
            par[block]={}
            fields[block]=[]
            blocks.append(block)
        line=fp.readline()
        if '' == line: break
        sp=line.split('=')
        if len(sp) >= 2:
            
            field=sp[0].strip()
            sp2="=".join(sp[1:]).split('#')
            value=sp2[0].strip()
            if len(sp2) == 2:
                comment=sp2[1].strip()
            else:
                comment=''
            par[block][field]=[value,comment]
            fields[block].append(field)


    par['par_end']=fp.tell()

    fp.close()

    return par,blocks,fields

def get_params(rstfile):
    par,blocks,fields = parse_par(rstfile)

    params={}
    param_blocks=['domain1','problem']
    if 'feedback' in blocks: param_blocks.append('feedback')
    for block in param_blocks:
        for key in par[block]:
            params[key]=float(par[block][key][0])
    if 'configure' in blocks:
        params['nscalars']=int(par['configure']['nscalars'][0])
    return params
