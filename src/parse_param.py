def parse_par(rstfile):

    fp=open(rstfile,'rb')
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


    par[block]=fp.tell()

    fp.close()

    return par,blocks,fields
