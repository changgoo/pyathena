import periodictable as pt
import pandas as pd
import numpy as np

def int_to_roman(input):
    """ Convert an integer to a Roman numeral. """

    if not isinstance(input, type(1)):
        raise TypeError( "expected integer, got %s" % type(input))
    if not 0 < input < 4000:
        raise ValueError( "Argument must be between 1 and 3999")
    ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
    nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
    result = []
    for i in range(len(ints)):
        count = int(input / ints[i])
        result.append(nums[i] * count)
        input -= ints[i] * count
    return ''.join(result)

class GF12_table(object):
    def __init__(self,data_path='./data/'):
        self.data_path=data_path
        self.elements = self.read_ion_abundance_table_()
        self.ion_frac = self.read_ion_frac_table_()
        self.temp = self.ion_frac['temp']
        self.cie_cooling = {}
        self.cie_cooling_per_ion = {}

    def get_total_cie_cooling(self,elements=['H','He','C','N','O','Ne','Mg','Si','S','Fe']):
        total_cooling = np.zeros_like(self.temp)
        for ion_name in elements:
            element=self.elements.loc[ion_name]
            nion=element['number']+1
            A=element['abundance']
            total_cooling+=self.get_cie_cooling(ion_name)*A
            
        return total_cooling
            
    def get_cie_cooling(self,ion_name):
        nion=self.elements.loc[ion_name]['number']+1
        if not (ion_name in self.cie_cooling):
            cie_cooling_ion,cie_cooling=self.read_cie_cooling_table_(ion_name)
            self.cie_cooling[ion_name]=cie_cooling
            self.cie_cooling_per_ion[ion_name]=cie_cooling_ion
        return self.cie_cooling[ion_name]
    
    def read_ion_abundance_table_(self):
        ion_abundance_file=self.data_path+'apjs420150t2_ascii.txt'

        fp=open(ion_abundance_file,'r')
        lines=fp.readlines()[4:-2]
        fp.close()

        element={}
        for l in lines:
            sp=l.split('\t')
            ename=pt.elements.name(sp[1].lower())
            element[ename.symbol]={}
            element[ename.symbol]['abundance']=eval(sp[3].replace(' x 10^','e'))
            element[ename.symbol]['number']=ename.number
            element[ename.symbol]['mass']=ename.mass
            element[ename.symbol]['datafile']='{}datafile{}.txt'.format(self.data_path,sp[2])

        return pd.DataFrame(element).T.sort_values('number')
    
    def read_ion_frac_table_(self):
        ion_frac_file=self.data_path+'tab2.txt'
        #ion_frac=pd.read_table(ion_frac_file,skiprows=125,delimiter=' ',)
        fp=open(ion_frac_file,'r')
        lines=fp.readlines()
        fp.close()
        nlines=len(lines[125:])

        fields=['temp']
        ion_frac={}
        import re
        for l in lines[9:121]:
            ion_name=l.split()[-5].split('{')[0]+l[l.rfind('{')+1:max(l.rfind('{')+2,l.rfind('}')-1)].replace('+','1')
            fields.append(ion_name)
        for f in fields:
            ion_frac[f]=np.empty(nlines)
        for iline,l in enumerate(lines[125:]):
            sp=l.split()
            nfields = len(sp)
            if nfields == len(fields):
                for ifield in range(nfields):
                    ion_frac[fields[ifield]][iline]=sp[ifield]
        return ion_frac
    
    def read_cie_cooling_table_(self,ion_name):
        element=self.elements.loc[ion_name]
        nion=element['number']+1
        nskip=nion+12
        fp=open(element['datafile'],'r')
        lines=fp.readlines()
        fp.close()
        #print lines[nskip]
        temp=self.temp

        cie=[]
        for l in lines[nskip:]:
            cie.append(l.split())
        cie=np.array(cie).astype('float')

        cie_new=np.empty((len(temp),nion))
        for i in range(nion):
            cie_new[:,i]=np.interp(temp,cie[:,0],cie[:,i+1])
        total_cooling_new=np.interp(temp,cie[:,0],cie[:,-1])
        
        return cie_new,total_cooling_new

class GF12_NIE(GF12_table):
    def __init__(self,data_path='./data/'):
        self.data_path=data_path
        self.elements = self.read_ion_abundance_table_()
        self.ion_frac = self.read_ion_frac_table_()
        self.temp = self.ion_frac['temp']

    def read_ion_frac_table_(self):
        ion_frac_file=self.data_path+'tab3.txt'
        #ion_frac=pd.read_table(ion_frac_file,skiprows=125,delimiter=' ',)
        fp=open(ion_frac_file,'r')
        lines=fp.readlines()
        fp.close()
        nlines=len(lines[124:][::-1])

        fields=['temp']
        ion_frac={}
        import re
        for l in lines[9:121]:
            ion_name=l.split()[-5].split('{')[0]+l[l.rfind('{')+1:max(l.rfind('{')+2,l.rfind('}')-1)].replace('+','1')
            fields.append(ion_name)
        for f in fields:
            ion_frac[f]=np.empty(nlines)
        for iline,l in enumerate(lines[124:][::-1]):
            sp=l.split()
            nfields = len(sp)
            if nfields == len(fields):
                for ifield in range(nfields):
                    ion_frac[fields[ifield]][iline]=sp[ifield]
        return ion_frac
    

