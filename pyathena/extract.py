import cPickle as p
import glob

files=glob.glob('*.p')

for f in files:
    slc=p.load(open(f,'rb'))
    p.dump(slc['z']['nH'],open(f.replace('slice.p','nH.zslice.p'),'wb'),p.HIGHEST_PROTOCOL)
