
from pyathena import preprocessing
import glob

base='/tigress/changgoo/'
dirs=glob.glob('{}/R*_rst_newacc/hst'.format(base))
ids=[]
for dd in dirs:
    ids.append(dd.split('/')[-2])
print ids
#ids=['R8_4pc_newacc']
for problem_id in ids: preprocessing.doall(base,problem_id,do_pickling=True)
