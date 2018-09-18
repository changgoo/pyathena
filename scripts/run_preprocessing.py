
from pyathena import preprocessing
import glob

base='/u/ckim14/'
dirs=glob.glob('{}/*_newacc*/id0'.format(base))
ids=[]
for dd in dirs:
    ids.append(dd.split('/')[-2])
print ids
ids=['R2_1pc_newacc']
for problem_id in ids: preprocessing.doall(base,problem_id,do_pickling=True)
