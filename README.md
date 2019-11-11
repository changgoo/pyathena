# pyathena
python scripts and example notebooks for Athena-TIGRESS

for cythonized modules, please run

% python setup.py build_ext --inplace

Dependent packages are speicified in `tigpy.yml`

For TIGRESS users, following instruction may help:

```
$ module load anaconda3
$ conda env create --prefix /tigress/$USER/tigpy-env -f tigpy-env.yml
$ conda activate /tigress/$USER/tigpy-env
```
