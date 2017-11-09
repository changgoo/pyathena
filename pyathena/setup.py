from distutils.core import setup
import numpy as np

from distutils.extension import Extension
from Cython.Build import cythonize

# We only have one extension for now:
extensions = [
    Extension("pyathena.synthetic_observations.get_los",
              ["synthetic_observations/get_los.pyx"],
              include_dirs=[np.get_include()])
]

setup(
    name="get_los",
    ext_modules=cythonize(extensions)
)
