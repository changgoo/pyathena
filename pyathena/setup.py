from distutils.core import setup
from Cython.Build import cythonize
setup(
    ext_modules=cythonize(["vtk_reader.pyx","synthetic_observations.pyx"])
)
