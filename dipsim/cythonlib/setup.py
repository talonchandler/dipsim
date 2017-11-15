from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import cython_gsl

ext_modules=[ Extension("calc_efficiencies",
                        ["calc_efficiencies.pyx"],
                        libraries=["m"] + cython_gsl.get_libraries(),
                        library_dirs=[cython_gsl.get_library_dir()],
                        include_dirs=[cython_gsl.get_cython_include_dir()],
                        extra_compile_args = ["-ffast-math"])]

setup(
  name = "calc_efficiencies",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules)
