Call

>> python3 setup.py build_ext --inplace

to "cythonize" the .pyx into a .c using the options in setup.py.

Call

>> cython *.pyx -a

to generate a .html file that shows bottlenecks. 