from distutils.core import setup, Extension
import numpy as np

ext1 = Extension('_wordcount', sources = ['wordcountmodule.c'], 
    include_dirs = [np.get_include()],
    extra_compile_args=["-O3"])

setup (name = 'wordcount',
    version = '0.1',
    description = 'Word counting/bias computing extension',
    author = 'Alexander Solovyov, Petr Sulc, Benjamin Greenbaum',
    scripts = ['get_complexity.py', 'compute_sliding_window_force.py'],
    py_modules = ['wordcount'],
    ext_modules = [ext1])
