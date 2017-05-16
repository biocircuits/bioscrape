from distutils.core import setup
from Cython.Build import cythonize
from numpy import get_include
from distutils.extension import Extension
import os

numpyInclude = [get_include(), '.']

# cyrandom module: has random number generation stuff

sourceFiles = ['random.pyx', 'types.pyx', 'simulator.pyx', 'inference.pyx']


ext_options = {}
ext_options['language'] = 'c++'
ext_options['include_dirs'] = numpyInclude
extra_compile_args = []


# building part

extensions = [Extension('bioscrape.'+s.split('.')[0],[s], **ext_options) for s in sourceFiles]

setup(
    name = 'bioscrape',
    packages = ['bioscrape'],
    package_dir = {'bioscrape' : '.'},
    ext_modules = cythonize(extensions)
)

