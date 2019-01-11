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
src_dir = 'bioscrape'
extensions = [Extension('bioscrape.'+s.split('.')[0],[src_dir+'/'+s], **ext_options) for s in sourceFiles]

setup(
    name = 'bioscrape',
    packages = ['bioscrape'],
    package_dir = {'bioscrape' : src_dir},
    package_data = {'bioscrape': ['*.pxd']},
    ext_modules = cythonize(extensions),
    zip_safe=False
)

