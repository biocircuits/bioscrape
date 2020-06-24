from distutils.core import setup
from Cython.Build import cythonize
from numpy import get_include
from distutils.extension import Extension
import platform
import os

numpyInclude = [get_include(), '.']

# cyrandom module: has random number generation stuff

sourceFiles = ['random.pyx', 'types.pyx', 'simulator.pyx', 'inference.pyx']


ext_options = {}
ext_options['language'] = 'c++'
ext_options['include_dirs'] = numpyInclude
if platform.system() == "Darwin":
    ext_options['extra_compile_args'] = ['-std=c++11', "-mmacosx-version-min=10.9"]
    ext_options['extra_link_args'] = ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
    print('Using macOS clang args')


# building part
src_dir = 'bioscrape'
extensions = [Extension('bioscrape.'+s.split('.')[0],[src_dir+'/'+s], **ext_options) for s in sourceFiles]

setup(
    name = 'bioscrape',
    packages = ['bioscrape'],
    package_dir = {'bioscrape' : src_dir},
    package_data = {'bioscrape': ['*.pxd']},
    ext_modules = cythonize(extensions),
    zip_safe=False,
    version = '1.0.0'
)

