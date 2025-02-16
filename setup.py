import setuptools
import platform
import os
import sys
import subprocess

from numpy import get_include
from Cython.Build import cythonize

# Set to true to enable line profiling
line_debug = False

#Load the readme as a long description
with open('README.md') as fp:
    long_description = fp.read()

# Compile Cython
try:
    numpyInclude = [get_include(), '.']

    #Install Bioscrape Core Package
    package_data = {'bioscrape': ['*.pxd', '*.pyx']}
    bioscrape_src_dir = 'bioscrape'

    ext_options = {}
    ext_options['language'] = 'c++'
    ext_options['include_dirs'] = numpyInclude
    if platform.system() == "Darwin":
        ext_options['extra_compile_args'] = ['-std=c++11', "-mmacosx-version-min=10.9"]
        ext_options['extra_link_args'] = ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
        print('Using macOS clang args')
    if line_debug:
        ext_options['define_macros'] = [('CYTHON_TRACE_NOGIL', '1')]

    #used to generate HTML annotations of the cython code for
    #optimization purposes.
    cythonize_options = {
    "include_path":[bioscrape_src_dir],
    "language_level":"2" #Language level 3 does not work yet
    }
    if "annotate" in sys.argv:
        cythonize_options['annotate'] = True
        sys.argv.remove("annotate")
    if line_debug:
        cythonize_options["compiler_directives"] = {
            'profile': True,
            'linetrace': True,
            'binding': True
        }

    # # Turn on to enable gdb debugging
    # cythonize_options["gdb_debug"] = True

    #Determine if we install bioscrape, lineage, or both
    install_bioscrape = False
    install_lineage = False
    if "bioscrape" not in sys.argv and "lineage" not in sys.argv:
        install_bioscrape = True
        install_lineage = False 
    if "bioscrape" in sys.argv:
        install_bioscrape = True
        sys.argv.remove("bioscrape")
    if "lineage" in sys.argv:
        install_lineage = True
        sys.argv.remove("lineage")

    # elif "bioscrape" not in sys.argv:
    #     install_lineage = False

    cython_extensions = []
    if install_bioscrape:
        print("Installing Bioscrape...")
        bioscrape_source_files = ['random.pyx',
                                  'types.pyx',
                                  'simulator.pyx',
                                  'inference.pyx']
        bioscrape_extensions = [
                setuptools.Extension(
                    name = 'bioscrape.'+s.split('.')[0],
                    sources = [bioscrape_src_dir+'/'+s],
                    **ext_options) for s in bioscrape_source_files
            ]
        cython_extensions += cythonize(bioscrape_extensions, **cythonize_options)
        print("Bioscrape Cythonized.")

    if install_lineage:
        package_data['lineage'] = ['*.pxd', '*.pyx']
        print("Installing Lineage...")
        lineage_src_dir = 'lineage'
        lineage_source_files = ['lineage.pyx']
        lineage_extensions = [
            setuptools.Extension(name = 'bioscrape.'+s.split('.')[0],
                sources = [lineage_src_dir+'/'+s],
                **ext_options) for s in lineage_source_files
        ]
        cython_extensions += cythonize(lineage_extensions, **cythonize_options)
        print("Lineage Cythonized.")

except Exception as e:
    print("Error occured during Cython Compilation. Check C++ Compiler and Cython Installation.")
    raise

setuptools.setup(
    long_description=long_description,
    # package_dir = {'bioscrape' : bioscrape_src_dir},
    package_data = package_data,
    ext_modules = cython_extensions,
    zip_safe=False,
)
