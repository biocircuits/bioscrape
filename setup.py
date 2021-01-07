from distutils.core import setup
from Cython.Build import cythonize
from numpy import get_include
from distutils.extension import Extension
import platform
import os
import sys



#Load the readme as a long description
with open('README.rst') as fp:
    long_description = fp.read()

print("sys.argv", sys.argv)
# Compile Cython
try:
    numpyInclude = [get_include(), '.']

    ext_options = {}
    ext_options['language'] = 'c++'
    ext_options['include_dirs'] = numpyInclude
    if platform.system() == "Darwin":
        ext_options['extra_compile_args'] = ['-std=c++11', "-mmacosx-version-min=10.9"]
        ext_options['extra_link_args'] = ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
        print('Using macOS clang args')

    #used to generate HTML annotations of the cython code for
    #optimization purposes.
    if "annotate" in sys.argv:
        ext_options['annotate'] = True
        sys.argv.remove("annotate")

    #Determine if we install bioscrape
    install_bioscrape = False
    if "bioscrape" in sys.argv:
        install_bioscrape = True
        sys.argv.remove("bioscrape")
    elif "lineage" not in sys.argv:
        install_bioscrape = True

    #Determine if we install lineage
    install_lineage = False
    if "lineage" in sys.argv:
        install_lineage = True
        sys.argv.remove("lineage")
    elif "bioscrape" not in sys.argv:
        install_lineage = True

    #Install Bioscrape Core Package
    package_data = {'bioscrape': ['*.pxd']}
    bioscrape_src_dir = 'bioscrape'

    cython_extensions = []
    if install_bioscrape:
        print("Installing Bioscrape...")
        bioscrape_source_files = ['random.pyx', 'types.pyx', 'simulator.pyx', 'inference.pyx']
        bioscrape_extensions = [
                Extension('bioscrape.'+s.split('.')[0],[bioscrape_src_dir+'/'+s], **ext_options) 
                for s in bioscrape_source_files
            ]
        cython_extensions += cythonize(bioscrape_extensions)
        print("Bioscrape compiled.")

    if install_lineage:
        package_data['lineage'] = ['*.pxd']
        print("Installing Lineage...")
        lineage_src_dir = 'lineage'
        lineage_source_files = ['lineage.pyx']
        lineage_extensions = [
            Extension('bioscrape.'+s.split('.')[0],[lineage_src_dir+'/'+s], **ext_options) 
            for s in lineage_source_files
        ]
        cython_extensions += cythonize(lineage_extensions)
        print("Lineage compiled.")
        
except Exception as e:
    print("Error occured during Cython Compilation. Check C++ Compiler and Cython Installation.")
    raise

setup(
    name = 'bioscrape',
    version = '1.0.0',
    author='Biocircuits',
    url='https://github.com/biocircuits/bioscrape/',
    description='Biological Stochastic Simulation of Single Cell Reactions and Parameter Estimation.',
    long_description=long_description,
    packages = ['bioscrape'],
    package_dir = {'bioscrape' : bioscrape_src_dir},
    package_data = package_data,
    ext_modules = cython_extensions,
    zip_safe=False,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
    ],
    install_requires=[
        "matplotlib",
        "pytest",
        "numpy",
        "scipy",
        "Cython",
        "python-libsbml",
        "beautifulsoup4",
        "sympy",
        "emcee",
        "pandas"
    ],
    python_requires='>=3.6',
    keywords="SBML synthetic biology modeling Chemical Reaction Network CRN simulator stochastic parameter inference",
    tests_require=["pytest"],
    project_urls={
    'Documentation': 'https://readthedocs.org/projects/biocrnpyler/',
    'Funding': 'http://www.cds.caltech.edu/~murray/wiki/DARPA_BioCon',
    'Source': 'https://github.com/biocircuits/bioscrape/',
    'Tracker': 'https://github.com/biocircuits/bioscrape/issues',
    }
)