Installation
============

``pip install bioscrape`` is the easiest way to install the software and
all its dependencies. However, this requires that you have a C++
compiler already configured appropriately with your Python environment
(Bioscrape is compiled on installation in order to ensure compatability
across different operating systems and hardware). If pip fails, please
see the below information on setting up the dependencies in your system.

You can also try to run the package without installing through Google
Colab. We have provided three Python notebooks to get started, however,
there are many more in the examples folder that is provided with the
repository.

-  Getting started with Bioscrape: |Bioscrape Core|

-  Bioscrape analysis features: |Bioscrape Sensitivity Analysis|

-  Parameter inference with Bioscrape: |Bioscrape Inference|

To run the example notebooks in the Github repository, you can either
clone the repository by running
``git clone https://github.com/biocircuits/bioscrape/`` from your
terminal/the Anaconda command prompt on Windows or create an online
`Github Codespace <https://github.com/features/codespaces>`__ from the
repository. To run the example jupyter notebooks, you can open the open
the ``.ipynb`` files in the examples/ directory. Note that a `Jupyter
notebook installation <https://jupyter.org/install>`__ is required to
setup a Jupyter notebook kernel (with a Python version >=3.8) that you
can use to run the example notebooks.

Prereqs and Dependencies
------------------------

Required dependencies
~~~~~~~~~~~~~~~~~~~~~

C++ Compiler
^^^^^^^^^^^^

This package relies heavily on Cython, which is an extension of Python
that converts Cython code to C/C++ and then compiles it to make it
extremely fast. As such, you will need a C++ compiler and build system
with the ability to compile C++ to dynamic libraries. This is different
on different operating systems.

Linux
'''''

The standard compiler is g++. If you type g++ at the terminal, you
should get an error message about no input files, and not an error
message that says "command not found". If the command isn't found, you
can install g++ (along with a bunch of other things) with the following
terminal command on Ubuntu, which will require your password. If you are
not using Ubuntu, then there should be a similar way to install
software.

::

   sudo apt-get install build-essential

Mac OS
''''''

On Mac, if you type g++ at the terminal, you should get an error that
says something like "no input files". If you get an error saying
something like "command not found," you will need to install
`Xcode <https://developer.apple.com/xcode/downloads/>`__. This will give
you all the tools you need.

If you are trying to install and are getting errors that suggest some
files are not being found, you most likely need to install the Xcode
Developer Tools in addition to just Xxode itself. To do this, running
the following command in the terminal should be sufficient.

::

   xcode-select --install

Windows
'''''''

If you already have Anaconda installed (and don't have a C++ compiler
which already works with python), you might need to `uninstall
Anaconda <https://docs.anaconda.com/anaconda/install/uninstall/>`__.

The standard compiler on Windows is the one that comes with Visual C++.
You can get the compiler from
`here <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`__.
When installing visual C++, you need to make sure that the following
features are downloaded and installed by selecting a custom install and
making sure the following boxes are selected (adding additional features
is okay as well):

1. Visual C++ Build tools core features [typically selected
   automatically]
2. VC++ 20xx toolset (x86, x64) [option name varies slightly and 
   there may be multiple items with the same/similar name; 
   install the older versions if the latest one does not work].
3. Visual C++ 20xx Redistributable Update [latest version should work].
4. Windows 11 SDK (10.0.16299.0) for Desktop C++ [option name may vary
   slightly]
5. Desktop development C++

Note: if you have errors installing Bioscrape in the later steps, such
as "cannot find vcvarsall.bat" or linkage errors (compilation warnings
are OK) you can also try installing Visual Studio Community from `this
link <https://docs.microsoft.com/en-us/cpp/build/vscpp-step-0-installation?view=vs-2019>`__.

Python dependencies
^^^^^^^^^^^^^^^^^^^

WARNING: It is very important that the appropriate C++ compiler is
installed BEFORE the python dependencies.

The easy way
''''''''''''

The easy way to get all the required Python dependencies is to use
`Anaconda <https://www.continuum.io/downloads>`__. Make sure you
download the Python 3 version for your machine (64 bit). If you already
have Anaconda, then no need to download and install again.

Once you have installed Anaconda, running the following commands will
get you ready to install and use this package. numpy,scipy,matplotlib,
cython, sympy, pandas, and emcee are absolutely required to use the
software. libsbml, lxml, and beautifulsoup4 are required for
loading/saving bioscrape XML and SBML model files. Running the following
commands should get everything installed that you need.

To make sure your anaconda is updated, run the following two commands.
These are unnecessary if you have just installed anaconda.

::

   conda update conda
   conda update anaconda

Then, in order to get the required dependencies, run the following
commands.

::

   conda install -c anaconda lxml beautifulsoup4 numpy scipy cython jupyter matplotlib sympy pandas
   conda install -c astropy emcee
   conda install -c sbmlteam python-libsbml

NOTE: if you are running python 3.8, you may need to install a different
version of libsbml via:

::

   pip install python-libsbml

The hard way (NOT RECOMMENDED)
''''''''''''''''''''''''''''''

Please use anaconda to ensure that bioscrape installs properly.

If you really really don't like Anaconda for some reason, you can also
use your own Python distribution and just make sure you have updated
versions of numpy, scipy, matplotlib, sympy, and cython. Also, make sure
you have the emcee package. You will also want jupyter in order to view
and run the example notebooks. libsbml is required if you want to import
SBML files.

Setup Instructions
------------------

Once you have Anaconda installed and updated along with the appropriate
C++ compiler for your platform as outlined above, the installation
process is the same on all operating systems. It is detailed below.

To install the bioscrape package, first check it out or download and
extract it from GitHub
`here <https://github.com/biocircuits/bioscrape>`__. Go into the
bioscrape directory and run

::

   pip install .

The installation will cause a bunch of C++ to be compiled. This may
create a bunch of warnings and text. If there are no **errors**, you are
fine and have probably completed a successful installation and can
attempt to run the testing notebook.

Note : When reinstalling, if you have any bioscrape related Jupyter
notebook open or other bioscrape code running on your system already,
then the installation is likely to fail. In this case, simply close all
terminals (including Jupyter notebooks) that might be running bioscrape
and then try installing again.

Uninstalling
~~~~~~~~~~~~

If you want to uninstall, simply go to your site-packages folder for
your Python distribution and delete the bioscrape directory as well as
the bioscrape .egg file if one exists. This will uninstall the package.

Installing locally (not recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you're a Python expert and want to use the package locally, only run
the first command.

::

   python setup.py build_ext

This will build the files into a local subdirectory called
build/lib/bioscrape, where the lib might have a suffix with the details
of your operating system. You can import bioscrape if you run Python
from the lib directory, IF you first copy the **init**.py file from the
main directory into the build/lib/bioscrape directory.

Testing and Examples
--------------------

To test the code and see some examples, run the one of the `example
notebooks <https://github.com/biocircuits/bioscrape/tree/master/examples>`__
and see if it works.

To run the testing notebook, from the main package directory, run the
command

::

   jupyter notebook

and you should have a browser window open with a file system. Click on
the file examples.ipynb to open the examples notebook. Hit Shift+Enter
to execute the cells that have code in them. Enjoy!

You can also run our unit test suite from the terminal by changing your
directory to the bioscrape root directory and running the command:

::

   pytest

If all the tests work fine, then you are successfully installed. Check
out the :doc:`usage page <examples>` to get going with actually using
the package.

.. |Bioscrape Core| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/examples/Basic%20Examples%20-%20START%20HERE.ipynb#scrollTo=Jmm8mTPfhMMS
.. |Bioscrape Sensitivity Analysis| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/examples/Sensitivity%20Analysis%20using%20Bioscrape.ipynb#scrollTo=Gu1r4H4ti_z7
.. |Bioscrape Inference| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/inference%20examples/Bioscrape%20Inference%20-%20Getting%20Started.ipynb#scrollTo=yvYVliBgjyzF
