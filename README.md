# single_cell_package
Stochastic Simulation and Bayesian Parameter Inference for Single Cell Data

Installation: 

Prereqs:

To install, make sure you have a C++ compiler (eg g++) as well as updated versions of Cython, Numpy, and Scipy. The best way to get all this working is to just use the Anaconda Python distribution.

Building:

From the same directory that this README file is in (i.e. the base directory for the project), run the following command.

python setup.py build_ext --inplace

This will build and install the libraries in the current directory.

Testing:

To run, run the command "jupyter notebook" from the directory that this README file is in. This should open up a Python notebook GUI in a web browser. Then, click on interactive_testing.ipynb to open the suite of examples and tests. You can navigate through the notebook and hit Shift+Enter to run different cells.
