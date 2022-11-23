# Bioscrape &mdash; Biological Stochastic Simulation of Single Cell Reactions and Parameter Estimation
## Python toolbox to simulate, analyze, and learn biological system models

[![Build Status](https://github.com/biocircuits/bioscrape/actions/workflows/bioscrape.yml/badge.svg)](https://github.com/biocircuits/bioscrape/actions/workflows/bioscrape.yml)
[![PyPI version](https://badge.fury.io/py/bioscrape.svg)](https://badge.fury.io/py/bioscrape)

* Getting started with Bioscrape: [![Bioscrape Core](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/examples/Basic%20Examples%20-%20START%20HERE.ipynb#scrollTo=Jmm8mTPfhMMS)

* Bioscrape analysis features: [![Bioscrape Sensitivity Analysis](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/examples/Sensitivity%20Analysis%20using%20Bioscrape.ipynb#scrollTo=Gu1r4H4ti_z7)

* Parameter inference with Bioscrape: [![Bioscrape Inference](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/inference%20examples/Bioscrape%20Inference%20-%20Getting%20Started.ipynb#scrollTo=yvYVliBgjyzF)

Bioscrape is a Systems Biology Markup Language (SBML) simulator written in Cython for speed and python compatibility. It can be used for deterministic, stochastic, or single cell simulation and also has parameter inference capabilities.

- **Mailing list:** [SBTools Google Group](https://groups.google.com/g/sbtools/) Email: sbtools@googlegroups.com
- **Source:** https://github.com/biocircuits/bioscrape
- **Preprint:** - [Fast and flexible simulation and parameter estimation for synthetic biology using bioscrape](https://www.biorxiv.org/content/10.1101/121152v3)
- **Bug reports:** https://github.com/biocircuits/bioscrape/issues
- **Slack** Join the #bioscrape channel on SBTools slack: Ask on the public SBTools Google group to be added or send a message to one of the maintainers. 

# Example 1: Simulating an SBML

Bioscrape allows for deterministic and stochastic simulation of SBML models:

```python
from bioscrape.types import Model
# Load an SBML file repressilator.xml 
# (you can find this file in `examples/models` directory)
M = Model(sbml_filename = 'repressilator_sbml.xml')
# Simulate the model
from bioscrape.simulator import py_simulate_model
import numpy as np
tp = np.linspace(0,256,100)
result = py_simulate_model(timepoints=tp, Model=M, stochastic=True)
# Plot the simulation result (the result is a Pandas dataframe)
import matplotlib.pyplot as plt
plt.plot(tp, result['X'])
```

# Example 2: Run Bayesian inference with Bioscrape 

Bioscrape can be used to identify model parameters using experimental data. In the example below, we show the user-friendly plug-and-play nature of bioscrape inference. We load the data as a Pandas dataframe and the model as an SBML file. The Bayesian inference is implemented as a wrapper for Python emcee that implements Markov Chain Monte Carlo (MCMC) sampler. Bioscrape inference provides various features such as: multiple data conditions, multiple data trajectories, deterministic inference, automatic visualization of posteriors, convergence checking tools, built-in and customizable priors, and lots more!

```python
from bioscrape.types import Model
import pandas as pd
from bioscrape.inference import py_inference

# Load an SBML model 
# (you can get this file in `inference examples/models/` directory)
M = Model(sbml_filename='toy_sbml_model.xml')

# Load experimental data 
# (you can find test data in `inference examples/data/` directory)
df = pd.read_csv('test_data.csv', delimiter = '\t', 
                 names = ['X','time'], skiprows = 1)

# Use built-in priors, 
# For 'd1': a Gaussian distribution of mean 0.2 and standard deviation of 20,
# while ensuring the parameter remains positive
# For 'k1': a Uniform distribution with minimum value 0 and maximum value 100

prior = {'d1' : ['gaussian', 0.2, 20, 'positive'], 'k1' : ['uniform', 0, 100]}

# Run Bayesian inference
sampler, pid = py_inference(Model = M, exp_data = df, measurements = ['X'], 
                            time_column = ['time'], nwalkers = 20, nsteps = 5500,
                            params_to_estimate = ['d1', 'k1'], prior = prior)
# A sampler object containing all samples is returned.
# The pid object consists of various utilities for further analysis.
# This will plot the resulting posterior parameter distributions as well.
```


All examples can be found in the [examples](https://github.com/biocircuits/bioscrape/tree/master/examples), the [inference examples](https://github.com/biocircuits/bioscrape/tree/master/inference%20examples), and the [lineage examples](https://github.com/biocircuits/bioscrape/tree/master/lineage%20examples) folders. If you prefer to run the package without installing the package, please use the Google Colab links above. If you want a local installation for bioscrape (recommended for faster speeds), follow the steps below: 

# Installation

Install the latest version of Bioscrape::

    $ pip install bioscrape
    

Please note that Bioscrape is a Cython extension module and requires a C++ compiler to be set up on your computer for installation.

Try online without installing, open self-explanatory jupyter notebooks with Google Colab (linked at the top of this README).

Further details about the installation process can be found in the [Bioscrape wiki](https://github.com/biocircuits/bioscrape/wiki#installation).

# Bugs and Contributing to Bioscrape

Please report any bugs that you find [here](https://github.com/biocircuits/bioscrape/issues).
Or, even better, fork the repository on [GitHub](https://github.com/biocircuits/bioscrape),
and create a pull request (PR). We welcome all changes, big or small, and we
will help you make the PR if you are new to `git` (just ask on the issue).

# Versions

Bioscrape versions:

* 1.1.1 (latest release): To install run `pip install bioscrape==1.1.1` 
* 1.1.0 (latest stable release): To install run `pip install bioscrape`
* 1.0.4 (beta release): To install run `pip install bioscrape==1.0.4`

# License
Released under the MIT License (see `LICENSE`)

Copyright (c) 2022, Biocircuits, California Institute of Technology. All rights reserved.
