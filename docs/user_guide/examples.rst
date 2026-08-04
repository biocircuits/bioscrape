Examples
========

Example 1: Simulating an SBML
-----------------------------

Bioscrape allows for deterministic and stochastic simulation of SBML
models:

.. code:: python

   from bioscrape.types import Model
   ## Load an SBML file repressilator.xml 
   ## (you can find this file in `examples/models` directory)
   M = Model(sbml_filename = 'repressilator_sbml.xml')
   ## Simulate the model
   from bioscrape.simulator import py_simulate_model
   import numpy as np
   tp = np.linspace(0,256,100)
   result = py_simulate_model(timepoints=tp, Model=M, stochastic=True)
   ## Plot the simulation result (the result is a Pandas dataframe)
   import matplotlib.pyplot as plt
   plt.plot(tp, result['X'])

Example 2: Run Bayesian inference with Bioscrape
------------------------------------------------

Bioscrape can be used to identify model parameters using experimental
data. In the example below, we show the user-friendly plug-and-play
nature of bioscrape inference. We load the data as a Pandas dataframe
and the model as an SBML file. The Bayesian inference is implemented as
a wrapper for Python emcee that implements Markov Chain Monte Carlo
(MCMC) sampler. Bioscrape inference provides various features such as:
multiple data conditions, multiple data trajectories, deterministic
inference, automatic visualization of posteriors, convergence checking
tools, built-in and customizable priors, and lots more!

.. code:: python

   from bioscrape.types import Model
   import pandas as pd
   from bioscrape.inference import py_inference

   ## Load an SBML model 
   ## (you can get this file in `inference examples/models/` directory)
   M = Model(sbml_filename='toy_sbml_model.xml')

   ## Load experimental data 
   ## (you can find test data in `inference examples/data/` directory)
   df = pd.read_csv('test_data.csv', delimiter = '\t', 
                    names = ['X','time'], skiprows = 1)

   ## Use built-in priors, 
   ## For 'd1': a Gaussian distribution of mean 0.2 and standard deviation of 20,
   ## while ensuring the parameter remains positive
   ## For 'k1': a Uniform distribution with minimum value 0 and maximum value 100

   prior = {'d1' : ['gaussian', 0.2, 20, 'positive'], 'k1' : ['uniform', 0, 100]}

   ## Run Bayesian inference
   sampler, pid = py_inference(Model = M, exp_data = df, measurements = ['X'], 
                               time_column = ['time'], nwalkers = 20, nsteps = 5500,
                               params_to_estimate = ['d1', 'k1'], prior = prior)
   ## A sampler object containing all samples is returned.
   ## The pid object consists of various utilities for further analysis.
   ## This will plot the resulting posterior parameter distributions as well.

All examples can be found in the
`examples <https://github.com/biocircuits/bioscrape/tree/master/examples>`__,
the `inference
examples <https://github.com/biocircuits/bioscrape/tree/master/inference%20examples>`__,
and the `lineage
examples <https://github.com/biocircuits/bioscrape/tree/master/lineage%20examples>`__
folders.
