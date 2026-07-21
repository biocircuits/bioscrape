Overview
========

Bioscrape has three primary features related to modeling, analysis, simulations, and parameter inference:

1. Models can be :doc:`imported from SBML files <sbml_support>`
   and modified using the internal API. Other than changes to species, parameters, and reactions, users can also add :doc:`delays <delays>` to reactions and add assignment rules. Bioscrape provides sensitivity analysis tools for these models to determine how sensitive the model is to local changes in parameters.

2. Bioscrape includes :doc:`fast deterministic and stochastic simulators <simulators>` written
   in `Cython <https://cython.org/>`__ that can be used to simulate the
   model. These simulators are flexible supporting many different
   :doc:`chemical reaction models <propensities>` and are much faster
   than other packages where the simulation is done in pure Python.
   Additionally, Bioscrape allows for the
   simulation of growing and dividing cells via its :doc:`cell lineage
   simulators <lineage_package>`. However, deterministic delayed
   reactions are not yet supported when simulating cell lineages
   (stochastic delay reactions are). It is also possible to write your
   own simulator, propensities, delays, or cell lineage model that
   easily integrate with Bioscrape.

3. Bioscrape provides a easy-to-use :doc:`Bayesian parameter
   inference <parameter_inference>` module that uses the fast simulators
   underneath as well as existing Markov Chain Monte Carlo (MCMC) libraries to
   do parameter inference given data. Bayesian inference sampler using
   Python `emcee <https://emcee.readthedocs.io/en/stable/>`__ is
   built-in with Bioscrape that allows for estimation of parameter
   distributions from time-series data. Multiple data conditions,
   multiple initial conditions, time points, and many other common
   experimental situations are already addressed with helpful bioscrape
   inference functions. Additionally, a deterministic parameter inference interface
   using Python `LMFit <https://lmfit.github.io/lmfit-py/index.html>`__
   is also built-in. 

