Overview
========

There are fundamentally 3 parts to the package. If you want to simply
*simulate* some biological circuits and maybe do parameter sweeps, you
only need to know about parts 1 and 2. If you want to do parameter
inference as well, then you'll need all 3 parts. If you have your own
preferred method for parameter inference, you can still use parts 1 and
2 and build your own inference method on top of them.

1. The first is a :doc:`reproducible and standardized <sbml_support>`
   support for modeling species, reactions, and reaction rates in a
   easily modifiable way using an API. This way, if you specify a model
   and want to tweak it later to get different behavior or better fits,
   you can just quickly edit the model using the API and proceed rather
   than re-write a ton of code/edit SBML directly.

2. The second is a set of :doc:`fast simulators <simulators>` written
   in `Cython <https://cython.org/>`__ that can be used to simulate the
   model. These simulators are flexible supporting many different
   :doc:`chemical reaction models <propensities>` and are much faster
   than other packages where the simulation is done in pure Python.
   There are deterministic and stochastic simulators which also support
   :doc:`delays <delays>`. Additionally, Bioscrape allows for the
   simulation of growing and dividing cells via its :doc:`cell lineage
   simulators <lineage_package>`. However, deterministic delayed
   reactions are not yet supported when simulating cell lineages
   (stochastic delay reactions are). It is also possible to write your
   own simulator, propensities, delays, or cell lineage model that
   easily integrate with bioSCRAPE.

3. The third piece is a set of methods for :doc:`Bayesian parameter
   inference <parameter_inference>` that use the fast simulators
   underneath as well as existing Markov chain Monte Carlo libraries to
   do parameter inference given data. In this scenario, you specify a
   model, you specify the parameters you would like to identify, and you
   specify the data that you have measured. The inference method then
   runs for a while (could be a loooong time) and spits out parameter
   estimates as well as distributions. Bayesian inference sampler using
   Python `emcee <https://emcee.readthedocs.io/en/stable/>`__ is
   built-in with Bioscrape that allows for estimation of parameter
   distributions from time-series data. Multiple data conditions,
   multiple initial conditions, time points, and many other common
   experimental situations are already addressed with helpful bioscrape
   inference functions. A deterministic parameter inference interface
   using Python `LMFit <https://lmfit.github.io/lmfit-py/index.html>`__
   is also built-in with Bioscrape inference. The deterministic
   inference solves a nonlinear optimization problem to find best
   estimates of parameters given a data trajectory. Other than built-in
   inference, Bioscrape allows for easy implementation of new parameter
   inference that use different cost functions or apply parameter
   identification from different kinds of data (such as distributional
   data, steady-state data etc.) Refer to the inference documentation
   for more information.

In addition to all of these, Bioscrape lends itself to easy extension of
analysis tools. For example, users may run sensitivity analysis on
Bioscrape models.
