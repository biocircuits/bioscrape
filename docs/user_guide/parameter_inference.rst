Bioscrape Inference
===================

Bioscrape can be used to run Bayesian parameter inference. To get
started, make sure you have Python ``emcee`` package installed in your
environment. ```emcee`` <https://emcee.readthedocs.io/en/stable/>`__ is
an implementation of a Markov chain Monte Carlo (MCMC) method that can
be used for parameter estimation.

Examples:
---------

Ready-to-run examples are available in bioscrape to help you learn the
various tools available with bioscrape inference package:

1. Fitting noisy data to a linear Model (y = mx + b) - This is a replica
   of an
   `example <https://emcee.readthedocs.io/en/stable/tutorials/line/>`__
   demonstrated in the emcee package documentation - but with bioscrape.
2. Parameter estimation of a Birth-death process.
3. Using multiple trajectories, time points, and initial conditions in
   your experimental data is very easy with bioscrape inference. An
   example demonstrating this is available in the inference examples
   directory.
4. Stochastic inference (to run stochastic simulations for likelihood
   computation) options are available as well - simply set the
   ``sim_type`` to stochastic and that's it. Check out the related
   example to learn how to use it.

Tutorial style Jupyter notebooks are available in `this
folder <https://github.com/biocircuits/bioscrape/tree/master/inference%20examples>`__
for all of the above examples.

Documentation:
--------------

Priors
~~~~~~

All prior types can be set to a inference object by using ``set_prior``
function or in your ``py_inference`` function call in the ``prior``
keyword argument. The expected signature is as follows (a dictionary
with keys as parameter names and a list that gives the parameters for
different prior types)

::

   list_corresponding_to_prior_type = ['prior-type-key', prior_parameter1, prior_parameter2,...,'positive',custom_function]
   ## Here 'positive' is an optional keyword in the list that can make any prior accept only positive parameter values.
   ## The `custom_function` is an optional keyword, which has to be the last element in the list. 
   ## This is used for a custom prior (see the last option in the list below)

   prior = {'param_name': list_corresponding_to_prior_type}
   pid.set_prior(prior)

-  Uniform: To give a uniform prior to a parameter use, prior-type-key:
   ``"uniform"``. The ``prior_parameter1`` is the lower bound of the
   uniform distribution and the ``prior_parameter2`` is the upper bound.

Example: ``{'k1': ['uniform', 0,10]}`` sets the prior for the parameter
``k1`` as uniform(0,10).

-  Gaussian: To give a gaussian prior to a parameter use,
   prior-type-key: ``"gaussian"``. The ``prior_parameter1`` is the mean
   of the gaussian distribution and the ``prior_parameter2`` is the
   standard deviation.

Example: ``{'k1': ['gaussian', 0, 10]}`` sets the prior for the
parameter ``k1`` as gaussian(0,10).

Example: ``{'k1': ['gaussian', 0, 10, 'positive']}`` sets the prior for
the parameter ``k1`` as gaussian(0,10) but with only positive values.

-  Exponential: To give a exponential prior to a parameter use,
   prior-type-key: ``"exponential"``. The ``prior_parameter1`` is the
   inverse of the mean of exponential distribution.

Example: ``{'k1': ['exponential', 5]}`` sets the prior for the parameter
``k1`` as exponential(5).

-  Log Uniform: To give a log uniform prior to a parameter use,
   prior-type-key: ``"log-uniform"``. The ``prior_parameter1`` is the
   lower bound of the log uniform probability distribution and the
   ``prior_parameter2`` is the upper bound.

Example: ``{'k1': ['log-uniform', 5, 10]}`` sets the prior for the
parameter ``k1`` as log-uniform(5, 10).

-  Log Gaussian: To give a log gaussian prior to a parameter use,
   prior-type-key: ``"log-gaussian"``. The ``prior_parameter1`` is the
   mean of the gaussian distribution and the ``prior_parameter2`` is the
   standard deviation.

Example: ``{'k1': ['log-gaussian', 0, 10]}`` sets the prior for the
parameter ``k1`` as log-gaussian(0,10).

-  Custom function: ``“custom”`` type with signature: (param_name (str),
   param_value(float)). Build your own custom prior function that
   returns the log-probability given the parameter name and the
   parameter value as the function arguments.

Example: ``{'k1': ['custom', custom_prior]}`` sets the prior for the
parameter ``k1`` as a custom function which is the callable
``custom_prior``.

::

   def custom_prior(param_name, param_value):
       '''
       Returns log-probability
       '''
       # sanity check
       return np.Inf
       # calculate probability p
       return p

For an example on how to run custom prior distribution, check out the
prior distribution notebook in ``inference examples`` directory.

Likelihood options
~~~~~~~~~~~~~~~~~~

MCMC parameter set-up options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``nwalkers`` : The number of Markov chains to initialize, called the
   number of walkers. Use ``set_nwalkers`` for the inference object.
-  ``init_seed``: Initial seeding for each Markov chain around the
   initial guesses for parameters. This argument is a float that sets
   the initial parameter guess to start MCMC as follows:

::

   ## p0 is the initial parameter vector set in MCMC algorithm
   p0 = np.array(params_values) + init_seed * np.random.randn(self.nwalkers, ndim)
   ## params_values are the list of all parameter values and ndim = number of parameters to estimate

Use ``set_init_seed`` for the inference object.

-  ``nsteps`` : Length of each Markov chain, called the number of
   walkers. Use ``set_nsteps`` for the inference object.
