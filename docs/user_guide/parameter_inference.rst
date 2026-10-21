Bioscrape Inference
===================

Bioscrape can be used to run Bayesian parameter inference. To get
started, make sure you have Python ``emcee`` package installed in your
environment. `emcee <https://emcee.readthedocs.io/en/stable/>`__ is
an implementation of a Markov chain Monte Carlo (MCMC) method that can
be used for parameter estimation.

Examples
---------

Ready-to-run examples are available in bioscrape to help you learn the
various tools available with bioscrape inference package:

1. Fitting noisy data to a linear Model (y = mx + b) - This is a replica
   of an `example
   <https://emcee.readthedocs.io/en/stable/tutorials/line/>`__
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

Tutorial style Jupyter notebooks are available in `this folder
<https://github.com/biocircuits/bioscrape/tree/master/\
inference%20examples>`__
for all of the above examples.

Inference Setup
--------------

Priors
~~~~~~

All prior types can be set to a inference object by using `set_prior`
function or in your `py_inference` function call in the ``prior``
keyword argument. The expected signature is as follows (a dictionary
with keys as parameter names and a list that gives the parameters for
different prior types)

::

   list_corresponding_to_prior_type = ['prior-type-key',
       prior_parameter1, prior_parameter2,...,'positive',
       custom_function]
   ## Here 'positive' is an optional keyword in the list that can make
   ## any prior accept only positive parameter values.
   ## The `custom_function` is an optional keyword, which has to be the
   ## last element in the list.
   ## This is used for a custom prior (see the last option in the list
   ## below)

   prior = {'param_name': list_corresponding_to_prior_type}
   pid.set_prior(prior)

-  Uniform: To give a uniform prior to a parameter use, prior-type-key:
   ``"uniform"``. The ``prior_parameter1`` is the lower bound of the
   uniform distribution and the ``prior_parameter2`` is the upper bound.

   Example: ``{'k1': ['uniform', 0,10]}`` sets the prior for the
   parameter ``k1`` as uniform(0,10).

-  Gaussian: To give a gaussian prior to a parameter use,
   prior-type-key: ``"gaussian"``. The ``prior_parameter1`` is the mean
   of the gaussian distribution and the ``prior_parameter2`` is the
   standard deviation.

   Example: ``{'k1': ['gaussian', 0, 10]}`` sets the prior for the
   parameter ``k1`` as gaussian(0,10).

   Example: ``{'k1': ['gaussian', 0, 10, 'positive']}`` sets the prior
   for the parameter ``k1`` as gaussian(0,10) but with only positive
   values.

-  Exponential: To give a exponential prior to a parameter use,
   prior-type-key: ``"exponential"``. The ``prior_parameter1`` is the
   inverse of the mean of exponential distribution.

   Example: ``{'k1': ['exponential', 5]}`` sets the prior for the
   parameter ``k1`` as exponential(5).

-  Gamma: To give a gamma prior to a parameter use, prior-type-key:
   ``"gamma"``. The ``prior_parameter1`` is the shape parameter (alpha)
   and the ``prior_parameter2`` is the rate parameter (beta) of the
   gamma distribution.

   Example: ``{'k1': ['gamma', 2, 3]}`` sets the prior for the parameter
   ``k1`` as gamma(alpha=2, beta=3).

-  Beta: To give a beta prior to a parameter use, prior-type-key:
   ``"beta"``. The ``prior_parameter1`` and ``prior_parameter2`` are the
   alpha and beta shape parameters of the beta distribution. The
   parameter value is expected to lie in ``[0, 1]``.

   Example: ``{'k1': ['beta', 2, 3]}`` sets the prior for the parameter
   ``k1`` as beta(alpha=2, beta=3).

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

-  Custom function: ``“custom”`` type with signature: (param_name
   (str), param_value(float)). Build your own custom prior function that
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

Joint custom priors
~~~~~~~~~~~~~~~~~~~

Use the optional ``custom_joint_prior`` argument to add a prior contribution
that depends on multiple inferred parameters. The callable receives a
dictionary containing every inferred parameter name and its proposed value,
and returns a log-prior contribution. A non-finite return value rejects the
complete parameter tuple. For example::

   def custom_joint_prior(params_dict):
       if params_dict["k1"] > params_dict["k2"]:
           return np.inf
       return 0.0

   sampler, pid = py_inference(
       Model=M,
       ...,
       prior=prior,
       custom_joint_prior=custom_joint_prior,
   )

This option is independent of the parameter-wise entries in ``prior``;
existing built-in and ``custom`` parameter priors are evaluated first.

Batched parallel inference
--------------------------

Use ``parallel_method='batched'`` with ``parallel=True`` to evaluate priors
in the parent process and distribute only the accepted parameter tuples::

    sampler, pid = py_inference(
        Model=model,
        params_to_estimate=parameter_names,
        exp_data=data,
        measurements=['GFP'],
        prior=prior,
        custom_joint_prior=custom_joint_prior,
        parallel=True,
        n_processes=4,
        parallel_method='batched',
    )

This method supports all existing built-in priors, parameter-wise custom
priors, and ``custom_joint_prior``. Prior callbacks keep their existing
signatures and return log-prior contributions. A non-finite log prior
rejects the tuple before it is sent to a worker. Finite contributions are
added to the returned log likelihood exactly once. A custom prior can,
for example, evaluate a trained surrogate model in the parent process.

Each worker constructs its model and deterministic or stochastic likelihood
once and reuses them across proposals. Tasks carry parameter vectors;
models and experimental data are sent during worker initialization.
Results are returned to emcee in proposal order. emcee's proposal and
acceptance rules are unchanged. The batch callback manages the likelihood
pool through emcee's ``vectorize=True`` interface; individual prior
callbacks do not need to accept arrays.

For an expensive joint prior that supports array evaluation, supply the
optional ``custom_joint_prior_batch`` callback::

    def joint_prior_batch(params):
        return np.where(params['k1'] <= params['k2'], 0.0, -np.inf)

    sampler, pid = py_inference(
        Model=model,
        ...,
        prior=prior,
        custom_joint_prior_batch=joint_prior_batch,
        parallel=True,
        n_processes=4,
        parallel_method='batched',
    )

The dictionary values are one-dimensional arrays of model parameter values,
with one entry per proposal. Return a one-dimensional array containing one
log-prior contribution per proposal in the same order. A non-finite entry
rejects that proposal. The parameter-wise priors run first; only tuples
passing them reach the batch joint callback. Empty and fully rejected
batches do not call it. The callback remains in the parent process and is
not sent to likelihood workers.

``custom_joint_prior`` and ``custom_joint_prior_batch`` are scalar and batch
implementations of the same contribution. If both are supplied, their
values must agree: scalar evaluation uses the scalar callback, whereas
batched evaluation uses the batch callback. Their contributions are never
added together. A batch-only callback also works with scalar inference,
using one-proposal batches. Existing scalar callbacks continue to work
without modification. Callback exceptions and invalid output shapes raise
errors. Parameter transformations and the inference target are unchanged.

``parallel_method='pool'`` is the default and uses the existing
per-proposal multiprocessing implementation. Both methods use
``n_processes`` workers. The batched method records evaluated proposal
values in the parent object's ``cost_progress`` and ``cost_params``.
As with serial inference, these histories grow with the number of evaluated
proposals. Exceptions from priors and likelihood workers propagate to the
caller.

On systems that spawn worker processes, including Windows, launch inference
from a script's ``if __name__ == '__main__':`` block. The verification
notebook ``inference_examples/individual examples/Batched parallel
inference.ipynb`` exercises numerical agreement, rejection, worker reuse,
and repeated-batch timing. Stochastic likelihoods retain their existing
simulation random-number behavior; changing task scheduling can change
their sampled trajectories.

Data Types
~~~~~~~~~~

Experimental data passed to inference is held in one of the
`~bioscrape.inference.Data` subclasses, each corresponding to a
different kind of measurement:

-  `~bioscrape.inference.BulkData`: a single output measured for a whole
   culture over time (e.g. total fluorescence or optical density from a
   plate reader), where the measured value is effectively summed or
   averaged over all the cells. This is a bulk time series.
-  `~bioscrape.inference.FlowData`: a *distribution* of single-cell
   outputs at one or more time points, with individual cells not tracked
   from one time point to the next (e.g. flow cytometry or mRNA FISH).
   This is a population snapshot.
-  `~bioscrape.inference.StochasticTrajectories`: one or more individual
   cells tracked over time (e.g. time-lapse microscopy). This is a
   single-cell time series.

When using `~bioscrape.inference.py_inference`, the ``sim_type``
argument selects the data type and likelihood together:
``sim_type="deterministic"`` (the default) uses
`~bioscrape.inference.BulkData`, and ``sim_type="stochastic"`` uses
`~bioscrape.inference.StochasticTrajectories`.
`~bioscrape.inference.FlowData` is used by constructing the likelihood
object directly rather than through `~bioscrape.inference.py_inference`.

Likelihood options
~~~~~~~~~~~~~~~~~~

Each type of experimental data is paired with a likelihood class that
computes the log-likelihood of a parameter set given that data:

-  `~bioscrape.inference.DeterministicLikelihood`: compares a single
   deterministic (ODE) simulation against bulk time-series data. Used
   when ``sim_type`` is ``"deterministic"`` (the default).
-  `~bioscrape.inference.StochasticTrajectoriesLikelihood`: runs several
   stochastic simulations (see the ``N_simulations`` option) and scores
   how close the simulated single-cell trajectories are to the measured
   ones. Used when ``sim_type`` is ``"stochastic"``.
-  `~bioscrape.inference.StochasticTrajectoryMomentLikelihood`: a
   variant of the previous likelihood that compares the *moments* of the
   trajectories -- their means with ``Moments=1``, or means and second
   moments with ``Moments=2`` -- rather than individual trajectories.
-  `~bioscrape.inference.StochasticStatesLikelihood`: compares the
   *distribution* of simulated states to a measured distribution at each
   time point, for population-snapshot data.

MCMC parameter set-up options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``nwalkers`` : The number of Markov chains to initialize, called the
   number of walkers. Use ``set_nwalkers`` for the inference object.
-  ``init_seed``: Initial seeding for each Markov chain around the
   initial guesses for parameters. This argument is a float that sets
   the initial parameter guess to start MCMC as follows:

::

   ## p0 is the initial parameter vector set in MCMC algorithm
   p0 = np.array(params_values) + init_seed * np.random.randn(
       self.nwalkers, ndim)
   ## params_values are the list of all parameter values and ndim =
   ## number of parameters to estimate

Use ``set_init_seed`` for the inference object.

-  ``nsteps`` : Length of each Markov chain, called the number of
   walkers. Use ``set_nsteps`` for the inference object.
