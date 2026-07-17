"""
Orchestrates Bayesian and least-squares parameter inference for
bioscrape models.

`InferenceSetup` wraps `emcee <https://emcee.readthedocs.io/>`_
(Bayesian/MCMC) and `lmfit <https://lmfit.github.io/lmfit-py/>`_
(least-squares) inference around a `~bioscrape.types.Model`, using the
`~bioscrape.pid_interfaces.PIDInterface` subclasses to build the
likelihood being sampled or minimized. Most users should use the
higher-level `~bioscrape.inference.py_inference` entry point, which
constructs an `InferenceSetup` and returns it (along with the `emcee`
sampler) rather than requiring one to be built directly.
"""
import numpy as np
import warnings
try:
    import pandas as pd
except:
    warnings.warn('Pandas package not found.')
import sys
import matplotlib.pyplot as plt
from bioscrape.types import Model
from bioscrape.sbmlutil import import_sbml as sbmlutil_import_sbml
from bioscrape.simulator import ModelCSimInterface, DeterministicSimulator, SSASimulator
from bioscrape.pid_interfaces import *

def initialize_inference(**kwargs):
    """
    Construct an `InferenceSetup`.

    Parameters
    ----------
    **kwargs
        Forwarded to `InferenceSetup`.

    Returns
    -------
    InferenceSetup
        The constructed inference object.
    """
    return InferenceSetup(**kwargs)

class InferenceSetup(object):
    """
    Sets up and runs parameter inference for a bioscrape model.

    Constructed (and returned) by `~bioscrape.inference.py_inference`;
    most users will not need to construct this directly. Holds the
    model, data, prior, and inference settings, and provides `run_mcmc`
    (Bayesian inference via `emcee`) and `run_lmfit` (least-squares fit
    via `lmfit`) to actually perform the inference.

    Parameters
    ----------
    Model : Model, optional
        The bioscrape model to fit.
    params_to_estimate : list of str, optional
        The names of the parameters in `Model` to estimate.
    exp_data : pandas.DataFrame or list of pandas.DataFrame, optional
        The experimental data to fit against.
    prior : dict, optional
        The prior for each parameter in `params_to_estimate`; see
        `set_prior`.
    measurements : list of str, optional
        The name(s) of the species in `Model` that are measured
        (default `['']`).
    time_column : str, default 'time'
        The column name in `exp_data` that holds the time points.
    timepoints : list or numpy.ndarray, optional
        Force the timepoints used for inference, instead of extracting
        them from `exp_data`.
    initial_conditions : dict or list of dict, optional
        Initial condition(s) for the simulated trajectories, in the
        same shape as `exp_data`. Defaults to `Model`'s own initial
        conditions.
    parameter_conditions : dict or list of dict, optional
        Per-trajectory parameter overrides, in the same shape as
        `exp_data`.
    sim_type : {'deterministic', 'stochastic'}, default 'deterministic'
        The type of simulation to use for the inference run.
    method : {'emcee', 'lmfit'}, default 'emcee'
        The inference package to use; see `set_method`.
    nwalkers : int, default 100
        The number of walkers for MCMC sampling; see
        `emcee.EnsembleSampler`.
    nsteps : int, default 1000
        The number of steps for MCMC sampling; see
        `emcee.EnsembleSampler`.
    init_seed : float or numpy.ndarray or list or 'prior', default 0.01
        Controls the initial parameter values for the MCMC sampler. A
        float `r` samples initial values from a Gaussian ball of
        radius `r` around the model's parameter values; an array or
        list of length `len(params_to_estimate)` is used directly as
        the initial values; 'prior' samples the initial values from
        `prior`.
    norm_order : int, default 2
        The norm used to compute the log-likelihood distance between
        simulated and measured trajectories; see
        `~bioscrape.inference.StochasticTrajectoriesLikelihood`.
    N_simulations : int, default 3
        The number of stochastic simulations to average per data point
        when `sim_type='stochastic'`.
    parallel : bool, default False
        If True, run MCMC sampling with a `multiprocessing.Pool`.
    debug : bool, default False
        If True, print verbose diagnostic messages while preparing and
        running inference.
    dimension : int, default 1
        Stored on the object but not currently used elsewhere.
    hmax : float, optional
        Stored on the object but not currently used elsewhere.

    Attributes
    ----------
    LL_data : numpy.ndarray
        The experimental data reshaped to `(N, T, M)` (trajectories by
        timepoints by measured species) by `extract_data`; set when
        `exp_data` is given at construction, or by `prepare_inference`.
    cost_progress : list
        The log-likelihood value computed at each `cost_function` call
        during a `run_mcmc`/`run_emcee` run.
    cost_params : list
        The parameter values passed to `cost_function` at each call
        during a `run_mcmc`/`run_emcee` run, in the same order as
        `cost_progress`.
    """
    def __init__(self, **kwargs):
        """See class docstring."""
        self.M = None
        self.M = kwargs.get('Model', None)
        self.params_to_estimate = kwargs.get('params_to_estimate', [])
        self.init_seed = kwargs.get('init_seed', 0.01)
        self.prior = kwargs.get('prior', None)
        self.nwalkers = kwargs.get('nwalkers', 100)
        self.nsteps = kwargs.get('nsteps', 1000)
        self.dimension = kwargs.get('dimension', 1)
        self.exp_data = kwargs.get('exp_data', None)
        self.sim_type = kwargs.get('sim_type', 'deterministic')
        self.method = kwargs.get('method', 'emcee')
        self.timepoints = kwargs.get('timepoints', None)
        self.time_column = kwargs.get('time_column', 'time')
        self.measurements = kwargs.get('measurements', [''])
        self.initial_conditions = kwargs.get('initial_conditions', None)
        self.parameter_conditions = kwargs.get('parameter_conditions', None)
        self.norm_order = kwargs.get('norm_order', 2)
        self.N_simulations = kwargs.get('N_simulations', 3)
        self.LL_data = None
        self.debug = kwargs.get('debug', False)
        self.cost_progress = []
        self.cost_params = []
        self.hmax = kwargs.get('hmax', None)
        self.parallel = kwargs.get('parallel', False)
        print("Received parallel as", self.parallel)
        if self.exp_data is not None:
            self.prepare_inference()
            self.setup_cost_function()
        return

    #Whenever new settings are updated in the constructor, please add them to getstate and setstate
    def __getstate__(self):
        """Return this object's state as a tuple of picklable objects."""
        return (
            self.M,
            self.params_to_estimate,
            self.init_seed,
            self.prior,
            self.nwalkers,
            self.nsteps,
            self.dimension,
            self.exp_data,
            self.sim_type,
            self.method,
            self.timepoints,
            self.time_column,
            self.measurements,
            self.initial_conditions,
            self.parameter_conditions,
            self.norm_order,
            self.N_simulations,
            self.LL_data,
            self.debug,
            self.cost_progress,
            self.cost_params,
            self.hmax,
            self.parallel
            )

    def __setstate__(self, state):
        """Restore this object's state from a tuple produced by `__getstate__`."""
        self.M = state[0]
        self.params_to_estimate = state[1]
        self.init_seed = state[2]
        self.prior = state[3]
        self.nwalkers = state[4]
        self.nsteps = state[5]
        self.dimension = state[6]
        self.exp_data = state[7]
        self.sim_type = state[8]
        self.method = state[9]
        self.timepoints = state[10]
        self.time_column = state[11]
        self.measurements = state[12]
        self.initial_conditions = state[13]
        self.parameter_conditions = state[14]
        self.norm_order = state[15]
        self.N_simulations = state[16]
        self.LL_data = state[17]
        self.debug = state[18]
        self.cost_progress = state[19]
        self.cost_params = state[20]
        self.hmax = state[21]
        self.parallel = state[22]
        if self.exp_data is not None:
            self.prepare_inference()
            self.setup_cost_function()


    def set_model(self, M):
        """
        Set the bioscrape model to fit.

        Parameters
        ----------
        M : Model
            The bioscrape model to fit.

        Returns
        -------
        bool
            True.
        """
        self.M = M
        return True

    def get_model(self):
        """
        The bioscrape model being fit.

        Returns
        -------
        Model
            The bioscrape model to fit.
        """
        return self.M

    def set_prior(self, prior):
        """
        Set the prior distribution for each parameter to estimate.

        Parameters
        ----------
        prior : dict
            Maps each parameter name in `self.params_to_estimate` to a
            list describing its prior distribution; see
            :doc:`../inference` for the supported prior types.

        Returns
        -------
        bool
            True.

        Raises
        ------
        ValueError
            If `prior` does not have exactly one entry per parameter
            in `self.params_to_estimate`.
        """
        self.prior = prior
        if len(list(self.prior.keys())) != len(self.params_to_estimate):
            raise ValueError('Prior keys length must be equal to the length of params_to_estimate.')
        return True

    def set_nsteps(self, nsteps: int):
        """
        Set the number of MCMC steps to run.

        Parameters
        ----------
        nsteps : int
            The number of steps for MCMC sampling; see
            `emcee.EnsembleSampler`.

        Returns
        -------
        bool
            True.
        """
        self.nsteps = nsteps
        return True

    def set_nwalkers(self, nwalkers: int):
        """
        Set the number of MCMC walkers to use.

        Parameters
        ----------
        nwalkers : int
            The number of walkers for MCMC sampling; see
            `emcee.EnsembleSampler`.

        Returns
        -------
        bool
            True.
        """
        self.nwalkers = nwalkers
        return True

    def set_dimension(self, dimension: int):
        """
        Set the `dimension` attribute.

        Parameters
        ----------
        dimension : int
            Stored on the object but not currently used elsewhere.

        Returns
        -------
        bool
            True.
        """
        self.dimension = dimension
        return True

    def set_init_seed(self, init_seed: float):
        """
        Set the initial parameter values for the MCMC sampler.

        Parameters
        ----------
        init_seed : float or numpy.ndarray or list or 'prior'
            Controls the initial parameter values for the sampler; see
            the class docstring.

        Returns
        -------
        bool
            True.
        """
        self.init_seed = init_seed
        return True

    def set_timepoints(self, timepoints):
        """
        Force the timepoints to use for inference.

        By default, the timepoints are extracted from `self.exp_data`.

        Parameters
        ----------
        timepoints : list or numpy.ndarray
            The timepoints to use.

        Returns
        -------
        bool
            True.
        """
        self.timepoints = timepoints
        return True

    def set_params_to_estimate(self, params_to_estimate: list):
        """
        Set the list of parameters to estimate.

        Parameters
        ----------
        params_to_estimate : list of str
            The names of the parameters in the model (attached via
            `set_model`) to estimate.

        Returns
        -------
        bool
            True.
        """
        self.params_to_estimate = params_to_estimate
        return True

    def set_sim_type(self, sim_type: str):
        """
        Set the type of simulation to use for inference.

        Parameters
        ----------
        sim_type : {'deterministic', 'stochastic'}
            The type of simulation to use for the inference run.

        Returns
        -------
        bool
            True.
        """
        self.sim_type = sim_type
        return True

    def set_method(self, method: str):
        """
        Set the parameter identification method to use.

        Parameters
        ----------
        method : {'emcee', 'lmfit'}
            The inference package to use: 'emcee' for Bayesian (MCMC)
            inference (see `emcee <https://emcee.readthedocs.io/en/stable/>`_),
            or 'lmfit' for non-linear least-squares (see
            `lmfit <https://lmfit.github.io/lmfit-py/>`_).

        Returns
        -------
        bool
            True.
        """
        self.method = method
        return True

    def set_N_simulations(self, N_simulations: int):
        """
        Set the number of stochastic simulations averaged per data point.

        Only applies when `self.sim_type` is 'stochastic'.

        Parameters
        ----------
        N_simulations : int
            The number of simulations to run for each condition.

        Returns
        -------
        bool or None
            True if set; otherwise None, with a warning, if
            `self.sim_type` is not 'stochastic'.
        """
        if self.sim_type == 'stochastic':
            self.N_simulations = N_simulations
            return True
        else:
            warnings.warn('N_simulations needs to be set only for stochastic inference which is not currently the case.')

    def set_initial_conditions(self, initial_conditions):
        """
        Set the initial condition(s) for the simulated trajectories.

        Parameters
        ----------
        initial_conditions : dict or list of dict
            A dictionary mapping species names to initial values, or a
            list of such dictionaries (one per trajectory in
            `self.exp_data`, in the same order). If None, defaults to
            the model's own initial conditions (see `set_model`).

        Returns
        -------
        bool
            True.

        Raises
        ------
        ValueError
            If `initial_conditions` is a non-empty list containing a
            non-dict entry.
        """
        # Get initial_conditions from the model if not given explicitly
        if initial_conditions is None:
            initial_conditions = self.M.get_species_dictionary()
        if type(initial_conditions) is list and len(initial_conditions):
            for curr_ic_index in range(len(initial_conditions)):
                ic = initial_conditions[curr_ic_index]
                if type(ic) is not dict:
                    raise ValueError('All entries in the initial condition list must be dictionaries.')
        self.initial_conditions = initial_conditions
        return True

    def set_parameter_conditions(self, parameter_conditions):
        """
        Set per-trajectory parameter overrides.

        Parameters
        ----------
        parameter_conditions : dict or list of dict
            A dictionary mapping parameter names (that change between
            measurements) to values, or a list of such dictionaries
            (one per trajectory in `self.exp_data`, in the same
            order).

        Returns
        -------
        bool
            True.
        """
        self.parameter_conditions = parameter_conditions
        return True


    def set_measurements(self, measurements: list):
        """
        Set the list of measured species to look for in `self.exp_data`.

        Parameters
        ----------
        measurements : list of str
            The name(s) of the species in the model (attached via
            `set_model`) that are measured.

        Returns
        -------
        bool
            True.
        """
        self.measurements = measurements
        return True

    def set_time_column(self, time_column: str):
        """
        Set the name of the time column in `self.exp_data`.

        Parameters
        ----------
        time_column : str
            The column name in `self.exp_data` that holds the time
            points.

        Returns
        -------
        bool
            True.
        """
        self.time_column = time_column
        return True

    def set_exp_data(self, exp_data):
        """
        Set the experimental data to fit against.

        Parameters
        ----------
        exp_data : pandas.DataFrame or list of pandas.DataFrame
            The experimental data.

        Returns
        -------
        bool
            True.

        Raises
        ------
        ValueError
            If `exp_data` is not a `pandas.DataFrame` or a list.
        """
        if isinstance(exp_data, (pd.DataFrame, list)):
            self.exp_data = exp_data
        else:
            raise ValueError('exp_data must be either a Pandas dataframe or a list of dataframes.')
        return True

    def set_norm_order(self, norm_order: int):
        """
        Set the norm order used to compute the log-likelihood.

        Parameters
        ----------
        norm_order : int
            The norm used to compute the log-likelihood distance
            between simulated and measured trajectories.

        Returns
        -------
        bool
            True.
        """
        self.norm_order = norm_order
        return True

    def set_parallel(self, parallel: bool):
        """
        Set whether to use parallel processing for MCMC sampling.

        Parameters
        ----------
        parallel : bool
            If True, run MCMC sampling with a `multiprocessing.Pool`.

        Returns
        -------
        bool
            True.
        """
        self.parallel = parallel
        return True

    def get_parameters(self):
        """
        The list of parameters to estimate.

        Returns
        -------
        list of str
            The names of the parameters in the model (attached via
            `set_model`) to estimate.
        """
        return self.params_to_estimate

    def run_mcmc(self, **kwargs):
        """
        Run Bayesian (MCMC) parameter inference via `emcee`.

        Prepares the inference (extracting data, setting up the
        likelihood) and then runs `emcee.EnsembleSampler` for
        `self.nsteps` steps with `self.nwalkers` walkers.

        Parameters
        ----------
        timepoints : list or numpy.ndarray, optional
            Force the timepoints used for inference; see
            `set_timepoints`.
        norm_order : int, optional
            The norm used to compute the log-likelihood; see
            `set_norm_order`.
        N_simulations : int, optional
            The number of stochastic simulations to average per data
            point; see `set_N_simulations`.
        init_seed : float or numpy.ndarray or list or 'prior', optional
            The initial parameter values for the sampler; see
            `set_init_seed`.
        debug : bool, optional
            If True, print verbose diagnostic messages.
        reuse_likelihood : bool, default False
            If True, reuse the likelihood set up by a previous call
            instead of rebuilding it.
        progress : bool, default True
            Passed to `emcee.EnsembleSampler.run_mcmc` to show (or
            hide) a progress bar.
        skip_initial_state_check : bool, default False
            Passed to `emcee.EnsembleSampler.run_mcmc`.
        convergence_check : bool, default False
            If True, compute the autocorrelation time for each
            parameter after sampling (stored as
            `autocorrelation_time`).
        convergence_diagnostics : bool, default `convergence_check`
            If True, also store a summary of the autocorrelation time
            and acceptance fraction (as `convergence_diagnostics`).
        filename_csv : str, default 'mcmc_results.csv'
            Where to write the flattened MCMC chain, in the current
            directory unless a path is given.
        filename_txt : str, default 'mcmc_results.txt'
            Where to write a text summary of the cost function
            progress (and convergence diagnostics, if computed).
        printout : bool, default True
            If True, print status messages while running.

        Returns
        -------
        emcee.EnsembleSampler
            The sampler, after running `self.nsteps` steps.
        """
        self.prepare_inference(**kwargs)
        sampler = self.run_emcee(**kwargs)
        return sampler

    def prepare_inference(self, **kwargs):
        """
        Prepare to run inference: extract data and set up conditions.

        Applies any of the `timepoints`/`norm_order`/`N_simulations`/
        `debug` keyword arguments given, then prepares the initial and
        parameter conditions and extracts `self.exp_data` into
        `self.LL_data`.

        Parameters
        ----------
        timepoints : list or numpy.ndarray, optional
            Force the timepoints used for inference; see
            `set_timepoints`.
        norm_order : int, optional
            The norm used to compute the log-likelihood; see
            `set_norm_order`.
        N_simulations : int, optional
            The number of stochastic simulations to average per data
            point; see `set_N_simulations`.
        debug : bool, optional
            If True, print verbose diagnostic messages.

        Raises
        ------
        ValueError
            If `timepoints` is given and is not a list or
            `numpy.ndarray`.
        """
        timepoints = kwargs.get('timepoints')
        norm_order = kwargs.get('norm_order')
        N_simulations = kwargs.get('N_simulations')
        debug = kwargs.get('debug')
        if N_simulations:
            self.set_N_simulations(N_simulations)
        if norm_order:
            # (integer) Which norm to use: 1-Norm, 2-norm, etc.
            self.set_norm_order(norm_order)
        if debug:
            self.debug = debug
        if timepoints is not None:
            if isinstance(timepoints, (list, np.ndarray)):
                self.set_timepoints(timepoints)
            else:
                raise ValueError('Expected type list or np.ndarray for timepoints.')
        self.prepare_initial_conditions()
        self.prepare_parameter_conditions()
        self.LL_data = self.extract_data()
        return

    def prepare_initial_conditions(self):
        """
        Expand `self.initial_conditions` to one dict per trajectory.

        Raises
        ------
        ValueError
            If `self.initial_conditions` is a list whose length does not
            match the number of trajectories in `self.exp_data`.
        """
        # Create initial conditions as required
        N = 1 if type(self.exp_data) is dict else len(self.exp_data)
        if type(self.initial_conditions) is dict:
            all_initial_conditions = [self.initial_conditions]*N
        elif type(self.initial_conditions) is list:
            if len(self.initial_conditions) != N:
                raise ValueError('For a list of initial conditions,'
                                    'each item must be a dictionary and'
                                    'the length of the list must be the'
                                    'same as the number of trajectories.')
            all_initial_conditions = self.initial_conditions
        self.initial_conditions = all_initial_conditions
        return

    def prepare_parameter_conditions(self):
        """
        Expand `self.parameter_conditions` to one dict per trajectory.

        Raises
        ------
        ValueError
            If `self.parameter_conditions` is a list whose length does
            not match the number of trajectories in `self.exp_data`.
        AssertionError
            If a parameter in `self.params_to_estimate` also appears
            in a per-trajectory parameter condition.
        """
        # Create parameter conditions as required
        N = 1 if type(self.exp_data) is dict else len(self.exp_data)
        if type(self.parameter_conditions) is dict:
            all_parameter_conditions = [self.parameter_conditions]*N
        elif type(self.parameter_conditions) is list:
            if len(self.parameter_conditions) != N:
                raise ValueError('For a list of parameter conditions,'
                                    'each item must be a dictionary and'
                                    'the length of the list must be the'
                                    'same as the number of trajectories.')
            all_parameter_conditions = self.parameter_conditions
        else:
            all_parameter_conditions = None
        self.parameter_conditions = all_parameter_conditions
        # Make sure that parameters to estimate do not intersect with parameters
        if self.parameter_conditions is not None:
            # that are changing through parameter conditions
            for param_condition in self.parameter_conditions:
                for param in self.params_to_estimate:
                    assert param not in param_condition.keys()
        return

    def extract_data(self):
        """
        Reshape `self.exp_data` into a single `(N, T, M)` array.

        Extracts the timepoints (from `self.time_column`) and measured
        species (from `self.measurements`) from `self.exp_data`, which
        may be a single `pandas.DataFrame` or a list of them (one per
        trajectory), into one array of measurements indexed by
        trajectory, timepoint, and measured species.

        Returns
        -------
        numpy.ndarray
            Array of shape `(N, T, M)`: `N` trajectories, `T`
            timepoints, `M` measured species.

        Raises
        ------
        TypeError
            If `self.time_column` is not set, if `self.exp_data` (or
            one of its entries, when a list) is not a
            `pandas.DataFrame`, or if `self.exp_data` is of an
            unrecognized type.
        """
        exp_data = self.exp_data
        # Get timepoints from given experimental data
        if isinstance(self.timepoints, (list, np.ndarray)) and self.debug:
            warnings.warn('Timepoints given by user, not using the data to extract the timepoints automatically.')
        M = len(self.measurements)# Number of measurements
        if type(exp_data) is list:
            if len(exp_data) == 1:
                exp_data = exp_data[0]
        # Multiple trajectories case
        if type(exp_data) is list:
            N = len(exp_data)# Number of trajectories
            data_list_final = []
            timepoints_list = []
            for df in exp_data:
                data_list = []
                if type(df) is not pd.DataFrame:
                    raise TypeError('All elements of exp_data attribute of an InferenceSetup object must be Pandas DataFrame objects.')
                # Extract timepoints
                if self.time_column:
                    timepoint_i = np.array(df.get(self.time_column), dtype='double').flatten()
                    timepoints_list.append(timepoint_i)
                else:
                    raise TypeError('time_column attribute of InferenceSetup object must be a string.')
                # Extract measurements
                if type(self.measurements) is list and len(self.measurements) == 1:
                    data_list.append(np.array(df.get(self.measurements[0]), dtype='double'))
                elif type(self.measurements) is list and len(self.measurements) > 1:
                    for m in self.measurements:
                        # Error in multiple measurements
                        data_list.append(np.array(df.get(m), dtype='double'))
                # Number of timepoints
                T = len(timepoints_list[0])
                if T != len(timepoint_i):
                    warnings.warn('The length of timepoints for all experimental trajectories must be the same, they can have different timepoints but not length of timepoints.')
                data_i = np.array(data_list)
                data_i = np.reshape(data_i, (T, M))
                data_list_final.append(data_i)
            data = np.array(data_list_final)
            self.timepoints = timepoints_list
            T = len(timepoints_list[0])
            data = np.reshape(data, (N,T,M))
            if self.debug:
                print('N (Number of trajectories) = {0}'.format(N))
                print('T (Length of timepoints) = {0}'.format(T))
                print('M (Number of measured species) = {0}'.format(M))
                print('The shape of data is {0}'.format(np.shape(data)))
            assert np.shape(data) == (N,T,M)
        elif type(exp_data) is pd.DataFrame:
            # Extract time
            if self.time_column:
                self.timepoints = np.array(exp_data.get(self.time_column), dtype='double').flatten()
            else:
                raise TypeError('time_column attribute of InferenceSetup object must be a string.')

            # Extract measurements
            if type(self.measurements) is list and len(self.measurements) == 1:
                data = np.array(exp_data.get(self.measurements[0]), dtype='double')
            elif type(self.measurements) is list and len(self.measurements) > 1:
                data_list = []
                for m in self.measurements:
                    data_list.append(np.array(exp_data.get(m), dtype='double'))
                data = np.array(data_list)
            else:
                raise ValueError('Something wrong with experimental data input to inference.')
            N = 1 # Number of trajectories
            T = len(self.timepoints) # Number of timepoints
            M = len(self.measurements)# Number of measurements
            data = np.reshape(data, (N,T,M))
        else:
            raise TypeError('exp_data attribute of InferenceSetup object must'
                            'be a list of Pandas DataFrames or a single Pandas DataFrame. ')
        return data

    def setup_cost_function(self, **kwargs):
        """
        Build the `PIDInterface`-derived likelihood for `self.sim_type`.

        Constructs a `~bioscrape.pid_interfaces.StochasticInference` or
        `~bioscrape.pid_interfaces.DeterministicInference` (depending
        on `self.sim_type`), and sets it up against `self.LL_data`.
        Stores the result as `self.pid_interface`, used by
        `cost_function`.

        Parameters
        ----------
        **kwargs
            Forwarded to the `PIDInterface` constructor and to
            `~bioscrape.pid_interfaces.PIDInterface.setup_likelihood_function`.
        """
        if self.sim_type == 'stochastic':
            self.pid_interface = StochasticInference(self.params_to_estimate, self.M, self.prior, **kwargs)
            self.pid_interface.setup_likelihood_function(self.LL_data, self.timepoints, self.measurements,
                                                         initial_conditions=self.initial_conditions,
                                                         parameter_conditions=self.parameter_conditions,
                                                         norm_order = self.norm_order,
                                                         N_simulations = self.N_simulations, **kwargs)
        elif self.sim_type == 'deterministic':
            self.pid_interface = DeterministicInference(self.params_to_estimate, self.M, self.prior, **kwargs)
            self.pid_interface.setup_likelihood_function(self.LL_data, self.timepoints, self.measurements,
                                                         initial_conditions=self.initial_conditions,
                                                         parameter_conditions=self.parameter_conditions,
                                                         norm_order=self.norm_order, **kwargs)

    def cost_function(self, params):
        """
        The log-likelihood of `params`, for use as an MCMC cost function.

        Evaluates `self.pid_interface`'s likelihood at `params`, and
        records both the value and `params` (in `self.cost_progress`/
        `self.cost_params`) for later inspection.

        Parameters
        ----------
        params : array_like
            The parameter values, in the order of
            `self.params_to_estimate`.

        Returns
        -------
        float
            The log-likelihood of `params`.

        Raises
        ------
        RuntimeError
            If `setup_cost_function` has not been called yet.
        """
        if self.pid_interface is None:
            raise RuntimeError("Must call InferenceSetup.setup_cost_function() \
                               before InferenceSetup.cost_function(params) can be used.")
        cost_value = self.pid_interface.get_likelihood_function(params)
        self.cost_progress.append(cost_value)
        self.cost_params.append(params)
        return cost_value

    def seed_parameter_values(self, **kwargs):
        """
        Sample the initial parameter values for each MCMC walker.

        Uses `self.init_seed` (optionally overridden by the `init_seed`
        keyword argument) to determine the initial values; see the
        class docstring for the supported forms. Values whose prior
        has a 'positive' flag are clipped to be non-negative (or a
        small positive epsilon, in log space).

        Parameters
        ----------
        init_seed : float or numpy.ndarray or list or 'prior', optional
            Overrides `init_seed`; see `set_init_seed`.

        Returns
        -------
        numpy.ndarray
            Array of shape `(self.nwalkers, len(self.params_to_estimate))`
            giving the initial parameter values for each walker.

        Raises
        ------
        ValueError
            If `init_seed` is not one of the supported forms.
        """
        if 'init_seed' in kwargs:
            self.set_init_seed(kwargs['init_seed'])
        ndim = len(self.params_to_estimate)
        params_values = []
        for p in self.params_to_estimate:
            value = self.M.get_parameter_dictionary()[p]
            params_values.append(value)
        #Sample a one percent ball around a given initial value
        if (isinstance(self.init_seed, np.ndarray) \
            or isinstance(self.init_seed, list)) \
            and len(self.init_seed) == ndim:
            p0 = np.array(self.init_seed) + 0.01*np.array(self.init_seed) \
                * np.random.randn(self.nwalkers, ndim)
        #Use this exact start value
        elif isinstance(self.init_seed, np.ndarray) \
            and self.init_seed.shape == (self.nwalkers, ndim):
            p0 =  np.array(self.init_seed)
        #Sample the Prior Distributions to determine initial values
        elif self.init_seed == "prior":
            p0 = np.zeros((self.nwalkers, ndim))
            for i, p in enumerate(self.params_to_estimate):
                prior = self.prior[p]
                if prior[0] == "uniform":
                    p0[:, i] = np.random.rand(self.nwalkers)*(prior[2]-prior[1])+prior[1]
                elif prior[0] == "gaussian":
                    p0[:, i] = prior[2]*np.random.randn(self.nwalkers)+prior[1]
                elif prior[0] == "log-uniform":
                    a = np.log(prior[1])
                    b = np.log(prior[2])
                    u = np.random.randn(self.nwalkers)*(b - a)+a
                    p0[:, i] = np.exp(u)
                else:
                    raise ValueError("Can only sample uniform and gaussian priors"
                                     "when 'init_seed' is set to prior. "
                                     "Try setting intial seed to a number [0, 1]"
                                     "to sample a gaussian ball around the model"
                                     "parameters instead.")
        #sample a gaussian ball around the initial model parameters
        elif isinstance(self.init_seed, float):
            p0 = np.array(params_values) + self.init_seed * np.array(params_values) * np.random.randn(self.nwalkers, ndim)
        else:
            raise ValueError("init_seed must be a float (will sample a gaussian ball"
                             " of this percent around the model initial condition), "
                             "array (of size parameters or walkers x parameters), "
                             "or the string 'prior' (will sample from uniform "
                             "and guassian priors)")

        #Ensure parameters are positive, if their priors are declared to be positive
        #When working in log space, a small non-zero value is used
        if hasattr(self.pid_interface, "log_space_parameters") and self.pid_interface.log_space_parameters:
            epsilon = 10**-8
        else:
            epsilon = 0

        for i, p in enumerate(self.params_to_estimate):
            if "positive" in self.prior[p]:
                p0[:, i] = p0[:, i]*(p0[:, i] > 0) + (p0[:, i] <= 0)*epsilon


        #convert to log space, if pid_interface.log_space_parameters
        if hasattr(self.pid_interface, "log_space_parameters") and self.pid_interface.log_space_parameters:
            p0 = np.log(p0)
        return p0

    def run_emcee(self, **kwargs):
        """
        Run `emcee.EnsembleSampler` and write out the results.

        Lower-level than `run_mcmc`: assumes `prepare_inference` has
        already been called (via `run_mcmc`, or `reuse_likelihood=True`
        to reuse an existing `pid_interface`). See `run_mcmc` for the
        accepted keyword arguments.

        Parameters
        ----------
        **kwargs
            See `run_mcmc`; also accepts `reuse_likelihood` (bool,
            default False) to skip rebuilding the likelihood via
            `setup_cost_function`.

        Returns
        -------
        emcee.EnsembleSampler
            The sampler, after running `self.nsteps` steps.

        Raises
        ------
        ImportError
            If the `emcee` package is not installed, or if
            `self.parallel` is True and the `multiprocessing` package
            is not available.
        """
        if kwargs.get("reuse_likelihood", False) is False:
            self.setup_cost_function(**kwargs)
        progress = kwargs.get('progress')
        convergence_check = kwargs.get('convergence_check', False)
        convergence_diagnostics = kwargs.get('convergence_diagnostics', convergence_check)
        skip_initial_state_check = kwargs.get('skip_initial_state_check', False)
        progress = kwargs.get('progess', True)
        fname_csv = kwargs.get('filename_csv', 'mcmc_results.csv')
        if 'results_filename' in kwargs:
            warnings.warn('The keyword results_filename is deprecated and'
                          'is being replaced by filename_csv and filename_txt where CSV is for'
                          'the MCMC samples and cost function progress respectively.', DeprecationWarning)
            fname_csv = kwargs.get('results_filename', 'mcmc_results.csv')
        fname_txt = kwargs.get('filename_txt', 'mcmc_results.txt')
        printout = kwargs.get('printout', True)

        try:
            import emcee
        except:
            raise ImportError('emcee package not installed.')
        ndim = len(self.params_to_estimate)
        p0 = self.seed_parameter_values(**kwargs)
        assert p0.shape == (self.nwalkers, ndim)
        if self.parallel:
            try:
                import multiprocessing
                pool = multiprocessing.Pool()
                if printout: print("Using {} cores for parallelization".format(multiprocessing.cpu_count()))
            except:
                pool = None
                raise ImportError('multiprocessing package not found. \
                                  Make sure to set parallel=False')
        else:
            pool = None
            if printout: print("creating an ensemble sampler without multiprocessing "\
                               "pool. Set parallel=True to use parallel processing.")
        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.cost_function, pool=pool)
        sampler.run_mcmc(p0, self.nsteps, progress=progress,
                         skip_initial_state_check=skip_initial_state_check)
        if self.parallel:
            pool.close()
            pool.join()
        if convergence_check:
            self.autocorrelation_time = sampler.get_autocorr_time()
        if convergence_diagnostics:
            if not convergence_check:
                warnings.warn('MCMC diagnostics cannot be printed when convergence check is False.')
                self.convergence_diagnostics = {}
            else:
                self.convergence_diagnostics = {'Autocorrelation time for each parameter':self.autocorrelation_time,
                                    'Acceptance fraction (fraction of steps that were accepted)':sampler.acceptance_fraction}
        # Write results
        import csv
        with open(fname_csv,'w', newline = "") as f:
            writer = csv.writer(f)
            writer.writerows(sampler.get_chain(flat = True))
            f.close()
        with open(fname_txt, 'w', newline = "") as f:
            f.write('\nCost function progress\n')
            f.write(str(self.cost_progress))
            if convergence_diagnostics:
                f.write('\nMCMC convergence diagnostics\n')
                f.write(str(self.convergence_diagnostics))
            f.close()
        if printout: print("Results written to" + fname_csv + " and " + fname_txt)
        if printout: print('Successfully completed MCMC parameter identification procedure. '
                           'Check the MCMC diagnostics to evaluate convergence.')
        return sampler

    def plot_mcmc_results(self, sampler, **kwargs):
        """
        Plot MCMC chain and posterior distributions for a sampled run.

        Parameters
        ----------
        sampler : emcee.EnsembleSampler
            A sampler that has already been run (e.g. the return value
            of `run_mcmc`).
        figsize : tuple, default (10, 7)
            Passed to `matplotlib.pyplot.subplots` for the chain plot.
        sharex : bool, default True
            Passed to `matplotlib.pyplot.subplots` for the chain plot.
        labels : list of str, optional
            Axis labels for each parameter; defaults to
            `self.params_to_estimate`.
        alpha : float, default 0.3
            Transparency of the individual walker traces in the chain
            plot.
        discard : int, default `2 * self.nwalkers`
            Number of initial steps to discard as burn-in before
            computing the posterior summary and corner plot.
        thin : int, default 1
            Keep only every `thin`-th post-burn-in step.
        flat : bool, default True
            Passed to `sampler.get_chain` when computing the posterior
            summary.
        percentiles : list of float, default [16, 50, 84]
            Percentiles of the post-burn-in samples to report for each
            parameter (as `param_report`).

        Returns
        -------
        param_report : dict
            Maps each parameter name plus ``'_true'`` to its 50th
            percentile value, and plus ``'_uncertainties'`` to the
            differences between consecutive `percentiles` (by default,
            the distances from the 16th/84th percentiles to the
            median).
        figure_objects : list
            `[fig, axes, corner_fig]`, where `fig`/`axes` are the
            per-parameter chain plot and `corner_fig` is the
            `corner.corner` posterior plot -- or just `[fig, axes]` if
            the `corner` package is not installed.
        """
        print('Parameter posterior distribution convergence plots:')
        figsize = kwargs.get('figsize', (10,7))
        sharex = kwargs.get('sharex', True)
        ndim = sampler.ndim
        fig, axes = plt.subplots(ndim, figsize=figsize, sharex=sharex)
        figure_objects = []
        samples = sampler.get_chain()
        labels = kwargs.get('labels', list(self.params_to_estimate))
        alpha = kwargs.get('alpha', 0.3)
        for i in range(ndim):
            if type(axes) is np.ndarray:
                ax = axes[i]
            else:
                ax = axes
            ax.plot(samples[:, :, i], "k", alpha=alpha)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
        if type(axes) is np.ndarray:
            axes[-1].set_xlabel("step number")
        else:
            axes.set_xlabel("step number")
        figure_objects.append(fig)
        figure_objects.append(axes)
        # arbitrarily discard 2*nwalkers steps
        discard = kwargs.get('discard', 2*self.nwalkers)
        thin = int(kwargs.get('thin', 1))
        flat = kwargs.get('flat', True)

        flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=flat)
        param_report = {}
        # Percentiles to compute for numpy.percentile
        percentiles = kwargs.get('percentiles', [16, 50, 84]) # Set percentiles to compute by default to q
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], q = percentiles)
            q = np.diff(mcmc)
            param_report[self.params_to_estimate[i] + '_true'] = mcmc[1]
            param_report[self.params_to_estimate[i] + '_uncertainties'] = q
            # uncertainty_list.append([-1.0*q[0], q[1]])
            # uncertainty_list.append(q)
        try:
            import corner
            corner_fig = corner.corner(
                flat_samples, labels=labels,
            )
            figure_objects.append(corner_fig)
        except:
            warnings.warn('corner package not found - cannot plot parameter distributions.')
        return param_report, figure_objects


    def run_lmfit(self, method = 'leastsq', plot_show = True, **kwargs):
        """
        Run least-squares parameter inference via `lmfit`.

        Fits each trajectory in `self.exp_data` independently with
        `lmfit.minimize`, using `~bioscrape.pid_interfaces.LMFitInference`.

        Parameters
        ----------
        method : str, default 'leastsq'
            Passed to `lmfit.minimize`; see the
            `lmfit docs <https://lmfit.github.io/lmfit-py/fitting.html>`_.
        plot_show : bool, default True
            If True, display the fit for each trajectory.
        **kwargs
            Forwarded to `~bioscrape.pid_interfaces.LMFitInference.get_minimizer_results`.

        Returns
        -------
        list of lmfit.minimizer.MinimizerResult
            One result per trajectory in `self.exp_data`.
        """
        self.prepare_inference(**kwargs)
        N = np.shape(self.LL_data)[0]
        if self.sim_type == 'stochastic':
            stochastic = True
        else:
            stochastic = False
        self.pid_interface = LMFitInference(self.params_to_estimate, self.M, self.prior)
        minimizer_result = [None]*N
        if N == 1:
            minimizer_result[0] = self.pid_interface.\
                get_minimizer_results(self.LL_data[0,:,:],
                                      self.timepoints, self.measurements,
                                      initial_conditions = self.initial_conditions,
                                      parameter_conditions = self.parameter_conditions,
                                      stochastic = stochastic, debug = self.debug,
                                      method = method, plot_show = plot_show, **kwargs)
        else:
            if self.parameter_conditions is None:
                self.parameter_conditions = [None]*N
            for i in range(N):
                minimizer_result[i] = self.pid_interface.\
                    get_minimizer_results(self.LL_data[i,:,:],
                                          self.timepoints[i], self.measurements,
                                          initial_conditions = self.initial_conditions[i],
                                          parameter_conditions = self.parameter_conditions[i],
                                          stochastic = stochastic, debug = self.debug,
                                          method = method, plot_show = plot_show, **kwargs)
        print('Successfully completed parameter identification'
              ' procedure using LMFit. Parameter values and fitness'
              ' reports written to lmfit_results.csv file. '
              'Check the minimizer_results object returned to '
              'further statistically evaluate the goodness of fit.')
        return minimizer_result

    def write_lmfit_results(self, minimizer_result, **kwargs):
        """
        Write `run_lmfit` results to a CSV file.

        Writes each trajectory's fitted parameter values to
        ``lmfit_results.csv`` in the current directory (one row per
        trajectory), optionally followed by each trajectory's full
        `lmfit` convergence report.

        Parameters
        ----------
        minimizer_result : list of lmfit.minimizer.MinimizerResult
            The result of `run_lmfit`.
        convergence_diagnostics : bool, default True
            If True, append each trajectory's `lmfit.fit_report` to
            ``lmfit_results.csv``.

        Raises
        ------
        ImportError
            If the `lmfit` package is not installed.
        """
        import csv
        try:
            from lmfit import fit_report
        except:
            raise ImportError('Package lmfit not found.')
        values_dict = {}
        for param_name in self.params_to_estimate:
            values_dict[param_name] = []
        for result in minimizer_result:
            for param_name in self.params_to_estimate:
                values_dict[param_name].append(dict(result.params.valuesdict())[param_name])
        df = pd.DataFrame.from_dict(values_dict)
        # df.to_csv('lmfit_results.csv', columns = self.params_to_estimate, index = False)
        with open('lmfit_results.csv','w') as f:
            f.write(str(df))

        convergence_diagnostics = kwargs.get('convergence_diagnostics', True)
        if convergence_diagnostics:
            count = 0
            for result in minimizer_result:
                with open('lmfit_results.csv','a') as f:
                    # writer = csv.writer(f)
                    f.write('\nFor trajectory: {0}\n'.format(count))
                    f.write(fit_report(result))
                    count += 1
                    f.close()
