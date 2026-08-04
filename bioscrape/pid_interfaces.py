"""
Parameter identification (PID) interfaces for bioscrape models.

Provides `PIDInterface` and built-in subclasses that connect a
`~bioscrape.types.Model`, a set of parameters to estimate, and a prior
specification to the appropriate likelihood function, for use with
`~bioscrape.inference.py_inference`.
"""
from bioscrape.inference import DeterministicLikelihood as DLL
from bioscrape.inference import StochasticTrajectoriesLikelihood as STLL
from bioscrape.inference import StochasticTrajectories
from bioscrape.inference import BulkData
from bioscrape.simulator import py_simulate_model
import matplotlib.pyplot as plt
import warnings
import numpy as np

class PIDInterface():
    """
    Base class for parameter identification (PID) interfaces.

    Connects a `~bioscrape.types.Model`, a set of parameters to estimate, and
    a prior specification to a likelihood function, for use with
    `~bioscrape.inference.py_inference`. Two subclasses are built in:
    `StochasticInference` and `DeterministicInference`. New interfaces can be
    added by subclassing `PIDInterface` with the desired log-likelihood
    functions; a subclass may also override `check_prior` to use custom priors
    instead of the built-in ones below.

    Parameters
    ----------
    params_to_estimate : list of str
        Names of the parameters to be estimated.
    M : Model
        The bioscrape model to use for inference.
    prior : dict
        A dictionary specifying the prior distribution for each parameter in
        `params_to_estimate`. Each entry is a list of the form
        `['prior_name', *distribution_params]`, e.g.
        `['uniform', lower_bound, upper_bound]` or
        `['gaussian', mean, standard_deviation]`. See `check_prior` for the
        full list of built-in prior types and their parameters.
    log_space_parameters : bool, default False
        Whether the parameters being estimated are in log space.
    debug : bool, default False
        If True, print additional diagnostic information during
        inference.
    """
    def __init__(self, params_to_estimate, M, prior, **kwargs):
        """See class docstring."""
        self.params_to_estimate = params_to_estimate
        self.M = M
        self.prior = prior
        self.default_parameters = dict(M.get_parameter_dictionary())
        self.log_space_parameters = kwargs.get('log_space_parameters', False)
        self.debug = kwargs.get('debug', False)
        return

    def check_prior(self, params_dict):
        """
        Compute the total log-prior probability for a parameter set.

        Dispatches each parameter in `params_dict` to the built-in prior
        function matching its type in `self.prior` (e.g. `uniform_prior`,
        `gaussian_prior`) and sums the results. A 'custom' prior type may
        supply its own callable directly in the prior dictionary (see
        `PIDInterface`); for anything more involved, override `check_prior` in
        a subclass.

        Parameters
        ----------
        params_dict : dict
            Dictionary mapping parameter name to a candidate value.

        Returns
        -------
        float
            The total log-prior probability, or `numpy.inf` if any parameter
            has a 'positive' prior flag and a negative value, or if any
            individual prior evaluates to an invalid probability.
        """
        lp = 0.0
        for key,value in params_dict.items():
            if 'positive' in self.prior[key] and value  < 0:
                return np.inf
            prior_type = self.prior[key][0]
            if prior_type == 'uniform':
                lp += self.uniform_prior(key, value)
            elif prior_type == 'gaussian':
                lp += self.gaussian_prior(key, value)
            elif prior_type == 'exponential':
                lp += self.exponential_prior(key, value)
            elif prior_type == 'gamma':
                lp += self.gamma_prior(key, value)
            elif prior_type == 'log-uniform':
                lp += self.log_uniform_prior(key, value)
            elif prior_type == 'log-gaussian':
                lp += self.log_gaussian_prior(key, value)
            elif prior_type == 'beta':
                lp += self.beta_prior(key, value)
            elif prior_type == 'custom':
                # The last element in the prior dictionary must be a callable function
                # The callable function shoud have the following signature :
                # Arguments: param_name (str), param_value(float)
                # Returns: log prior probability (float or numpy inf)
                custom_fuction = self.prior[key][-1]
                lp += custom_fuction(key, value)
            else:
                raise ValueError(f'Prior type undefined: recieved prior {value[0]} for param {key}.')
        return lp

    def uniform_prior(self, param_name, param_value):
        """
        Log-prior probability of `param_value` under a uniform prior.

        Parameters
        ----------
        param_name : str
            The parameter name; used to look up the prior's lower and upper
            bounds in `self.prior`.
        param_value : float
            The candidate parameter value.

        Returns
        -------
        float
            `numpy.log(1 / (upper_bound - lower_bound))` if `param_value` is
            within the prior's bounds, otherwise `numpy.inf` (used as a
            rejection sentinel by `check_prior`).

        Raises
        ------
        ValueError
            If no prior was specified for this `PIDInterface` (`self.prior`
            is None).
        """
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        lower_bound = prior_dict[param_name][1]
        upper_bound = prior_dict[param_name][2]
        if param_value > upper_bound or param_value < lower_bound:
            return np.inf
        else:
            return np.log( 1/(upper_bound - lower_bound) )

    def gaussian_prior(self, param_name, param_value):
        """
        Log-prior probability of `param_value` under a Gaussian prior.

        Parameters
        ----------
        param_name : str
            The parameter name; used to look up the prior's mean and standard
            deviation in `self.prior`.
        param_value : float
            The candidate parameter value.

        Returns
        -------
        float
            The log of the Gaussian probability density (with the mean and
            standard deviation from `self.prior`) at `param_value`.

        Raises
        ------
        ValueError
            If no prior was specified for this `PIDInterface` (`self.prior`
            is None), or if the prior's standard deviation is negative.
        """
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        mu = prior_dict[param_name][1]
        sigma = prior_dict[param_name][2]
        if sigma < 0:
            raise ValueError('The standard deviation must be positive.')
        # Using probability density function for normal distribution
        # Using scipy.stats.norm has overhead that affects speed up to 2x
        prob = 1/(np.sqrt(2*np.pi) * sigma) * np.exp(-0.5*(param_value - mu)**2/sigma**2)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking Gaussian prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

    def exponential_prior(self, param_name, param_value):
        """
        Log-prior probability of `param_value` under an exponential prior.

        Parameters
        ----------
        param_name : str
            The parameter name; used to look up the prior's rate parameter in
            `self.prior`.
        param_value : float
            The candidate parameter value.

        Returns
        -------
        float
            The log of the exponential probability density at `param_value`,
            or `numpy.inf` (with a warning) if the computed probability is
            negative.

        Raises
        ------
        ValueError
            If no prior was specified for this `PIDInterface` (`self.prior`
            is None).
        """
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        lambda_p = prior_dict[param_name][1]

        prob = lambda_p * np.exp(-lambda_p * param_value)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking Exponential prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

    def gamma_prior(self, param_name, param_value):
        """
        Log-prior probability of `param_value` under a gamma prior.

        Parameters
        ----------
        param_name : str
            The parameter name; used to look up the prior's shape (`alpha`)
            and rate (`beta`) parameters in `self.prior`.
        param_value : float
            The candidate parameter value.

        Returns
        -------
        float
            The log of the gamma probability density at `param_value`, or
            `numpy.inf` (with a warning) if the computed probability is
            negative.

        Raises
        ------
        ValueError
            If no prior was specified for this `PIDInterface` (`self.prior`
            is None).
        """
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        alpha = prior_dict[param_name][1]
        beta = prior_dict[param_name][2]
        from scipy.special import gamma
        prob = (beta**alpha)/gamma(alpha) * param_value**(alpha - 1) * np.exp(-1 * beta*param_value)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking Exponential prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

    def beta_prior(self, param_name, param_value):
        """
        Log-prior probability of `param_value` under a beta prior.

        Parameters
        ----------
        param_name : str
            The parameter name; used to look up the prior's `alpha` and `beta`
            shape parameters in `self.prior`.
        param_value : float
            The candidate parameter value (expected in `[0, 1]`).

        Returns
        -------
        float
            The log of the beta probability density at `param_value`, or
            `numpy.inf` (with a warning) if the computed probability is
            negative.

        Raises
        ------
        ValueError
            If no prior was specified for this `PIDInterface` (`self.prior`
            is None).
        """
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        alpha = prior_dict[param_name][1]
        beta = prior_dict[param_name][2]
        from scipy import special
        prob = (param_value**(alpha-1) * (1 - param_value)**(beta - 1) )/special.beta(alpha, beta)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking Exponential prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

    def log_uniform_prior(self, param_name, param_value):
        """
        Log-prior probability of `param_value` under a log-uniform prior.

        Parameters
        ----------
        param_name : str
            The parameter name; used to look up the prior's lower and upper
            bounds in `self.prior`.
        param_value : float
            The candidate parameter value.

        Returns
        -------
        float
            The log of the log-uniform probability density if `param_value` is
            within the prior's (positive) bounds, otherwise `numpy.inf` (used
            as a rejection sentinel by `check_prior`).

        Raises
        ------
        ValueError
            If no prior was specified for this `PIDInterface` (`self.prior`
            is None), or if the prior's lower or upper bound is negative.
        """
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        lower_bound = prior_dict[param_name][1]
        upper_bound = prior_dict[param_name][2]

        if lower_bound < 0 or upper_bound < 0:
            raise ValueError('Upper and lower bounds for log-uniform prior must be positive.')

        if param_value > upper_bound or param_value < lower_bound:
            return np.inf

        prob = 1/(param_value* (np.log(upper_bound) - np.log(lower_bound)))
        if prob < 0:
            warnings.warn('Probability less than 0 while checking Log-Uniform prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

    def log_gaussian_prior(self, param_name, param_value):
        """
        Log-prior probability of `param_value` under a log-normal prior.

        Parameters
        ----------
        param_name : str
            The parameter name; used to look up the prior's mean and standard
            deviation (of the underlying normal distribution) in
            `self.prior`.
        param_value : float
            The candidate parameter value.

        Returns
        -------
        float
            The log of the log-normal probability density at `param_value`.

        Raises
        ------
        ValueError
            If no prior was specified for this `PIDInterface` (`self.prior`
            is None), or if the prior's standard deviation is negative.
        """
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        mu = prior_dict[param_name][1]
        sigma = prior_dict[param_name][2]
        if sigma < 0:
            raise ValueError('The standard deviation must be positive.')
        # Using probability density function for log-normal distribution
        prob = 1/(param_value * np.sqrt(2*np.pi) * sigma) * np.exp((-0.5 * (np.log(param_value) - mu)**2)/sigma**2)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking log-normal prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

# Add a new class similar to this to create new interfaces.
class StochasticInference(PIDInterface):
    """
    PID interface for Bayesian inference of stochastic models.

    Fits `~bioscrape.types.Model` parameters to single-cell trajectory data
    (which may come from experiments or from another simulation) using
    `~bioscrape.inference.StochasticTrajectoriesLikelihood`, comparing
    simulated and measured trajectory statistics rather than individual
    trajectories.

    Parameters
    ----------
    params_to_estimate : list of str
        Names of the parameters to be estimated.
    M : Model
        The bioscrape model to use for inference.
    prior : dict
        A dictionary specifying the prior distribution for each parameter;
        see `PIDInterface`.
    log_space_parameters : bool, default False
        Whether the parameters being estimated are in log space; see
        `PIDInterface`.
    debug : bool, default False
        If True, print additional diagnostic information during
        inference.
    """
    def __init__(self, params_to_estimate, M, prior, **kwargs):
        """See class docstring."""
        self.LL_stoch = None
        self.dataStoch = None
        if 'debug' in kwargs:
            self.debug = kwargs.get('debug')
        super().__init__(params_to_estimate, M, prior, **kwargs)
        return

    def setup_likelihood_function(self, data, timepoints, measurements,
                                  initial_conditions, parameter_conditions,
                                  norm_order=2, N_simulations=3,
                                  **kwargs):
        """
        Build the likelihood function used by `get_likelihood_function`.

        Must be called once before `get_likelihood_function`. Wraps `data` in
        a `~bioscrape.inference.StochasticTrajectories` object and builds a
        `~bioscrape.inference.StochasticTrajectoriesLikelihood` from it,
        comparing `N_simulations` stochastic simulations per data point
        against the trajectory statistics in `data`.

        Parameters
        ----------
        data : array_like
            Experimental trajectory data, of shape
            `(N, len(timepoints), len(measurements))` for `N` trajectories.
        timepoints : array_like
            The time points at which `data` was measured.
        measurements : list of str
            The names of the measured species/outputs in `data`.
        initial_conditions : dict or list of dict
            Initial condition(s) for the simulated trajectories. A single dict
            if `data` has one trajectory, or a list of `N` dicts for multiple
            trajectories.
        parameter_conditions : dict or list of dict or None
            Per-trajectory parameter overrides, in the same shape as
            `initial_conditions`, or None to use the model's default
            parameters for every trajectory.
        norm_order : int, optional
            The norm order used to compare simulated and experimental
            trajectory statistics (default 2).
        N_simulations : int, optional
            The number of stochastic simulations to average per likelihood
            evaluation (default 3).
        **kwargs
            Additional keyword arguments passed to
            `~bioscrape.inference.StochasticTrajectoriesLikelihood`.
        """
        N = np.shape(data)[0]
        self.dataStoch = StochasticTrajectories(np.array(timepoints), data, measurements, N)
        if self.debug:
            print('Stochastic inference attributes:')
            print('The timepoints shape is {0}'.format(np.shape(timepoints)))
            print('The data shape is {0}'.format(np.shape(data)))
            print('The measurmenets is {0}'.format(measurements))
            print('The N is {0}'.format(N))
            print('Using the initial conditions: {0}'.format(initial_conditions))
            print('Using the parameter conditions: {0}'.format(parameter_conditions))
        #If there are multiple initial conditions in a data-set,
        # should correspond to multiple initial conditions for inference.
        # Note len(initial_conditions) must be equal to the number of trajectories N
        # Same holds for parameter_conditions
        if parameter_conditions is not None:
            self.LL_stoch = STLL(model = self.M, init_state = initial_conditions,
                                init_params = parameter_conditions,
                                data = self.dataStoch, N_simulations = N_simulations,
                                norm_order = norm_order, **kwargs)
        else:
            self.LL_stoch = STLL(model = self.M, init_state = initial_conditions,
                                data = self.dataStoch, N_simulations = N_simulations,
                                norm_order = norm_order, **kwargs)


    def get_likelihood_function(self, params):
        """
        Compute the log posterior probability for a parameter set.

        Combines the log-prior from `check_prior` with the log-likelihood from
        the `~bioscrape.inference.StochasticTrajectoriesLikelihood` built by
        `setup_likelihood_function`. Intended to be passed as the
        log-probability function for an MCMC sampler such as `emcee`.

        Parameters
        ----------
        params : array_like
            Candidate values for the parameters in `self.params_to_estimate`,
            in the same order. If `self.log_space_parameters` is True,
            these are interpreted as log-space values and exponentiated
            before use.

        Returns
        -------
        float
            The log posterior probability, or `-numpy.inf` if the prior
            rejects `params`.

        Raises
        ------
        RuntimeError
            If called before `setup_likelihood_function`.
        """
        # Set params here and return the likelihood object.
        if self.LL_stoch is None:
            raise RuntimeError("Must call StochasticInference.setup_likelihood_function before using StochasticInference.get_likelihood_function.")

        #Set params
        params_dict = {}
        for key, p in zip(self.params_to_estimate, params):
            if self.log_space_parameters:
                params_dict[key] = np.exp(p)
            else:
                params_dict[key] = p

        #Prior
        lp = 0
        lp = self.check_prior(params_dict)
        if not np.isfinite(lp):
            return -np.inf
        else:
            # Reset to default
            self.LL_stoch.set_init_params(self.default_parameters)
            self.LL_stoch.set_init_params(params_dict)
            if self.debug:
                print('current sample:', params_dict)
            LL_stoch_cost = self.LL_stoch.py_log_likelihood()
            ln_prob = lp + LL_stoch_cost
            if self.debug:
                print('current cost total:', ln_prob)
            return ln_prob

# Add a new class similar to this to create new interfaces.
class DeterministicInference(PIDInterface):
    """
    PID interface for Bayesian inference of deterministic models.

    Fits `~bioscrape.types.Model` parameters to bulk time-series data
    (which may come from experiments or from another simulation) using
    `~bioscrape.inference.DeterministicLikelihood`.

    Parameters
    ----------
    params_to_estimate : list of str
        Names of the parameters to be estimated.
    M : Model
        The bioscrape model to use for inference.
    prior : dict
        A dictionary specifying the prior distribution for each parameter;
        see `PIDInterface`.
    log_space_parameters : bool, default False
        Whether the parameters being estimated are in log space; see
        `PIDInterface`.
    debug : bool, default False
        If True, print additional diagnostic information during
        inference.
    """
    def __init__(self, params_to_estimate, M, prior, **kwargs):
        """See class docstring."""
        self.LL_det = None
        self.dataDet = None
        self.debug = None
        if 'debug' in kwargs:
            self.debug = kwargs.get('debug')
        super().__init__(params_to_estimate, M, prior, **kwargs)
        return

    def setup_likelihood_function(self, data, timepoints, measurements,
                                  initial_conditions, parameter_conditions,
                                  norm_order = 2, **kwargs):
        """
        Build the likelihood function used by `get_likelihood_function`.

        Must be called once before `get_likelihood_function`. Wraps `data` in
        a `~bioscrape.inference.BulkData` object and builds a
        `~bioscrape.inference.DeterministicLikelihood` from it.

        Parameters
        ----------
        data : array_like
            Experimental trajectory data, of shape
            `(N, len(timepoints), len(measurements))` for `N` trajectories.
        timepoints : array_like
            The time points at which `data` was measured.
        measurements : list of str
            The names of the measured species/outputs in `data`.
        initial_conditions : dict or list of dict
            Initial condition(s) for the simulated trajectories. A single dict
            if `data` has one trajectory, or a list of `N` dicts for multiple
            trajectories.
        parameter_conditions : dict or list of dict or None
            Per-trajectory parameter overrides, in the same shape as
            `initial_conditions`, or None to use the model's default
            parameters for every trajectory.
        norm_order : int, optional
            The norm order used to compare simulated and experimental
            trajectories (default 2).
        **kwargs
            Additional keyword arguments passed to
            `~bioscrape.inference.DeterministicLikelihood`.
        """
        N = np.shape(data)[0]
        #Create a data Objects
        # In this case the timepoints should be a list of timepoints vectors for each iteration
        self.dataDet = BulkData(np.array(timepoints), data, measurements, N)
        #If there are multiple initial conditions in a data-set,
        # should correspond to multiple initial conditions for inference.
        #Note len(initial_conditions) must be equal to the number of trajectories N
        # Similarly, if parameter_conditions are provided then the length must equal
        # number of trajectories N
        if self.debug:
            print('The deterministic inference attributes:')
            print('The timepoints shape is {0}'.format(np.shape(timepoints)))
            print('The data shape is {0}'.format(np.shape(data)))
            print('The measurmenets is {0}'.format(measurements))
            print('The N is {0}'.format(N))
            print('Using the initial conditions: {0}'.format(initial_conditions))
            print('Using the parameter conditions: {0}'.format(parameter_conditions))
        #Create Likelihood object
        if parameter_conditions is not None:
            self.LL_det = DLL(model = self.M, init_state = initial_conditions,
                              init_params = parameter_conditions,
                              data = self.dataDet, norm_order = norm_order, **kwargs)
        else:
            self.LL_det = DLL(model = self.M, init_state = initial_conditions,
                              data = self.dataDet, norm_order = norm_order, **kwargs)

    def get_likelihood_function(self, params):
        """
        Compute the log posterior probability for a parameter set.

        Combines the log-prior from `check_prior` with the log-likelihood from
        the `~bioscrape.inference.DeterministicLikelihood` built by
        `setup_likelihood_function`. Intended to be passed as the
        log-probability function for an MCMC sampler such as `emcee`.

        Parameters
        ----------
        params : array_like
            Candidate values for the parameters in `self.params_to_estimate`,
            in the same order. If `self.log_space_parameters` is True,
            these are interpreted as log-space values and exponentiated
            before use.

        Returns
        -------
        float
            The log posterior probability, or `-numpy.inf` if the prior
            rejects `params`.

        Raises
        ------
        RuntimeError
            If called before `setup_likelihood_function`.
        """
        if self.LL_det is None:
            raise RuntimeError("Must call DeterministicInference.setup_likelihood_function before using DeterministicInference.get_likelihood_function.")
        #this part is the only part that is called repeatedly
        params_dict = {}
        for key, p in zip(self.params_to_estimate, params):
            if self.log_space_parameters:
                params_dict[key] = np.exp(p)
            else:
                params_dict[key] = p

        # Check prior
        lp = 0
        lp = self.check_prior(params_dict)
        if not np.isfinite(lp):
            return -np.inf
        else:
            # Reset to default
            self.LL_det.set_init_params(self.default_parameters)
            # Set new sampler parameter
            self.LL_det.set_init_params(params_dict)
            if self.debug:
                print('current sample:', params_dict)
            #apply cost function
            LL_det_cost = self.LL_det.py_log_likelihood()
            # if self.debug:
            #     print('current cost:', LL_det_cost)
            ln_prob = lp + LL_det_cost
            if self.debug:
                print('current cost total:', ln_prob)
            #print("params", params, "ln_prob", ln_prob)
            return ln_prob

class LMFitInference(PIDInterface):
    """
    PID interface for least-squares/maximum-likelihood parameter fitting.

    Provides a point-estimate alternative to the MCMC-based
    `StochasticInference`/`DeterministicInference`, built on
    `lmfit <https://lmfit.github.io/lmfit-py/>`_, for cases where a full
    posterior is not needed.

    Parameters
    ----------
    params_to_estimate : list of str
        Names of the parameters to be estimated.
    M : Model
        The bioscrape model to use for inference.
    prior : dict
        A dictionary specifying the prior distribution for each parameter;
        used to set the lower/upper bounds passed to `lmfit`. See
        `PIDInterface`.
    log_space_parameters : bool, default False
        Whether the parameters being estimated are in log space; see
        `PIDInterface`.
    debug : bool, default False
        If True, print additional diagnostic information during
        inference.
    """
    def __init__(self, params_to_estimate, M, prior, **kwargs):
        """See class docstring."""
        self.residual_function = None
        self.dataLMFit = None
        super().__init__(params_to_estimate, M, prior, **kwargs)
        return

    def get_minimizer_results(self, data, timepoints, measurements, initial_conditions,
                              parameter_conditions, method = 'leastsq',
                              stochastic = False, debug = False,
                              plot_show = True, **kwargs):
        """
        Fit `self.params_to_estimate` to `data` via residual minimization.

        Builds an `lmfit.Parameters` object from `self.params_to_estimate`
        and `self.prior` (using each parameter's prior bounds, or a `[0, inf)`
        bound if its prior has a 'positive' flag), then calls `lmfit.minimize`
        on a residual function that simulates the model at each candidate
        parameter set and compares it against `data`. If a candidate parameter
        set is rejected by `check_prior`, the residual is set to NaN for that
        evaluation (handled via `nan_policy='propagate'`).

        Parameters
        ----------
        data : array_like
            Experimental trajectory data, of shape
            `(len(timepoints), len(measurements))`.
        timepoints : array_like
            The time points at which `data` was measured.
        measurements : list of str
            The names of the measured species/outputs in `data`.
        initial_conditions : dict
            Initial conditions for the simulated trajectory.
        parameter_conditions : dict or None
            Parameter overrides to apply before simulating, or None to use the
            model's default parameters.
        method : str, optional
            The `lmfit` minimization method to use (default 'leastsq'); see
            the `lmfit docs <https://lmfit.github.io/lmfit-py/fitting.html>`_
            for the full list.
        stochastic : bool, optional
            If True, simulate stochastically instead of deterministically
            (default False).
        debug : bool, optional
            If True, print additional diagnostic information (default False).
        plot_show : bool, optional
            If True, plot the fitted simulation against `data` for each
            measured species using matplotlib (default True).
        **kwargs
            Additional keyword arguments passed to `lmfit.minimize`.

        Returns
        -------
        lmfit.minimizer.MinimizerResult
            The result object returned by `lmfit.minimize`, including the
            fitted parameter values.
        """
        # In this case the timepoints should be a list of timepoints vectors for each iteration
        #If there are multiple initial conditions in a data-set,
        # should correspond to multiple initial conditions for inference.
        #Note len(initial_conditions) must be equal to the number of trajectories N
        # Same holds for parameter_conditions
        try:
            import lmfit
        except:
            raise ImportError('LMFit package not found.')
        if debug:
            print('The LMFit inference attributes:')
            print('The timepoints shape is {0}'.format(np.shape(timepoints)))
            print('The data shape is {0}'.format(np.shape(data)))
            print('The measurmenets is {0}'.format(measurements))
            print('The N is {0}'.format(N))
            print('Using the initial conditions: {0}'.format(initial_conditions))
            print('Using the parameter conditions: {0}'.format(parameter_conditions))

        def residual_function(params, data = data):
            # Integrate ODEs
            params_values_dict = {}
            for p in dict(params).keys():
                params_values_dict[params[p].name] = params[p].value
                self.M.set_parameter(p, params[p].value)
            # Check prior
            lp = 0
            lp = self.check_prior(params_values_dict)
            if not np.isfinite(lp):
                nans_array = np.array([np.nan]*len(timepoints))
                return nans_array
            self.M.set_species(initial_conditions)
            if parameter_conditions is not None:
                self.M.set_params(parameter_conditions)
            model_sim = py_simulate_model(timepoints, self.M,
                                          stochastic = stochastic,
                                          delay = self.M.has_delays())
            residual_value = np.zeros(len(timepoints))
            measurements_counter = 0
            for species in measurements:
                residual_value += lp + np.reshape(np.array(model_sim[species]), (len(timepoints))) - data[:,measurements_counter]
                measurements_counter += 1
            return residual_value
        model_param_dict = self.M.get_parameter_dictionary()
        # Create LMFit parameter objects
        params = lmfit.Parameters()
        for param_name in self.params_to_estimate:
            if self.prior[param_name][0] == 'uniform':
                param_min = self.prior[param_name][1]
                param_max = self.prior[param_name][2]
            elif self.prior[param_name][-1] == 'positive':
                param_min = 0
                param_max = np.inf
            else:
                param_min = -np.inf
                param_max = np.inf
            params[param_name] = lmfit.Parameter(name=param_name,
                                                 value=model_param_dict[param_name],
                                                 min = param_min, max = param_max)
        # Use the residual function above (that needs to be minimized)
        result = lmfit.minimize(residual_function, params, kws = {'data': data},
                                method = method, nan_policy = 'propagate', **kwargs)
        if plot_show:
            self.M.set_species(initial_conditions)
            if parameter_conditions is not None:
                self.M.set_params(parameter_conditions)
            self.M.set_params(dict(result.params.valuesdict()))
            model_sim_fit = py_simulate_model(timepoints, Model = self.M,
                                              stochastic = stochastic,
                                              delay = self.M.has_delays())
            measurements_counter = 0
            for species in measurements:
                plt.plot(timepoints, model_sim_fit[species], color = 'orange', alpha = 0.5)
                plt.plot(timepoints, data[:,measurements_counter], color = 'black', ls = 'dotted')
                measurements_counter += 1
            plt.xlabel('Time', fontsize = 14)
            plt.ylabel('Species', fontsize = 14)
        return result
