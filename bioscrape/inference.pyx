"""
Bayesian parameter inference for bioscrape models.

Provides `Data` (`BulkData`, `FlowData`, `StochasticTrajectories`) and
`Likelihood` (`DeterministicLikelihood`, `StochasticTrajectoriesLikelihood`,
`StochasticTrajectoryMomentLikelihood`, `StochasticStatesLikelihood`) class
hierarchies for comparing simulated and experimental data, along with the
user-facing `py_inference` entry point. See also
`~bioscrape.pid_interfaces`, which wraps these classes into higher-level
parameter identification interfaces.
"""
cimport numpy as np
import numpy as np
from libc.math cimport log

from types import Model
from types cimport Model
from simulator cimport CSimInterface, RegularSimulator, ModelCSimInterface, DeterministicSimulator, SSASimulator
from simulator import CSimInterface, RegularSimulator, ModelCSimInterface, DeterministicSimulator, SSASimulator
from simulator cimport DelaySimulator, DelaySSASimulator, ArrayDelayQueue, SSAResult
from simulator import DelaySimulator, DelaySSASimulator, ArrayDelayQueue, SSAResult
from inference_setup import initialize_inference
import sys

import emcee

##################################################                ####################################################
######################################              DISTRIBUTION                      ################################
#################################################                     ################################################


cdef class Distribution:
    """
    Base class for probability distributions used as priors.

    Subclasses implement `unprob` (unnormalized probability), `prob`
    (normalized probability), and `dim` (dimensionality); these are
    `cdef` methods, called from Python via `py_unprob`, `py_prob`, and
    `py_dim` below.
    """
    cdef double unprob(self,np.ndarray a):
        raise NotImplementedError("Calling unprob() on Distribution")
    cdef double prob(self,np.ndarray a):
        raise NotImplementedError("Calling prob() on Distribution")
    cdef unsigned dim(self):
        raise NotImplementedError("Calling dim() on Distribution")

    def py_unprob(self, np.ndarray a):
        """
        Unnormalized probability of a point.

        Parameters
        ----------
        a : numpy.ndarray
            The point at which to evaluate the distribution.

        Returns
        -------
        float
            The unnormalized probability of `a`.
        """
        return self.unprob(a)

    def py_prob(self, np.ndarray a):
        """
        Normalized probability of a point.

        Parameters
        ----------
        a : numpy.ndarray
            The point at which to evaluate the distribution.

        Returns
        -------
        float
            The normalized probability of `a`.
        """
        return self.prob(a)

    def py_dim(self):
        """
        The dimensionality of the distribution.

        Returns
        -------
        int
            The number of dimensions.
        """
        return self.dim()


cdef class UniformDistribution(Distribution):
    """
    A uniform distribution over an axis-aligned box.

    Parameters
    ----------
    lb : numpy.ndarray
        The lower bound in each dimension.
    ub : numpy.ndarray
        The upper bound in each dimension; must be the same shape as `lb`.
    """
    def __init__(self, np.ndarray lb, np.ndarray ub):
        """See class docstring."""
        self.dimension = lb.shape[0]
        self.lower_bounds = lb.copy()
        self.upper_bounds = ub.copy()
        self.volume = 1.0

        if lb.shape[0] != ub.shape[0]:
            raise SyntaxError("Wrong bound array sizes")

        cdef unsigned i
        for i in range(self.dimension):
            if self.lower_bounds[i] > self.upper_bounds[i]:
                raise SyntaxError("Lower bound > Upper bound")
            self.volume *= (self.upper_bounds[i]-self.lower_bounds[i])

    cdef double unprob(self,np.ndarray[np.double_t, ndim=1] a):
        cdef unsigned i
        for i in range(self.dimension):
            if a[i] > self.upper_bounds[i] or a[i] < self.lower_bounds[i]:
                return 0

        return 1

    cdef double prob(self,np.ndarray a):
        cdef unsigned i
        for i in range(self.dimension):
            if a[i] > self.upper_bounds[i] or a[i] < self.lower_bounds[i]:
                return 0

        return 1/self.volume


    cdef unsigned dim(self):
        return self.dimension


##################################################                ####################################################
######################################              DATA                              ################################
#################################################                     ################################################
#Top level Data Object Class.
cdef class Data():
    """
    Base class for experimental data used by a `Likelihood`.

    Holds a set of measurements of `measured_species` at `timepoints`.
    Subclasses (`BulkData`, `FlowData`, `StochasticTrajectories`) fix the
    expected array shapes for a particular kind of experimental data; see
    each subclass's `set_data` for its specific shape conventions.

    Parameters
    ----------
    timepoints : numpy.ndarray, optional
        The time points at which the data was measured.
    measurements : numpy.ndarray, optional
        The measured values.
    measured_species : list of str, optional
        The names of the measured species.
    N : int, optional
        The number of samples/trajectories (default 1).
    """
    def __init__(self, np.ndarray timepoints = None, np.ndarray measurements = None, list measured_species = [], unsigned N = 1):
        """See class docstring."""
        self.set_data(timepoints, measurements, measured_species, N)

    def set_data(self,np.ndarray timepoints, np.ndarray measurements, list measured_species, unsigned N):
        """
        Store `timepoints`, `measurements`, and `measured_species`.

        Called by `__init__`; subclasses override this to validate and
        reshape `measurements` according to their specific data layout.

        Parameters
        ----------
        timepoints : numpy.ndarray
            The time points at which the data was measured.
        measurements : numpy.ndarray
            The measured values.
        measured_species : list of str
            The names of the measured species.
        N : int
            The number of samples/trajectories.
        """
        self.timepoints = timepoints
        self.measurements = measurements
        self.measured_species = measured_species
        self.M = len(measured_species)
        self.N = N

    cdef np.ndarray get_measurements(self):
        return self.measurements

    def py_get_measurements(self):
        """
        The stored measurement array.

        Returns
        -------
        numpy.ndarray
            The measurements, in the shape set by `set_data`.
        """
        return self.get_measurements()

    def py_get_timepoints(self):
        """
        The stored time points.

        Returns
        -------
        numpy.ndarray
            The time points at which the data was measured.
        """
        return self.get_timepoints()

    def py_get_measured_species(self):
        """
        The names of the measured species.

        Returns
        -------
        list of str
            The measured species names.
        """
        return self.get_measured_species()

#Data consists of a single timecourse at T points gathered across M measurements at timempoints timepoints.
#Measurements are assumed to correspond to species names in measured_species.
#Data Dimensions:
# timepoints: T
# Measurements: T x M
# Measured Species: M
cdef class BulkData(Data):
    """
    Bulk (deterministic-style) trajectory data.

    `N` samples, each with measurements of `M` species at `T` time points.
    `measurements` is reshaped to `(N, T, M)`. `timepoints` may be a single
    length-`T` vector shared by all samples, or an `(N, T)` array giving each
    sample its own time points.
    """
    def set_data(self,np.ndarray timepoints, np.ndarray measurements, list measured_species, unsigned N):
        """
        Validate and store bulk trajectory data.

        Parameters
        ----------
        timepoints : numpy.ndarray
            A length-`T` vector (shared across samples) or, if `N > 1`, an
            `(N, T)` array of per-sample time points.
        measurements : numpy.ndarray
            The measured values; reshaped to `(N, T, M)`.
        measured_species : list of str
            The names of the measured species.
        N : int
            The number of samples/trajectories.

        Raises
        ------
        ValueError
            If `N > 1` and the shapes of `timepoints` and `measurements` are
            inconsistent with `N` samples.
        """
        if N > 1:
            if timepoints.ndim == 1:
                self.multiple_timepoints = False
                self.nT = timepoints.shape[0]
            elif timepoints.ndim == 2 and (timepoints.shape[0] == measurements.shape[0] == N):
                self.multiple_timepoints = True
                self.nT = timepoints.shape[1]
            else:
                raise ValueError("For N > 1 samples timepoints must be N x T, measurements N x T x M. Dimensiosn do not match.")
        else:
            self.multiple_timepoints = False
            self.nT = timepoints.shape[0]
        self.timepoints = timepoints
        self.N = N
        self.M = len(measured_species)
        self.measured_species = measured_species
        self.measurements = np.reshape(measurements, (self.N, self.nT, self.M))

    def has_multiple_timepoints(self):
        """
        Whether each sample has its own time points.

        Returns
        -------
        bool
            True if `timepoints` was given as an `(N, T)` array (one vector
            per sample), False if a single shared length-`T` vector was used.
        """
        return self.multiple_timepoints

    def get_minimum_dt(self):
        """
        The minimum spacing between consecutive measured time points.

        Returns
        -------
        float
            The smallest gap between consecutive time points, across all
            samples if `has_multiple_timepoints` is True.
        """
        if not self.multiple_timepoints:
            dt = self.timepoints[self.nT-1] - self.timepoints[0]
            for i in range(1, self.nT-1):
                dt = min([dt, self.timepoints[i] - self.timepoints[i-1]])
        else:
            dt = self.timepoints[0, self.nT-1] - self.timepoints[0, 0]
            for j in range(self.N):
                for i in range(1, self.nT):
                    dt = min([dt, self.timepoints[j, i] - self.timepoints[j, i-1]])

        return dt



#Data consists of a set of N flow samples which each contain measurements for M measured_species
#Data Dimensions:
# timepoints: None
# Measurements: N x M
# Measured Species: M
cdef class FlowData(Data):
    """
    Flow-cytometry-style data.

    `N` samples, each a single-timepoint measurement of `M` species. Unlike
    `BulkData` and `StochasticTrajectories`, there is no time dimension --
    `measurements` has shape `(N, M)`.
    """
    def set_data(self,np.ndarray timepoints, np.ndarray measurements, list measured_species, unsigned N):
        """
        Validate and store flow data.

        Parameters
        ----------
        timepoints : None
            Must be None -- flow data is assumed to be collected at a
            single, unspecified time point.
        measurements : numpy.ndarray
            The measured values, of shape `(N, M)`.
        measured_species : list of str
            The names of the measured species.
        N : int
            Unused; `N` is instead taken from `measurements.shape[0]`.

        Raises
        ------
        ValueError
            If `timepoints` is not None, or if the second dimension of
            `measurements` does not match `len(measured_species)`.
        """
        if timepoints is not None:
            raise ValueError("Flow Data is assumed to be collected at a single timepoint")

        self.measured_species = measured_species
        self.M = len(self.measured_species)

        if measurements.shape[1] != self.M:
            raise ValueError("Second dimension of measurments must be the same length as measured_species")
        else:
            self.measurements = measurements
            self.N = measurements.shape[0]


#Data consists of a set of N stochastic trajectories at T timepoints which each contain measurements of M measured_species.
#Data Dimensions:
# timepoints: N x T
# Measurements: N x T x M
# Measured Species: M
cdef class StochasticTrajectories(Data):
    """
    Single-cell stochastic trajectory data.

    `N` trajectories, each with measurements of `M` species at `T` time
    points. Unlike `BulkData`, `measurements` must always be given in full
    `(N, T, M)` form (no implicit `N = 1` case).
    """

    cdef  np.ndarray get_measurements(self):
        return np.reshape(self.measurements, (self.N, self.nT, self.M))

    def set_data(self, np.ndarray timepoints, np.ndarray measurements, list measured_species, unsigned N):
        """
        Validate and store stochastic trajectory data.

        Parameters
        ----------
        timepoints : numpy.ndarray
            A length-`T` vector (shared across trajectories) or an `(N, T)`
            array of per-trajectory time points.
        measurements : numpy.ndarray
            The measured values, of shape `(N, T, M)`.
        measured_species : list of str
            The names of the measured species.
        N : int
            Unused; `N` is instead taken from `measurements.shape[0]`.

        Raises
        ------
        ValueError
            If the third dimension of `measurements` does not match
            `len(measured_species)`, or if the shape of `timepoints` is
            neither a single length-`T` vector nor an `(N, T)` array.
        """
        self.measured_species = measured_species
        self.M = len(self.measured_species)

        #if (self.M == 1 and measurements.ndim == 2) or (measurements.shape[2] == self.M):
        # ( M will always be shape N X T X M (instead of having a different shape of T X M for a single trajectory case ))
        if measurements.shape[2] == self.M:
            self.nT = measurements.shape[1]
            self.N = measurements.shape[0]
            self.measurements = np.reshape(measurements, (self.N, self.nT, self.M))
        else:
            raise ValueError("Third dimension of the data array must be the same length as measured_species")


        if timepoints.shape[0] == self.nT:
            self.timepoints = timepoints
            self.multiple_timepoints = False
        elif timepoints.shape[0] == self.N and self.nT == timepoints.shape[1]:
            self.timepoints = timepoints
            self.multiple_timepoints = True
        else:
            raise ValueError("timepoints must be a single vector of length T or N (# of samples) vectors each of length T")

    def has_multiple_timepoints(self):
        """
        Whether each trajectory has its own time points.

        Returns
        -------
        bool
            True if `timepoints` was given as an `(N, T)` array (one vector
            per trajectory), False if a single shared length-`T` vector was
            used.
        """
        return self.multiple_timepoints

##################################################                ####################################################
######################################              LIKELIHOOD                        ################################
#################################################                     ################################################

cdef class Likelihood:
    """
    Base class for log-likelihood functions.

    Subclasses implement `get_log_likelihood` (a `cdef` method), called from
    Python via `py_log_likelihood` below.
    """
    cdef double get_log_likelihood(self):
        raise NotImplementedError("Calling get_log_likelihood() for Likelihood")

    def py_log_likelihood(self):
        """
        The log-likelihood of the currently-set model parameters and state.

        Returns
        -------
        float
            The log-likelihood value.
        """
        return self.get_log_likelihood()

cdef class ModelLikelihood(Likelihood):
    """
    Base class for likelihoods comparing model simulations to data.

    Compares simulations of a `~bioscrape.types.Model` against a `Data`
    object. Abstract: `set_likelihood_options` and `set_data` must be
    implemented by subclasses (`DeterministicLikelihood`,
    `StochasticTrajectoriesLikelihood`, `StochasticStatesLikelihood`).

    Parameters
    ----------
    model : Model
        The bioscrape model to simulate.
    init_state : dict or list of dict, optional
        Initial condition(s) for the simulated trajectories. A single dict
        for one trajectory, or a list of dicts for multiple trajectories.
    init_params : dict or list of dict, optional
        Per-trajectory parameter overrides, in the same shape as
        `init_state`.
    interface : CSimInterface, optional
        A pre-built simulation interface for `model`; built automatically if
        not given.
    simulator : RegularSimulator, optional
        The simulator to use; the specific subclass chooses an appropriate
        default (e.g. deterministic or stochastic) if not given.
    data : Data, optional
        The experimental data to compare simulations against. If not given,
        `set_data` must be called separately before use.
    **keywords
        Additional keyword arguments, including `hmax` (maximum integrator
        step size, deterministic simulation only) and `atol`/`rtol`
        (absolute/relative integrator error tolerances).

    Raises
    ------
    ValueError
        If `init_state` is neither a dict nor a list of dicts.
    """
    def __init__(self, model = None, init_state = {}, init_params = {},
                 interface = None, simulator = None, data = None,
                 **keywords):
        """See class docstring."""
        self.set_model(model, simulator, interface)
        self.set_likelihood_options(**keywords)
        if isinstance(init_state, dict):
            self.Nx0 = 1
            self.set_init_species([init_state])

        elif isinstance(init_state, list):
            self.Nx0 = len(init_state)
            self.set_init_species(init_state)

        else:
            raise ValueError("Init_state must either be a dictionary or a list of dictionaries.")

        if isinstance(init_params, dict):
            self.Nx0 = 1
            if init_params:
                self.initial_parameters = [init_params]
        elif isinstance(init_params, list):
            self.Nx0 = len(init_params)
            self.initial_parameters = init_params
        else:
            self.initial_parameters = None
        
        if data is not None:
            self.set_data(data)
        else:
            print("Warning: No Data passed into likelihood constructor. Must call set_data seperately.")

        #hmax is used in deterministic simulation to set the maximum step size of the integrator
        if 'hmax' in keywords and hasattr(self, 'py_set_hmax'):
            self.py_set_hmax(keywords['hmax'])

        #atol and rtol are the absolute and relative error tollerances used by the integrator
        if ('rtol' in keywords or 'atol' in keywords) and hasattr(self.propagator, "py_set_tolerance"):
            if 'atol' not in keywords:
                print(f"{self} only recieved 'rtol' keyword. Setting 'atol' to the 'rtol' value.")
                atol = keywords['rtol']
                rtol = keywords['rtol']
            elif 'rtol' not in keywords:
                print(f"{self} only recieved 'atol' keyword. Setting 'rtol' to the 'atol' value.")
                atol = keywords['atol']
                rtol = keywords['atol']
            else:
                atol = keywords['atol']
                rtol = keywords['rtol']

            self.propagator.py_set_tolerance(atol, rtol)

    def set_model(self, Model m, RegularSimulator prop, CSimInterface csim = None, DelaySimulator prop_delay = None):
        """
        Attach a model, simulator, and simulation interface.

        Parameters
        ----------
        m : Model
            The bioscrape model to simulate.
        prop : RegularSimulator
            The simulator to use.
        csim : CSimInterface, optional
            A pre-built simulation interface for `m`; built automatically
            with `~bioscrape.simulator.ModelCSimInterface` if not given.
        prop_delay : DelaySimulator, optional
            A delay-aware simulator, for subclasses that support delayed
            reactions.
        """
        if csim is None:
            self.csim = ModelCSimInterface(m)
        self.m = m
        self.propagator = prop
        self.propagator_delay = prop_delay
        self.set_default_params(dict(self.m.get_parameter_dictionary()))
        self.set_default_species(self.m.get_species_array())


    def set_default_species(self, species):
        """
        Set the default species values.

        Used to fill in unspecified initial conditions.

        Parameters
        ----------
        species : numpy.ndarray
            The default species value array, indexed as in the model
            attached via `set_model`.
        """
        self.default_species = species

    def set_default_params(self, params):
        """
        Set (or update) the default parameter values.

        Used to fill in unspecified per-trajectory parameter conditions.

        Parameters
        ----------
        params : dict
            Parameter name/value pairs. Merged into any existing defaults
            rather than replacing them.
        """
        if hasattr(self, "default_params"):
            self.default_params.update(params)
        else:
            self.default_params = dict(params)

    def set_init_species(self, list sds):
        """
        Set the initial species values for each trajectory.

        Parameters
        ----------
        sds : list of dict
            One dictionary of initial species values per trajectory (length
            `self.Nx0`); species not given a value fall back to the default
            set by `set_default_species`.
        """
        self.initial_states = np.zeros((self.Nx0, self.m.get_number_of_species()))
        species2index = self.m.get_species2index()
        
        for i in range(self.Nx0):
            for s in species2index:
                j  = species2index[s]
                if s in sds[i]:
                    self.initial_states[i, j] = sds[i][s]
                else:
                    self.initial_states[i, j] = self.default_species[j]

    def set_init_params(self, dict pd):
        """
        Set model parameter values for the current trajectory.

        Parameters
        ----------
        pd : dict
            Parameter name/value pairs to set on the model and simulation
            interface.
        """
        # print("before set", self.m.get_parameter_dictionary())
        # print("asking to be set to", pd)
        self.m.set_params(pd)
        self.csim.py_set_param_values(self.m.get_params_values())
        # print("After set", self.m.get_parameter_dictionary())
        #pairs = pd.items()
        #self.init_param_indices = np.zeros(len(pairs), dtype=int)
        #self.init_param_vals = np.zeros(len(pairs),)

        #index = 0
        #for (key,val) in pairs:
        #    self.init_param_indices[index] = self.m.get_param_index( key  )
        #    self.init_param_vals[index] = val
        #    index +=  1

    cdef np.ndarray get_initial_state(self, int n):
        # Return the initial state for the nth trajectory 
        # (this may correspond to multiple measurements m)
        return self.initial_states[n, :]

    cdef dict get_initial_params(self, int n):
        # Return the initial params dict for the nth trajectory 
        # (this may correspond to multiple measurements m)
        if self.initial_parameters is None:
            return None
        return self.initial_parameters[n]

    def set_likelihood_options(self, **keywords):
        """
        Set likelihood-specific options.

        Must be implemented by subclasses; see e.g.
        `DeterministicLikelihood.set_likelihood_options`.

        Parameters
        ----------
        **keywords
            Subclass-specific option names and values.

        Raises
        ------
        NotImplementedError
            Always, unless overridden by a subclass.
        """
        raise NotImplementedError("set_likelihood_options must be implemented in subclasses of ModelLikelihood")

    def set_data(self, **keywords):
        """
        Attach experimental data to compare simulations against.

        Must be implemented by subclasses; see e.g.
        `DeterministicLikelihood.set_data`.

        Parameters
        ----------
        **keywords
            Subclass-specific arguments (typically a single `Data` object).

        Raises
        ------
        NotImplementedError
            Always, unless overridden by a subclass.
        """
        raise NotImplementedError("set_data must be implemented in subclasses of ModelLikelihood")

cdef class DeterministicLikelihood(ModelLikelihood):
    """
    Log-likelihood of `BulkData` under deterministic simulation.

    Compares deterministically-simulated trajectories to `BulkData` using a
    `norm_order`-norm distance; see `set_likelihood_options`.
    """
    def set_model(self, Model m, RegularSimulator prop = None, CSimInterface csim = None):
        """
        Attach a model, defaulting to a deterministic simulator.

        Parameters
        ----------
        m : Model
            The bioscrape model to simulate.
        prop : RegularSimulator, optional
            The simulator to use (default a new
            `~bioscrape.simulator.DeterministicSimulator`).
        csim : CSimInterface, optional
            A pre-built simulation interface for `m`; built automatically if
            not given.
        """
        if prop is None:
            prop = DeterministicSimulator()
        ModelLikelihood.set_model(self, m, prop, csim)
        self.csim.py_prep_deterministic_simulation()

        #self.m = m
        #if csim is None:
        #    csim = ModelCSimInterface(m)

        #self.csim = csim

        #self.propagator = prop


    def set_data(self, BulkData bd):
        """
        Attach bulk trajectory data to compare simulations against.

        Parameters
        ----------
        bd : BulkData
            The experimental data.

        Raises
        ------
        ValueError
            If the number of samples in `bd` and the number of initial
            conditions set via `set_init_species` are incompatible (neither
            equal nor either equal to 1).
        """
        #Set bulk data
        self.bd = bd
        self.N = self.bd.get_N() #Number of samples
        self.propagator.py_set_hmax(self.bd.get_minimum_dt()) 

        #Get species indices in Model
        species_list = bd.get_measured_species()
        self.M = len(species_list) #Number of measured species
        self.meas_indices = np.zeros(len(species_list), dtype=int)
        for i in range(self.M):
            self.meas_indices[i] = self.m.get_species_index(species_list[i])
        

        if not ((self.N == 1 or self.Nx0 == 1) or (self.Nx0 == self.N)):
            print("self.N", self.N, "self.Nx0", self.Nx0)
            raise ValueError("Either the number of samples and the number of"
                             "initial conditions match or one of them must be 1")

    def py_get_data(self):
        """
        The attached data.

        Returns
        -------
        BulkData
            The data set by `set_data`.
        """
        return self.bd

    def set_likelihood_options(self, norm_order = 1, **keywords):
        """
        Set the norm order used to compare trajectories.

        Parameters
        ----------
        norm_order : int, optional
            The order of the norm used in the log-likelihood distance
            (default 1).
        """
        self.norm_order = norm_order

    def py_get_norm_order(self):
        """
        The norm order used to compare trajectories.

        Returns
        -------
        int
            The norm order set by `set_likelihood_options`.
        """
        return self.norm_order

    def py_set_hmax(self, hmax):
        """
        Set the integrator's maximum step size.

        Parameters
        ----------
        hmax : float
            The maximum step size for the deterministic integrator.
        """
        self.hmax = hmax
        self.propagator.py_set_hmax(hmax)

    

    cdef double get_log_likelihood(self):
        # Write in the specific parameters and species values.
        cdef np.ndarray[np.double_t, ndim = 1] species_vals = self.m.get_species_values()
        cdef np.ndarray[np.double_t, ndim = 1] param_vals = self.m.get_params_values()
        cdef np.ndarray[np.double_t, ndim = 2] ans
        cdef np.ndarray[np.double_t, ndim = 1] timepoints
        cdef unsigned i, t
        cdef double error = 0.0
        cdef double dif = 0
        cdef np.ndarray[np.double_t, ndim = 3] measurements = self.bd.get_measurements()
        #cdef SSAResult res

        #Go through trajectories (which may have unique initial states)
        for n in range(self.N):
            #Set Timepoints
            if self.bd.has_multiple_timepoints():
                timepoints = self.bd.get_timepoints()[n, :]
            else:
                timepoints = self.bd.get_timepoints()

            #Reset csim interface - SOLVES ISSUES WITH RULES
            #self.csim = ModelCSimInterface(self.m)
            #self.csim.py_prep_deterministic_simulation()

            self.csim.set_initial_state(self.get_initial_state(n))
            # print('current params conditions', self.csim.py_get_param_values())
            # print('lets set for new traj to', self.get_initial_params(n))
            if self.get_initial_params(n) is not None:
                self.set_init_params(self.get_initial_params(n))
            ans = self.propagator.simulate(self.csim, timepoints).get_result()

            #print("simulation", t12-t11)
            # Compare the data using norm and return the likelihood.
            #Cycle through all the measurements of each trajectory
            for i in range(self.M):
                for t in range(len(timepoints)):
                    dif = measurements[n, t, i] - ans[t,self.meas_indices[i]]
                    if dif < 0:
                        dif = -dif
                    error += dif**self.norm_order

        error = error**(1./self.norm_order)

        if np.isnan(error):
            return -np.inf
        else:
            return -error

cdef class StochasticTrajectoriesLikelihood(ModelLikelihood):
    """
    Log-likelihood of `StochasticTrajectories` data.

    Compares the average of `N_simulations` stochastic simulations per
    trajectory to `StochasticTrajectories` data using a `norm_order`-norm
    distance; see `set_likelihood_options`.
    """
    def set_model(self, Model m, RegularSimulator prop = None, CSimInterface csim = None, DelaySimulator prop_delay = None):
        """
        Attach a model, defaulting to an SSA-based simulator.

        Defaults to `~bioscrape.simulator.SSASimulator`, or a
        `~bioscrape.simulator.DelaySSASimulator` if `m` has delayed
        reactions.

        Parameters
        ----------
        m : Model
            The bioscrape model to simulate.
        prop : RegularSimulator, optional
            The simulator to use; chosen automatically based on
            `m.has_delay` if not given.
        csim : CSimInterface, optional
            A pre-built simulation interface for `m`; built automatically if
            not given.
        prop_delay : DelaySimulator, optional
            The delay-aware simulator to use if `m` has delayed reactions.
        """
        if prop is None:
            if m.has_delay:
                if prop_delay is None:
                    prop_delay = DelaySSASimulator()
                self.has_delay = True
            else:
                prop = SSASimulator()
                self.has_delay = False
        ModelLikelihood.set_model(self, m, prop, csim, prop_delay)

        #self.m = m
        #self.has_delay = False
        #if csim is None:
        #    csim = ModelCSimInterface(m)
        
        #else:
        #    self.propagator = prop
        #self.csim = csim

    def set_data(self, StochasticTrajectories sd):
        """
        Attach stochastic trajectory data to compare simulations against.

        Parameters
        ----------
        sd : StochasticTrajectories
            The experimental data.

        Raises
        ------
        ValueError
            If the number of samples in `sd` and the number of initial
            conditions set via `set_init_species` are incompatible (neither
            equal nor either equal to 1).
        """
        self.sd = sd #Holds Data Object
        self.N = sd.get_N() #Number of samples

        #Get Species Model Indices
        species_list = sd.get_measured_species()
        self.M = len(species_list) #Number of
        self.meas_indices = np.zeros(len(species_list), dtype=int)

        self.M = len(species_list)
        for i in range(len(species_list)):
            self.meas_indices[i] = self.m.get_species_index(species_list[i])

        #The number of intitial conditions in the likelihood model must either:
        # 1: Same initial condition used for all samples
        # N: Unique initial condition used for each sample
        if not ((self.N == 1 or self.Nx0 == 1) or (self.Nx0 == self.N)):
            raise ValueError("Either the number of samples and the number of"
                             "initial conditions match or one of them must be 1")
    def py_get_data(self):
        """
        The attached data.

        Returns
        -------
        StochasticTrajectories
            The data set by `set_data`.
        """
        return self.sd

    def set_likelihood_options(self, N_simulations = None, norm_order = None, **keywords):
        """
        Set the number of repeated simulations and the norm order.

        Parameters
        ----------
        N_simulations : int, optional
            The number of stochastic simulations to average per trajectory
            (default 1).
        norm_order : int, optional
            The order of the norm used in the log-likelihood distance
            (default 1).
        """
        if norm_order is not None:
            self.norm_order = norm_order
        if norm_order in [None, 0]:
            self.norm_order = 1

        if N_simulations is not None:
            self.N_simulations = N_simulations
        if self.N_simulations in [None, 0]:
            self.N_simulations = 1

    def py_get_N_simulations(self):
        """
        The number of repeated simulations averaged per trajectory.

        Returns
        -------
        int
            The value set by `set_likelihood_options`.
        """
        return self.N_simulations
    def py_get_norm_order(self):
        """
        The norm order used to compare trajectories.

        Returns
        -------
        int
            The value set by `set_likelihood_options`.
        """
        return self.norm_order

    cdef double get_log_likelihood(self):
        # Write in the specific parameters and species values.
        cdef np.ndarray[np.double_t, ndim = 1] species_vals = self.m.get_species_values()
        cdef np.ndarray[np.double_t, ndim = 1] param_vals = self.m.get_params_values()
        cdef np.ndarray[np.double_t, ndim = 1] timepoints

        cdef unsigned i
        cdef unsigned n
        cdef unsigned s
        cdef double error = 0.0
        cdef double dif = 0
        cdef np.ndarray[np.double_t, ndim = 2] ans
        cdef np.ndarray[np.double_t, ndim = 3] measurements = self.sd.get_measurements()

        #Go through trajectories (which may have unique initial states)
        for n in range(self.N):
            #Set Timepoints
            if self.sd.has_multiple_timepoints():
                timepoints = self.sd.get_timepoints()[n, :]
            else:
                timepoints = self.sd.get_timepoints()
            # Set initial state conditions
            self.csim.set_initial_state(self.get_initial_state(n))
            # Set init params
            if self.get_initial_params(n) is not None:
                self.set_init_params(self.get_initial_params(n))
            # Do N*N_simulations simulations of the model with time points specified by the data.
            for s in range(self.N_simulations):
                # Initial parameters are set for each trajectory passed in.
                # Initial conditions are set for each trajectory passsed in.
                if self.has_delay:
                    q = ArrayDelayQueue.setup_queue(self.csim.py_get_num_reactions(), len(timepoints),timepoints[1]-timepoints[0])
                    ans = self.propagator_delay.delay_simulate(self.csim, q, timepoints).get_result()
                elif not self.has_delay:
                    ans = self.propagator.simulate(self.csim, timepoints).get_result()
                for i in range(self.M):
                    # Compare the data using norm and return the likelihood.
                    for t in range(len(timepoints)):
                        dif = measurements[n, t, i] - ans[t,self.meas_indices[i]]
                        if dif < 0:
                            dif = -dif
                        error += dif**self.norm_order
        error = error**(1./self.norm_order)
        error = -1.0*error/(1.0*self.N_simulations)

        if np.isnan(error):
            return -np.inf
        else:
            return error

cdef class StochasticTrajectoryMomentLikelihood(StochasticTrajectoriesLikelihood):
    """
    `StochasticTrajectoriesLikelihood` with moment-matching options.

    Adds `initial_state_matching` and `Moments` options to
    `set_likelihood_options`; the log-likelihood computation itself is
    inherited unchanged from `StochasticTrajectoriesLikelihood`.
    """
    def set_likelihood_options(self, N_simulations = 1, initial_state_matching = False, Moments = 2, **keywords):
        """
        Set the moment-matching options for this likelihood.

        Parameters
        ----------
        N_simulations : int, optional
            Present for interface compatibility with
            `StochasticTrajectoriesLikelihood.set_likelihood_options`;
            always set to 1 regardless of the value passed.
        initial_state_matching : bool, optional
            Stored for use by subclasses (default False).
        Moments : int, optional
            The number of moments to match: 1 (averages) or 2 (averages and
            second moments); default 2.

        Raises
        ------
        ValueError
            If `Moments` is greater than 2.
        """
        self.N_simulations = 1
        self.initial_state_matching = initial_state_matching
        if Moments > 2:
            raise ValueError("Moments must be 1 (Averages) or 2 (Averages and 2nd Moments)")
        self.Moments = Moments

cdef class StochasticStatesLikelihood(ModelLikelihood):
    """
    Log-likelihood of `FlowData` under stochastic simulation.

    Compares single-timepoint state distributions from stochastic
    simulation against `FlowData`.
    """
    def set_model(self, Model m, RegularSimulator prop = None, CSimInterface csim = None, DelaySimulator prop_delay = None):
        """
        Attach a model, defaulting to an SSA-based simulator.

        Defaults to `~bioscrape.simulator.SSASimulator`, or a
        `~bioscrape.simulator.DelaySSASimulator` if `m` has delayed
        reactions.

        Parameters
        ----------
        m : Model
            The bioscrape model to simulate.
        prop : RegularSimulator, optional
            The simulator to use; chosen automatically based on
            `m.has_delay` if not given.
        csim : CSimInterface, optional
            A pre-built simulation interface for `m`; built automatically if
            not given.
        prop_delay : DelaySimulator, optional
            The delay-aware simulator to use if `m` has delayed reactions.
        """
        self.m = m

        if csim is None:
            csim = ModelCSimInterface(m)
        if prop is None:
            if m.has_delay:
                if prop_delay is None:
                    prop_delay = DelaySSASimulator()
                self.has_delay = True
                self.propagator_delay = prop_delay
            else:
                prop = SSASimulator()
                self.propagator = prop
        else:
            self.propagator = prop

        self.csim = csim

    def set_data(self, FlowData fd):
        """
        Attach flow data to compare simulations against.

        Parameters
        ----------
        fd : FlowData
            The experimental data.
        """
        self.fd = fd
        the_list = fd.get_measured_species()
        self.meas_indices = np.zeros(len(the_list), dtype=int)
        for i in range(len(the_list)):
            self.meas_indices[i] = self.m.get_species_index(the_list[i])


def py_inference(Model = None, params_to_estimate = None, exp_data = None, initial_conditions = None,
                 parameter_conditions = None, measurements = None, time_column = None, nwalkers = None, 
                 nsteps = None, init_seed = None, prior = None, sim_type = None, inference_type = 'emcee',
                 method = 'mcmc', plot_show = True, parallel = None, **kwargs):
    """
    User-level interface for Bayesian parameter inference.

    Fits `params_to_estimate` to `exp_data` via Markov Chain Monte
    Carlo (MCMC) sampling, using
    `emcee <https://emcee.readthedocs.io/>`_ under the hood. Returns
    the `emcee` sampler (holding the full set of posterior samples)
    together with the `~bioscrape.inference_setup.InferenceSetup`
    object that orchestrated the run; unless `plot_show=False`, it
    also plots the resulting posterior parameter distributions, and
    always writes the raw samples and a fit summary to
    `filename_csv`/`filename_txt`.

    Parameters
    ----------
    Model : Model
        The bioscrape model to fit.
    params_to_estimate : list of str
        The names of the parameters in `Model` to estimate.
    exp_data : pandas.DataFrame or list of pandas.DataFrame
        The experimental data to fit against.
    initial_conditions : dict or list of dict, optional
        Initial condition(s) for the simulated trajectories, in the same
        shape as `exp_data`. Defaults to the model's own initial conditions.
    parameter_conditions : dict or list of dict, optional
        Per-trajectory parameter overrides, in the same shape as `exp_data`.
    measurements : list of str or str, optional
        The name(s) of the species in `Model` that are measured.
    time_column : str, optional
        The column name in `exp_data` that holds the time points, for
        time-series inference.
    nwalkers : int, optional
        The number of walkers for MCMC sampling; see
        `emcee.EnsembleSampler`.
    nsteps : int, optional
        The number of steps for MCMC sampling; see `emcee.EnsembleSampler`.
    init_seed : float or numpy.ndarray or list or 'prior', optional
        Controls the initial parameter values for the sampler. A float `r`
        samples initial values from a Gaussian ball of radius `r` around
        the model's parameter values; an array or list of length
        `len(params_to_estimate)` is used directly as the initial values;
        'prior' samples the initial values from `prior`.
    prior : dict, optional
        The prior for each parameter in `params_to_estimate`. For example::

            {'parameter1': ['uniform', min_value, max_value],
             'parameter2': ['gaussian', mean, std, 'positive']}

        The 'positive' flag rejects negative parameter values. See the
        `bioscrape wiki <https://github.com/biocircuits/bioscrape/wiki>`_
        for the full list of supported priors.
    sim_type : {'deterministic', 'stochastic'}, optional
        The type of simulation to use for the inference run.
    inference_type : {'emcee', 'lmfit'}, optional
        The inference package to use (default 'emcee').
    method : str, optional
        Ignored for `inference_type='emcee'`. For `inference_type='lmfit'`,
        passed as the `method` keyword to `lmfit.minimize`; see the
        `lmfit docs <https://lmfit.github.io/lmfit-py/fitting.html>`_.
    plot_show : bool, optional
        If True, display the plots generated by the inference run (default
        True).
    parallel : bool, optional
        If True, run the sampler with a `multiprocessing.Pool` passed to
        `emcee.EnsembleSampler` for parallel processing. If False (default),
        multiprocessing is not used.
    **kwargs
        Additional keyword arguments passed into the inference setup.

    Returns
    -------
    tuple or lmfit.minimizer.MinimizerResult
        For `inference_type='emcee'`: a tuple `(sampler, pid)` of the
        `emcee.EnsembleSampler` and the bioscrape PID object. For
        `inference_type='lmfit'`: an `lmfit.minimizer.MinimizerResult`.

    Raises
    ------
    ValueError
        If `Model` is None, or if `inference_type` is not one of the
        supported values ('emcee', 'lmfit').
    """
    if Model is None:
        raise ValueError('Model object cannot be None.')
        
    pid = initialize_inference(Model = Model, **kwargs)
    if exp_data is not None:
        pid.set_exp_data(exp_data)
    if measurements is not None:
        pid.set_measurements(measurements)
    if initial_conditions is None:
        initial_conditions = dict(Model.get_species_dictionary())
    pid.set_initial_conditions(initial_conditions)
    pid.set_parameter_conditions(parameter_conditions)
    if time_column is not None:
        pid.set_time_column(time_column)
    if nwalkers is not None:
        pid.set_nwalkers(nwalkers)
    if init_seed is not None:
        pid.set_init_seed(init_seed)
    if nsteps is not None:
        pid.set_nsteps(nsteps)
    if sim_type is not None:
        pid.set_sim_type(sim_type)
    if parallel is not None:
        pid.set_parallel(parallel)
    if params_to_estimate is not None:
        pid.set_params_to_estimate(params_to_estimate)
    if prior is not None:
        pid.set_prior(prior)
    if inference_type == 'emcee' and method == 'mcmc':
        sampler = pid.run_mcmc(**kwargs)
        if plot_show:
            pid.plot_mcmc_results(sampler, **kwargs)
        return sampler, pid
    elif inference_type == 'lmfit':
        minimizer_result = pid.run_lmfit(method = method,
                                         plot_show = plot_show, **kwargs)
        pid.write_lmfit_results(minimizer_result, **kwargs)
        return minimizer_result
    else:
        raise ValueError("Set inference_type keyword argument to your"
                         "preferred inference package name. Currently"
                         "emcee and lmfit are supported.")
