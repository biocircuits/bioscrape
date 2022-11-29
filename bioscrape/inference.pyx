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
    cdef double unprob(self,np.ndarray a):
        raise NotImplementedError("Calling unprob() on Distribution")
    cdef double prob(self,np.ndarray a):
        raise NotImplementedError("Calling prob() on Distribution")
    cdef unsigned dim(self):
        raise NotImplementedError("Calling dim() on Distribution")

    def py_unprob(self, np.ndarray a):
        return self.unprob(a)

    def py_prob(self, np.ndarray a):
        return self.prob(a)

    def py_dim(self):
        return self.dim()


cdef class UniformDistribution(Distribution):
    def __init__(self, np.ndarray lb, np.ndarray ub):
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
    def __init__(self, np.ndarray timepoints = None, np.ndarray measurements = None, list measured_species = [], unsigned N = 1):
        self.set_data(timepoints, measurements, measured_species, N)

    def set_data(self,np.ndarray timepoints, np.ndarray measurements, list measured_species, unsigned N):
        self.timepoints = timepoints
        self.measurements = measurements
        self.measured_species = measured_species
        self.M = len(measured_species)
        self.N = N

    cdef np.ndarray get_measurements(self):
        return self.measurements

    def py_get_measurements(self):
        return self.get_measurements()

    def py_get_timepoints(self):
        return self.get_timepoints()

    def py_get_measured_species(self):
        return self.get_measured_species

#Data consists of a single timecourse at T points gathered across M measurements at timempoints timepoints. 
#Measurements are assumed to correspond to species names in measured_species.
#Data Dimensions:
# timepoints: T
# Measurements: T x M
# Measured Species: M
cdef class BulkData(Data):
    def set_data(self,np.ndarray timepoints, np.ndarray measurements, list measured_species, unsigned N):
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
        return self.multiple_timepoints

    def get_minimum_dt(self):
        """
        Returns the minimum delta between measured timepoints
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
    def set_data(self,np.ndarray timepoints, np.ndarray measurements, list measured_species, unsigned N):
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

    cdef  np.ndarray get_measurements(self):
        return np.reshape(self.measurements, (self.N, self.nT, self.M))

    def set_data(self, np.ndarray timepoints, np.ndarray measurements, list measured_species, unsigned N):
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
        return self.multiple_timepoints

##################################################                ####################################################
######################################              LIKELIHOOD                        ################################
#################################################                     ################################################

cdef class Likelihood:
    cdef double get_log_likelihood(self):
        raise NotImplementedError("Calling get_log_likelihood() for Likelihood")

    def py_log_likelihood(self):
        return self.get_log_likelihood()

cdef class ModelLikelihood(Likelihood):
    def __init__(self, model = None, init_state = {}, init_params = {},
                 interface = None, simulator = None, data = None,
                 **keywords):
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
        if csim is None:
            self.csim = ModelCSimInterface(m)
        self.m = m
        self.propagator = prop
        self.propagator_delay = prop_delay
        self.set_default_params(dict(self.m.get_parameter_dictionary()))
        self.set_default_species(self.m.get_species_array())


    def set_default_species(self, species):
        self.default_species = species
       
    def set_default_params(self, params):
        if hasattr(self, "default_params"):
            self.default_params.update(params)
        else:
            self.default_params = dict(params)

    def set_init_species(self, list sds):
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
        raise NotImplementedError("set_likelihood_options must be implemented in subclasses of ModelLikelihood")

    def set_data(self, **keywords):
        raise NotImplementedError("set_data must be implemented in subclasses of ModelLikelihood")

cdef class DeterministicLikelihood(ModelLikelihood):
    def set_model(self, Model m, RegularSimulator prop = None, CSimInterface csim = None):
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
        return self.bd

    def set_likelihood_options(self, norm_order = 1, **keywords):
        self.norm_order = norm_order

    def py_get_norm_order(self):
        return self.norm_order

    def py_set_hmax(self, hmax):
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
    def set_model(self, Model m, RegularSimulator prop = None, CSimInterface csim = None, DelaySimulator prop_delay = None):
        if prop is None:
            if m.has_delay:
                prop = DelaySSASimulator()
                self.has_delay = True
                self.propagator_delay = prop_delay
                prop = None
            else:
                prop = SSASimulator()
        ModelLikelihood.set_model(self, m, prop, csim)

        #self.m = m
        #self.has_delay = False
        #if csim is None:
        #    csim = ModelCSimInterface(m)
        
        #else:
        #    self.propagator = prop
        #self.csim = csim

    def set_data(self, StochasticTrajectories sd):

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
        return self.sd

    def set_likelihood_options(self, N_simulations = None, norm_order = None, **keywords):
        if norm_order is not None:
            self.norm_order = norm_order
        if norm_order in [None, 0]:
            self.norm_order = 1

        if N_simulations is not None:
            self.N_simulations = N_simulations
        if self.N_simulations in [None, 0]:
            self.N_simulations = 1

    def py_get_N_simulations(self):
        return self.N_simulations
    def py_get_norm_order(self):
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
    def set_likelihood_options(self, N_simulations = 1, initial_state_matching = False, Moments = 2, **keywords):
        self.N_simulations = 1
        self.initial_state_matching = initial_state_matching
        if Moments > 2:
            raise ValueError("Moments must be 1 (Averages) or 2 (Averages and 2nd Moments)")
        self.Moments = Moments

cdef class StochasticStatesLikelihood(ModelLikelihood):
    def set_model(self, Model m, RegularSimulator prop = None, CSimInterface csim = None, DelaySimulator prop_delay = None):
        self.m = m

        if csim is None:
            csim = ModelCSimInterface(m)
        if prop is None:
            if m.has_delay:
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
        self.fd = fd
        the_list = fd.get_measured_species()
        self.meas_indices = np.zeros(len(the_list), dtype=int)
        for i in range(len(the_list)):
            self.meas_indices[i] = self.m.get_species_index(the_list[i])


def py_inference(Model = None, params_to_estimate = None, exp_data = None, initial_conditions = None,
                 parameter_conditions = None, measurements = None, time_column = None, nwalkers = None, 
                 nsteps = None, init_seed = None, prior = None, sim_type = None, inference_type = 'emcee',
                 method = 'mcmc', plot_show = True, **kwargs):
    
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
    if params_to_estimate is not None:
        pid.set_params_to_estimate(params_to_estimate)
    if prior is not None:
        pid.set_prior(prior)
    if inference_type == 'emcee' and method == 'mcmc':
        sampler = pid.run_mcmc(plot_show = plot_show, **kwargs)
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
