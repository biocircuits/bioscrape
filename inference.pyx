cimport numpy as np
import numpy as np
from libc.math cimport log

from types import Model
from types cimport Model
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

cdef class BulkData:
    def __init__(self):
        self.measured_species = []
        self.timepoints = None
        self.measurements = None

    def set_data(self,np.ndarray timepoints, np.ndarray measurements, list measured_species):
        self.timepoints = timepoints
        self.measurements = measurements
        self.measured_species = measured_species



##################################################                ####################################################
######################################              LIKELIHOOD                        ################################
#################################################                     ################################################

cdef class Likelihood:
    cdef double get_log_likelihood(self):
        raise NotImplementedError("Calling get_log_likelihood() for Likelihood")

cdef class DeterministicLikelihood(Likelihood):
    def __init__(self):
        self.m = None
        self.init_state_vals = np.zeros(0)
        self.init_param_vals = np.zeros(0)
        self.init_param_indices = np.zeros(0,dtype=int)
        self.init_state_indices = np.zeros(0,dtype=int)
        self.csim = None
        self.propagator = None


    def set_model(self, Model m, CSimInterface csim = None, RegularSimulator prop = None):
        self.m = m
        self.csim = csim
        self.propagator = prop

        if csim is None:
            csim = ModelCSimInterface(m)
        if prop is None:
            prop = DeterministicSimulator()


    def set_data(self, BulkData bd):
        self.bd = bd
        the_list = bd.get_measured_species()
        self.meas_indices = np.zeros(len(the_list), dtype=int)
        for i in range(len(the_list)):
            self.meas_indices[i] = self.m.get_species_index( the_list[i] )


    def set_init_species(self, dict sd):
        pairs = sd.items()
        self.init_state_indices = np.zeros(len(pairs), dtype=int)
        self.init_state_vals = np.zeros(len(pairs),)

        index = 0
        for (key,val) in  pairs:
            self.init_state_indices[index] = self.m.get_species_index( key  )
            self.init_state_vals[index] = val
            index +=  1


    def set_init_params(self, dict pd):
        pairs = pd.items()
        self.init_param_indices = np.zeros(len(pairs), dtype=int)
        self.init_param_vals = np.zeros(len(pairs),)

        index = 0
        for (key,val) in pairs:
            self.init_param_indices[index] = self.m.get_param_index( key  )
            self.init_param_vals[index] = val
            index +=  1


    cdef double get_log_likelihood(self):
        # Write in the specific parameters and species values.
        cdef np.ndarray species_vals = self.m.get_species_values()
        cdef np.ndarray param_vals = self.m.get_params_values()

        cdef unsigned i
        for i in range(self.init_state_indices.shape[0]):
            species_vals[ self.init_state_indices[i] ] = self.init_state_vals[i]
        for i in range(self.init_param_indices.shape[0]):
            param_vals[ self.init_param_indices[i] ] = self.init_param_vals[i]

        # Do a simulation of the model with time points specified by the data.
        cdef np.ndarray ans = self.propagator.simulate(self.csim, self.bd.get_timepoints()).get_result()

        # Compare the data using l_1 norm and return the likelihood.
        cdef double error = 0.0
        for i in range(self.meas_indices.shape[0]):
            error += np.linalg.norm(  self.bd.get_measurements()[:,i] - ans[:,self.meas_indices[i]] )

        return -error

##################################################                ####################################################
######################################              INFERENCE                         ################################
#################################################                     ################################################


cdef class DeterministicInference:
    def __init__(self):
        self.m = None
        self.params_to_estimate = np.zeros(0,dtype=int)
        self.global_init_params = np.zeros(0,)
        self.global_init_state  = np.zeros(0,)
        self.prior = None
        self.num_walkers = 500
        self.num_iterations = 100
        self.likelihoods = []
        self.dimension = 0
        self.sigma = 10

    def set_model(self, Model m):
        self.m = m
        self.global_init_params = m.get_params_values().copy()
        self.global_init_state = m.get_species_values().copy()

    def set_sigma(self, double s):
        self.sigma = s

    def set_mcmc_params(self, unsigned walkers, unsigned iterations):
        self.num_walkers = walkers
        self.num_iterations = iterations

    def set_prior(self, Distribution prior):
        self.prior = prior

    def set_params_to_estimate(self, list params):
        self.params_to_estimate = np.zeros(len(params),dtype=int)
        cdef unsigned i
        for i in range(len(params)):
            self.params_to_estimate[i] =  self.m.get_param_index(params[i])

        self.dimension = len(params)

    def add_to_likelihood(self,list likelihoods):
        self.likelihoods.extend(likelihoods)

    def likelihood_function(self, np.ndarray params):
        cdef unsigned i
        cdef unsigned j
        cdef double answer = self.prior.unprob(params)

        if answer == 0.0:
            return -np.inf

        answer = log(answer)

        cdef np.ndarray actual_species = self.m.get_species_values()
        cdef np.ndarray actual_params = self.m.get_params_values()

        for i in range(len(self.likelihoods)):
            # Copy in the global initial conditions.
            np.copyto(actual_species, self.global_init_state)
            np.copyto(actual_params, self.global_init_params)

            # Then copy in the specified parameters
            for j in range(self.params_to_estimate.shape[0]):
                actual_params[ self.params_to_estimate[j] ]  = params[j]

            # The run the likelihood
            answer += (<Likelihood> self.likelihoods[i]).get_log_likelihood()

        return answer / self.sigma**2

    def py_likelihood_function(self, np.ndarray params):
        return self.likelihood_function(params)

    def run_mcmc(self, p0 = None):
        def lnprob(np.ndarray theta):
            loglhood = self.likelihood_function(np.exp(theta))
            if np.isnan(loglhood):
                sys.stderr.write('Parameters returned NaN likelihood: ' + str(np.exp(theta)) + '\n')
                sys.stderr.flush()
                return -np.Inf
            return loglhood

        sampler = emcee.EnsembleSampler(self.num_walkers, self.dimension, lnprob)
        if p0 is None:
            p0 = np.random.randn(self.dimension*self.num_walkers).reshape((self.num_walkers,self.dimension)) / 20.0
        sampler.run_mcmc(p0, self.num_iterations)

        return sampler



