cimport numpy as np
from libcpp cimport bool
from types import Model
from types cimport Model
from simulator cimport CSimInterface, RegularSimulator, ModelCSimInterface, DeterministicSimulator, SSASimulator
from simulator import CSimInterface, RegularSimulator, ModelCSimInterface, DeterministicSimulator, SSASimulator
from simulator cimport DelaySimulator, DelaySSASimulator, ArrayDelayQueue
from simulator import DelaySimulator, DelaySSASimulator, ArrayDelayQueue


##################################################                ####################################################
######################################              DISTRIBUTION                      ################################
#################################################                     ################################################


cdef class Distribution:
    cdef double unprob(self,np.ndarray a)
    cdef double prob(self,np.ndarray a)
    cdef unsigned dim(self)


cdef class UniformDistribution(Distribution):
    cdef np.ndarray lower_bounds
    cdef np.ndarray upper_bounds
    cdef unsigned dimension
    cdef double volume

    cdef double unprob(self,np.ndarray a)
    cdef double prob(self,np.ndarray a)
    cdef unsigned dim(self)


##################################################                ####################################################
######################################              DATA                              ################################
#################################################                     ################################################
cdef class Data():
    cdef list measured_species
    cdef np.ndarray timepoints
    cdef np.ndarray measurements
    cdef unsigned N #Number of samples
    cdef unsigned M #Number of measurements

    cdef np.ndarray get_measurements(self)

    cdef inline list get_measured_species(self):
        return self.measured_species

    cdef inline np.ndarray get_timepoints(self):
        return self.timepoints

    cdef inline unsigned get_N(self):
        return self.N

cdef class BulkData(Data):
    cdef unsigned multiple_timepoints
    cdef unsigned nT #Number of timepoints

    
cdef class FlowData(Data):
    pass

cdef class StochasticTrajectories(Data):
    cdef unsigned multiple_timepoints
    cdef unsigned nT #Number of timepoints

    cdef np.ndarray get_measurements(self)

##################################################                ####################################################
######################################              LIKELIHOOD                        ################################
#################################################                     ################################################

cdef class Likelihood:
    cdef double get_log_likelihood(self)

cdef class ModelLikelihood(Likelihood):
    cdef Model m
    cdef CSimInterface csim
    cdef RegularSimulator propagator
    cdef DelaySimulator propagator_delay
    cdef np.ndarray meas_indices
    cdef np.ndarray init_state_indices
    cdef np.ndarray init_state_vals
    cdef np.ndarray init_param_indices
    cdef np.ndarray init_param_vals
    cdef np.ndarray initial_states
    cdef list initial_parameters
    cdef unsigned Nx0 #number of initial conditions
    cdef unsigned N #number of samples
    cdef unsigned M #number of measurements
    cdef dict default_params
    cdef np.ndarray default_species

    cdef double get_log_likelihood(self)
    cdef np.ndarray get_initial_state(self, int n)
    cdef dict get_initial_params(self, int n)

cdef class DeterministicLikelihood(ModelLikelihood):
    cdef BulkData bd
    cdef unsigned norm_order
    cdef double get_log_likelihood(self)
    cdef double hmax


cdef class StochasticTrajectoriesLikelihood(ModelLikelihood):
    cdef StochasticTrajectories sd
    cdef unsigned initial_state_matching
    cdef unsigned N_simulations
    cdef unsigned norm_order
    cdef np.ndarray timepoints
    cdef double get_log_likelihood(self)
    cdef bool has_delay 
    

cdef class StochasticTrajectoryMomentLikelihood(StochasticTrajectoriesLikelihood):
    cdef unsigned Moments
    cdef double get_log_likelihood(self)
    

cdef class StochasticStatesLikelihood(ModelLikelihood):
    cdef FlowData fd
    cdef unsigned N_simulations
    cdef unsigned Moments
    cdef double get_log_likelihood(self)

##################################################                ####################################################
######################################              PRIORS        ######################################
#################################################                     ################################################

cdef class Prior(Likelihood):
    cdef Model m
    cdef CSimInterface csim
    cdef list parameters
    cdef Distribution dist
    cdef double get_log_likelihood(self)
   
##################################################                ####################################################
######################################              INFERENCE                         ################################
#################################################                     ################################################


cdef class DeterministicInference:
    cdef Model m
    cdef np.ndarray params_to_estimate #indices of parameters to estimate
    cdef np.ndarray global_init_state
    cdef np.ndarray global_init_params

    cdef Distribution prior

    cdef unsigned num_walkers
    cdef unsigned num_iterations
    cdef unsigned dimension

    cdef list likelihoods

    cdef double sigma



