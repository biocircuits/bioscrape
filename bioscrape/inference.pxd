cimport numpy as np
from types import Model
from types cimport Model
from simulator cimport CSimInterface, RegularSimulator, ModelCSimInterface, DeterministicSimulator, SSASimulator
from simulator import CSimInterface, RegularSimulator, ModelCSimInterface, DeterministicSimulator, SSASimulator


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
    cdef unsigned M #Number of measured species
    cdef list measured_species
    cdef np.ndarray timepoints
    cdef np.ndarray measurements

    cdef inline np.ndarray get_measurements(self):
        return self.measurements

    cdef inline list get_measured_species(self):
        return self.measured_species

    cdef inline np.ndarray get_timepoints(self):
        return self.timepoints

cdef class BulkData(Data):
    cdef unsigned N #Number of Samples (Assumed to be 1 currently)

cdef class FlowData(Data):
    cdef unsigned N #Number of Samples

cdef class StochasticTrajectories(Data):
    cdef unsigned N #Number of Samples


##################################################                ####################################################
######################################              LIKELIHOOD                        ################################
#################################################                     ################################################

cdef class Likelihood:
    cdef double get_log_likelihood(self)

cdef class ModelLikelihood(Likelihood):
    cdef Model m
    cdef CSimInterface csim
    cdef RegularSimulator propagator

    cdef double get_log_likelihood(self)

cdef class DeterministicLikelihood(ModelLikelihood):
    cdef np.ndarray init_state_indices
    cdef np.ndarray init_state_vals
    cdef np.ndarray init_param_indices
    cdef np.ndarray init_param_vals

    cdef BulkData bd
    cdef np.ndarray meas_indices

    cdef double get_log_likelihood(self)


cdef class StochasticTrajectoriesLikelihood(ModelLikelihood):

    cdef double get_log_likelihood(self)

cdef class StochasticStatesLikelihood(ModelLikelihood):

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



