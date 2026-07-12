# NOTE: class/method docstrings live in simulator.pyx, not here. If
# both files define a docstring for the same class/method, the .pyx
# one silently wins with no warning at compile time.

cimport numpy as np
from types cimport  Model, Volume, Schnitz, Lineage
from vector cimport vector

##################################################                ####################################################
######################################              DELAY QUEUE TYPES                   ##############################
#################################################                     ################################################

cdef class DelayQueue:
    cdef void add_reaction(self, double time, unsigned rxn_id, double amount)
    cdef double get_next_queue_time(self)
    cdef void get_next_reactions(self, double *rxn_array)
    cdef void advance_time(self)
    cdef void set_current_time(self, double t)
    cdef DelayQueue copy(self)
    cdef DelayQueue clear_copy(self)

    cdef np.ndarray binomial_partition(self, double p)

cdef class ArrayDelayQueue(DelayQueue):
    cdef double dt
    cdef double next_queue_time
    cdef unsigned num_cols
    cdef unsigned num_reactions
    cdef unsigned start_index
    cdef np.ndarray queue

    cdef void add_reaction(self, double time, unsigned rxn_id, double amount)
    cdef double get_next_queue_time(self)
    cdef void get_next_reactions(self, double *rxn_array)
    cdef void advance_time(self)
    cdef DelayQueue copy(self)
    cdef DelayQueue clear_copy(self)

    cdef np.ndarray binomial_partition(self, double p)

##################################################                ####################################################
######################################              SIMULATION INTERFACES               ##############################
#################################################                     ################################################


cdef class CSimInterface:
    cdef np.ndarray update_array
    cdef np.ndarray delay_update_array

    cdef unsigned delay_flag
    cdef np.ndarray initial_state
    cdef double initial_time
    cdef double dt

    cdef unsigned num_reactions
    cdef unsigned num_species

    # Stuff for deterministic simulation to be fast
    cdef np.ndarray propensity_buffer
    cdef vector[vector[int]] S_indices
    cdef vector[vector[int]] S_values

    cdef void prep_deterministic_simulation(self)
    cdef void calculate_deterministic_derivative(self, double *x, double *dxdt, double t)
    # end of deterministic simulation stuff

    #method meant to be overwritten to check if interfaces/models are correct. called by simulators.
    cdef void check_interface(self)

    cdef np.ndarray get_update_array(self)
    cdef np.ndarray get_delay_update_array(self)
    cdef double compute_delay(self, double *state, unsigned rxn_index)
    cdef void compute_propensities(self, double *state, double *propensity_destination, double time)
    cdef void compute_volume_propensities(self, double *state, double *propensity_destination, double volume, double time)
    cdef void compute_stochastic_propensities(self, double *state, double *propensity_destination, double time)
    cdef void compute_stochastic_volume_propensities(self, double *state, double *propensity_destination, double volume, double time)
    cdef unsigned requires_delay(self)

    cdef void apply_repeated_rules(self, double *state, double time, unsigned rule_step)
    cdef void apply_repeated_volume_rules(self, double *state, double volume, double time, unsigned rule_step)
    cdef unsigned get_number_of_rules(self)

    cdef np.ndarray get_initial_state(self)
    cdef void set_initial_state(self, np.ndarray a)

    cdef double get_initial_time(self)
    cdef void set_initial_time(self, double t)

    cdef unsigned get_num_reactions(self)
    cdef unsigned get_num_species(self)

    cdef void set_dt(self, double dt)
    cdef double get_dt(self)

    cdef double* get_param_values(self)
    cdef unsigned get_num_parameters(self)


cdef class ModelCSimInterface(CSimInterface):
    cdef Model model
    cdef vector[void*] *c_delays
    cdef vector[void*] *c_propensities
    cdef vector[void*] *c_repeat_rules
    cdef double *c_param_values
    cdef np.ndarray np_param_values

    cdef double compute_delay(self, double *state, unsigned rxn_index)
    cdef void compute_propensities(self, double *state, double *propensity_destination, double time)
    cdef void compute_volume_propensities(self, double *state, double *propensity_destination, double volume, double time)
    cdef void compute_stochastic_propensities(self, double *state, double *propensity_destination, double time)
    cdef void compute_stochastic_volume_propensities(self, double *state, double *propensity_destination, double volume, double time)
    cdef np.ndarray get_initial_state(self)
    cdef void set_param_values(self, np.ndarray params)

    cdef void apply_repeated_rules(self, double *state,double time, unsigned rule_step)
    cdef unsigned get_number_of_rules(self)
    cdef unsigned get_number_of_species(self)
    cdef unsigned get_number_of_reactions(self)

    cdef double* get_param_values(self)
    cdef unsigned get_num_parameters(self)


cdef class SafeModelCSimInterface(ModelCSimInterface):
    #Two not instantiate these everytime a propensity is computed
    cdef unsigned rxn_ind
    cdef unsigned s_ind
    cdef unsigned prop_is_0
    cdef int[:, :, :] reaction_input_indices
    cdef double max_species_count
    cdef double max_volume

    cdef void initialize_reaction_inputs(self)
    cdef void compute_stochastic_propensities(self, double *state, double *propensity_destination, double time)
    cdef void compute_stochastic_volume_propensities(self, double *state, double *propensity_destination, double volume, double time)
    cdef void check_count_function(self, double *state, double volume)
    
# Simulation output values here
cdef class SSAResult:
    cdef np.ndarray timepoints
    cdef np.ndarray simulation_result

    cdef np.ndarray empirical_distribution(self, double burn_in, list species_inds, double tend, list max_counts_list)
    cdef np.ndarray first_moment(self, double start_time, double final_time, list species_inds)
    cdef np.ndarray standard_deviation(self, double start_time, double final_time, list species_inds)
    cdef np.ndarray second_moment(self, double start_time, double final_time, list species_inds1, list species_inds2)
    cdef np.ndarray correlations(self, double start_time, double final_time, list species_inds1, list species_inds2)

    cdef inline np.ndarray get_timepoints(self):
        return self.timepoints

    cdef inline  np.ndarray get_result(self):
        return self.simulation_result

cdef class DelaySSAResult(SSAResult):
    cdef DelayQueue final_delay_queue

    cdef inline DelayQueue get_delay_queue(self):
        return self.final_delay_queue

cdef class VolumeSSAResult(SSAResult):
    cdef np.ndarray volume
    cdef unsigned cell_divided_flag
    cdef unsigned cell_dead_flag
    cdef Volume volume_object

    cdef inline unsigned cell_divided(self):
        return self.cell_divided_flag
    cdef inline unsigned cell_dead(self):
        return self.cell_dead_flag
    cdef inline void set_cell_dead(self, unsigned dead):
        self.cell_dead_flag = dead
        
    cdef inline np.ndarray get_volume(self):
        return self.volume

    cdef inline Volume get_volume_object(self):
        return self.volume_object

    cdef inline void set_volume_object(self, Volume v):
        self.volume_object = v

    cdef VolumeCellState get_final_cell_state(self)
    cdef Schnitz get_schnitz(self)


cdef class DelayVolumeSSAResult(VolumeSSAResult):
    cdef DelayQueue final_delay_queue

    cdef inline DelayQueue get_delay_queue(self):
        return self.final_delay_queue

    cdef VolumeCellState get_final_cell_state(self)


# Cell state classes

cdef class CellState:
    cdef np.ndarray state
    cdef double time

    cdef inline void set_state(self, np.ndarray state):
        self.state = state

    cdef inline np.ndarray get_state(self):
        return self.state

    cdef inline void set_time(self, double t):
        self.time = t

    cdef inline double get_time(self):
        return self.time

cdef class DelayCellState(CellState):
    cdef DelayQueue delay_queue

    cdef inline DelayQueue get_delay_queue(self):
        return self.delay_queue

    cdef inline void set_delay_queue(self, DelayQueue q):
        self.delay_queue = q

cdef class VolumeCellState(CellState):
    cdef double volume
    cdef Volume volume_object
    
    cdef inline void set_volume(self, double volume):
        self.volume = volume

    cdef inline double get_volume(self):
        return self.volume

    cdef inline void set_volume_object(self, Volume vol):
        self.volume_object = vol

    cdef inline Volume get_volume_object(self):
        return self.volume_object


cdef class DelayVolumeCellState(VolumeCellState):
    cdef DelayQueue delay_queue

    cdef inline DelayQueue get_delay_queue(self):
        return self.delay_queue

    cdef inline void set_delay_queue(self, DelayQueue q):
        self.delay_queue = q

##################################################                ####################################################
######################################              SPLITTERS                         ################################
#################################################                     ################################################

cdef class VolumeSplitter:
    cdef np.ndarray partition(self, VolumeCellState parent)

cdef class DelayVolumeSplitter:
    cdef np.ndarray partition(self, DelayVolumeCellState parent)



cdef class PerfectBinomialVolumeSplitter(VolumeSplitter):
    cdef np.ndarray partition(self, VolumeCellState parent)


cdef class GeneralVolumeSplitter(VolumeSplitter):
    cdef vector[int] binomial_indices
    cdef vector[int] perfect_indices
    cdef vector[int] duplicate_indices
    cdef double partition_noise

    cdef np.ndarray partition(self, VolumeCellState parent)


cdef class PerfectBinomialDelayVolumeSplitter(DelayVolumeSplitter):
    cdef np.ndarray partition(self, DelayVolumeCellState parent)

cdef class CustomSplitter(VolumeSplitter):
    cdef object split_function
    cdef np.ndarray partition(self, VolumeCellState parent)


##################################################                ####################################################
######################################              SIMULATORS                        ################################
#################################################                     ################################################

# Regular simulations with no volume or delay involved.
cdef class RegularSimulator:
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)

cdef class DeterministicSimulator(RegularSimulator):
    cdef double atol
    cdef double rtol
    cdef double hmax
    cdef unsigned mxstep

    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)

    cdef inline set_tolerance(self, double atol, double rtol):
        self.atol = atol
        self.rtol = rtol


cdef class SSASimulator(RegularSimulator):
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)


cdef class DelaySimulator:
    cdef DelaySSAResult delay_simulate(self, CSimInterface sim, DelayQueue dq, np.ndarray timepoints)

cdef class DelaySSASimulator(DelaySimulator):
    cdef DelaySSAResult delay_simulate(self, CSimInterface sim, DelayQueue dq, np.ndarray timepoints)



cdef class VolumeSimulator:
    cdef VolumeSSAResult volume_simulate(self, CSimInterface sim, Volume v, np.ndarray timepoints)


cdef class VolumeSSASimulator(VolumeSimulator):
    cdef VolumeSSAResult volume_simulate(self, CSimInterface sim, Volume v, np.ndarray timepoints)

cdef class DelayVolumeSimulator:
    cdef DelayVolumeSSAResult delay_volume_simulate(self, CSimInterface sim, DelayQueue q,
                                                    Volume v, np.ndarray timepoints)

cdef class DelayVolumeSSASimulator(DelayVolumeSimulator):
    cdef DelayVolumeSSAResult delay_volume_simulate(self, CSimInterface sim, DelayQueue q,
                                                    Volume v, np.ndarray timepoints)
