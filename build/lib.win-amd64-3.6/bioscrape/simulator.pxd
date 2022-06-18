cimport numpy as np
from types cimport  Model, Volume, Schnitz, Lineage
from vector cimport vector

##################################################                ####################################################
######################################              DELAY QUEUE TYPES                   ##############################
#################################################                     ################################################

cdef class DelayQueue:
    """
    An interface/virtual class for keeping track of queued delayed reactions that will resolve at some future time.

    This must be subclassed by implementations of delay queues.
    """
    cdef void add_reaction(self, double time, unsigned rxn_id, double amount)
    cdef double get_next_queue_time(self)
    cdef void get_next_reactions(self, double *rxn_array)
    cdef void advance_time(self)
    cdef void set_current_time(self, double t)
    cdef DelayQueue copy(self)
    cdef DelayQueue clear_copy(self)

    cdef np.ndarray binomial_partition(self, double p)

cdef class ArrayDelayQueue(DelayQueue):
    """
    A class implementing the DelayQueue interface with an array of future times book-keeping how many of each rection
    occurs at each future time with resolution dt.

    Attributes:
        dt (double): The time resolution
        next_queue_time (double): The next gridded time point. These are spaced by dt.
        num_cols (unsigned): The number of cols in the delay queue array. This is the number of dt grid points ahead
        num_reactions (unsigned): The number of reactions in the system. Also the number of rows in the queue array.
        start_index (unsigned): The index in the queue array corresponding to the next queue time. This cycles around.
        queue (np.ndarray): 2-D array where the columns are future time points spaced by dt and the rows are future rxn
                            occurrences for each reaction. The start index is the very next time and each successive
                            index corresponds to dt more than the previous one including looping around the end of the
                            array. Thus, the max future time that can be remembered is dt*num_cols

    """
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
    """
    An interface for keeping track of the stoichiometric matrix and delay
    """
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

    cdef void apply_repeated_rules(self, double *state, double time)
    cdef void apply_repeated_volume_rules(self, double *state, double volume, double time)
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

    cdef void apply_repeated_rules(self, double *state,double time)
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
    cdef int max_species_count
    cdef int max_volume

    cdef void initialize_reaction_inputs(self)
    cdef void compute_stochastic_propensities(self, double *state, double *propensity_destination, double time)
    cdef void compute_stochastic_volume_propensities(self, double *state, double *propensity_destination, double volume, double time)
    cdef void check_count_function(self, double *state, double volume)
    
# Simulation output values here
cdef class SSAResult:
    """
    A class for keeping track of the result from a regular simulation (timepoints and species over time).
    """
    cdef np.ndarray timepoints
    cdef np.ndarray simulation_result

    cdef inline np.ndarray get_timepoints(self):
        return self.timepoints

    cdef inline  np.ndarray get_result(self):
        return self.simulation_result

cdef class DelaySSAResult(SSAResult):
    """
    A class for keeping track of the result from a delay simulation (timepoints and species and the final set
    of  queued reactions).
    """
    cdef DelayQueue final_delay_queue

    cdef inline DelayQueue get_delay_queue(self):
        return self.final_delay_queue

cdef class VolumeSSAResult(SSAResult):
    """
    A class for keeping track of the result from a volume simulation (timepoints and species/volume over time).
    """
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
    """
    A class for keeping track of the result from a volume simulation with delay (timepoints and species/volume over
    time and the final set of queued reactions).
    """
    cdef DelayQueue final_delay_queue

    cdef inline DelayQueue get_delay_queue(self):
        return self.final_delay_queue

    cdef VolumeCellState get_final_cell_state(self)


# Cell state classes

cdef class CellState:
    """
    A class for keeping track of a simple cell state (species and time).
    """
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
    """
    A class for keeping track of cell state in a delay system (species, time, and queued reactions).
    """
    cdef DelayQueue delay_queue

    cdef inline DelayQueue get_delay_queue(self):
        return self.delay_queue

    cdef inline void set_delay_queue(self, DelayQueue q):
        self.delay_queue = q

cdef class VolumeCellState(CellState):
    """
    A class for keeping track of cell state in a system with volume (species, time, and volume).
    """
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
    """
    A class for keeping track of cell state in a system with delay and volume (species, time, volume, and queued rxns)
    """
    cdef DelayQueue delay_queue

    cdef inline DelayQueue get_delay_queue(self):
        return self.delay_queue

    cdef inline void set_delay_queue(self, DelayQueue q):
        self.delay_queue = q

##################################################                ####################################################
######################################              SPLITTERS                         ################################
#################################################                     ################################################

cdef class VolumeSplitter:
    """
    Interface class for defining a method to partition cells for the case with no delay and only volume involved.
    """
    cdef np.ndarray partition(self, VolumeCellState parent)

cdef class DelayVolumeSplitter:
    """
    Interface class for defining a method to partition cells for the case with delay AND volume involved.
    """
    cdef np.ndarray partition(self, DelayVolumeCellState parent)



cdef class PerfectBinomialVolumeSplitter(VolumeSplitter):
    """
    A volume splitting class that splits the cell into two equal halves and splits the molecuels binomially w/ p = 0.5
    """
    cdef np.ndarray partition(self, VolumeCellState parent)


cdef class GeneralVolumeSplitter(VolumeSplitter):
    """
    A volume splitting class that splits the cell into two cells and can split species
    binomially, perfectly, or by duplication.
    """
    cdef vector[int] binomial_indices
    cdef vector[int] perfect_indices
    cdef vector[int] duplicate_indices
    cdef double partition_noise

    cdef np.ndarray partition(self, VolumeCellState parent)


cdef class PerfectBinomialDelayVolumeSplitter(DelayVolumeSplitter):
    """
    A class which splits a cell with delays into two equal halves and partitions molecules and queued reactions
    binomially with p=0.5.

    WARNING: If the queued reactions have negative species updates, this can lead to possible negative species count.
    """
    cdef np.ndarray partition(self, DelayVolumeCellState parent)

cdef class CustomSplitter(VolumeSplitter):
    """
    A volume splitting class that splits the cell into two equal halves and splits the molecuels binomially w/ p = 0.5
    """
    cdef object split_function
    cdef np.ndarray partition(self, VolumeCellState parent)


##################################################                ####################################################
######################################              SIMULATORS                        ################################
#################################################                     ################################################

# Regular simulations with no volume or delay involved.
cdef class RegularSimulator:
    """
    Interface class for defining simulators for the regular case with no delay/volume involved
    """
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)

cdef class DeterministicSimulator(RegularSimulator):
    """
    A class for implementing a deterministic simulator.
    """
    cdef double atol
    cdef double rtol
    cdef unsigned mxstep

    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)

    cdef inline set_tolerance(self, double atol, double rtol):
        self.atol = atol
        self.rtol = rtol


cdef class DeterministicDilutionSimulator(RegularSimulator):
    """
    A class for implementing a deterministic simulator.
    """
    cdef double atol
    cdef double rtol
    cdef unsigned mxstep
    cdef double dilution_rate

    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)

    cdef inline set_tolerance(self, double atol, double rtol):
        self.atol = atol
        self.rtol = rtol


cdef class SSASimulator(RegularSimulator):
    """
    A class for implementing a stochastic SSA simulator.
    """
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)

cdef class SafeModeSSASimulator(RegularSimulator):
    """
    A class for implementing a stochastic SSA simulator.
    """
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)

cdef class TimeDependentSSASimulator(RegularSimulator):
    """
    A class for implementing a stochastic SSA simulator.
    """
    cdef SSAResult simulate(self, CSimInterface sim, np.ndarray timepoints)



cdef class DelaySimulator:
    """
    Interface class for defining a simulator with delay.
    """
    cdef DelaySSAResult delay_simulate(self, CSimInterface sim, DelayQueue dq, np.ndarray timepoints)

cdef class DelaySSASimulator(DelaySimulator):
    """
    A class for doing delay simulations using the stochastic simulation algorithm.
    """
    cdef DelaySSAResult delay_simulate(self, CSimInterface sim, DelayQueue dq, np.ndarray timepoints)



cdef class VolumeSimulator:
    """
    Interface class for doing volume simulations.
    """
    cdef VolumeSSAResult volume_simulate(self, CSimInterface sim, Volume v, np.ndarray timepoints)


cdef class VolumeSSASimulator(VolumeSimulator):
    """
    Volume SSA implementation.
    """
    cdef VolumeSSAResult volume_simulate(self, CSimInterface sim, Volume v, np.ndarray timepoints)

cdef class DelayVolumeSimulator:
    """
    Interface class for doing simulations with delay and volume.
    """
    cdef DelayVolumeSSAResult delay_volume_simulate(self, CSimInterface sim, DelayQueue q,
                                                    Volume v, np.ndarray timepoints)

cdef class DelayVolumeSSASimulator(DelayVolumeSimulator):
    """
    SSA implementation for doing simulations with delay and volume.
    """
    cdef DelayVolumeSSAResult delay_volume_simulate(self, CSimInterface sim, DelayQueue q,
                                                    Volume v, np.ndarray timepoints)



# Simulation functions for doing cell division related stuff.
cdef Lineage simulate_cell_lineage(CSimInterface sim, Volume v, np.ndarray timepoints,
                                    VolumeSimulator vsim, VolumeSplitter vsplit)
cdef Lineage simulate_delay_cell_lineage(CSimInterface sim, DelayQueue q, Volume v, np.ndarray timepoints,
                                   DelayVolumeSimulator dvsim, DelayVolumeSplitter dvsplit)


cdef list propagate_cell(ModelCSimInterface sim, VolumeCellState cell, double end_time,
                               VolumeSimulator vsim, VolumeSplitter vsplit)
cdef list propagate_cells(ModelCSimInterface sim, list cells, double end_time,
                          VolumeSimulator vsim, VolumeSplitter vsplit)
