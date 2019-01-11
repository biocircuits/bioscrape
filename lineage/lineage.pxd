# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

cimport numpy as np
import numpy as np
from bioscrape.vector cimport vector
#from vector cimport vector 


#from bioscrape.simulator cimport VolumeCellState, VolumeSSAResult, DelayVolumeSSAResult, VolumeSplitter, CSimInterface, ModelCSimInterface, DelayVolumeSimulator, VolumeSimulator
#from bioscrape.simulator import VolumeCellState, VolumeSSAResult, DelayVolumeSSAResult, VolumeSplitter, CSimInterface, ModelCSimInterface

#from simulator import VolumeSSAResult, DelayVolumeSSAResult, VolumeSplitter, CSimInterface, ModelCSimInterface, DelayVolumeSimulator, VolumeSimulator

import warnings

#Volume Events are stochastic events which alter a single cell's volume.
cdef class Event():
	#cdef Propensity propensity
	#cdef Propensity get_propensity(self)
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices)
	cdef double evaluate_event(self, double *state, double *params, double volume, double time)

cdef class VolumeEvent(Event):
	#cdef Propensity propensity
	#cdef Propensity get_propensity(self)
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices)
	cdef double evaluate_event(self, double *state, double *params, double volume, double time)
	cdef double get_volume(self, double *state, double *params, double volume, double time)

cdef class LinearVolumeEvent(VolumeEvent):
	#cdef Propensity propensity
	cdef unsigned growth_rate_ind
	#cdef Propensity get_propensity(self)
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices)
	cdef double get_volume(self, double *state, double *params, double volume, double time)
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices):


cdef class MultiplicativeVolumeEvent(VolumeEvent):
	#cdef Propensity propensity
	cdef unsigned growth_rate_ind
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices)
	cdef double get_volume(self, double *state, double *params, double volume, double time)

cdef class GeneralVolumeEvent(VolumeEvent):
	#cdef Propensity propensity
	cdef Term volume_equation
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices)
	cdef double get_volume(self, double *state, double *params, double volume, double time)

#Division Events are stochastic events which cause a cell to divide using a particular VolumeSplitter
cdef class DivisionEvent(Event):

	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices)
	cdef double evaluate_event(self, double *state, double *params, double volume, double time)

#Death Events are stochastic events which cause a cell to die
cdef class DeathEvent(Event):
	#cdef Propensity propensity
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices)
	cdef double evaluate_event(self, double *state, double *params, double volume, double time)

#A superclass for all lineage rules - mostly used to keep code cleaner
cdef class LineageRule()

#Volume Rules occur every dt (determined by the simulation timepoints) and update the volume
cdef class VolumeRule(LineageRule):
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt)

cdef class LinearVolumeRule(VolumeRule):
	cdef unsigned growth_rate_ind
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt)

cdef class MultiplicativeVolumeRule(VolumeRule):
	cdef unsigned growth_rate_ind
	cdef double get_volume(self, double *state, double *params, double volume, double time, doublt dt)

cdef class AssignmentVolumeRule(VolumeRule):
	cdef Term volume_equation
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt)

cdef class ODEVolumeRule(VolumeRule):
	cdef Term volume_equation
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt)

#Division rules are checked at the beginning of each simulation loop to see if a cell should divide
cdef class DivisionRule(LineageRule):

	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume)

#A division rule where division occurs after some amount of time (with an optional noise term)
cdef class TimeDivisionRule(DivisionRule):
	cdef unsigned threshold_ind
	cdef unsigned threshold_noise_ind
	cdef unsigned has_noise
	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume)


#A division rule where division occurs at some volume threshold (with an optional noise term)
cdef class VolumeDivisionRule(DivisionRule):
	cdef unsigned threshold_ind
	cdef unsigned threshold_noise_ind
	cdef unsigned has_noise
	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume)


#A division rule where division occurs after the cell has grown by some amount delta (with an optional noise term)
cdef class DeltaDivisionRule(DivisionRule):
	cdef unsigned threshold_ind
	cdef unsigned threshold_noise_ind
	cdef unsigned has_noise
	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume)


#A general division rule
cdef class GeneralDivisionRule(DivisionRule):
	cdef Term equation
	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume)


#Death rules are checked at the beginning of each simulation loop to see if a cell should die
cdef class DeathRule(LineageRule):
	cdef unsigned check_dead(self, double *state, double *params, double time, initial_time, initial_volume)

#A death rule where death occurs when some species is greater than, less than, or equal to a given threshold
cdef class SpeciesDeathRule(DeathRule):
	cdef double threshold
	cdef unsigned species_ind
	cdef unsigned threshold_ind
	cdef int comparison
	cdef unsigned threshold_noise_ind
	cdef unsigned has_noise
	cdef unsigned check_dead(self, double *state, double *params, double time, double volume, initial_time, initial_volume)

cdef class ParamDeathRule(DeathRule):
	cdef double threshold
	cdef unsigned param_ind
	cdef unsigned threshold_ind
	cdef int comparison
	cdef unsigned threshold_noise_ind
	cdef unsigned has_noise
	cdef unsigned check_dead(self, double *state, double *params, double time, double volume, initial_time, initial_volume)

#A general death rule
cdef class GeneralDeathRule(DeathRule):
	cdef Term equation
	cdef unsigned check_dead(self, double *state, double *params, double time, double volume, initial_time, initial_volume)

#A subclass of Model with lineage functionality. Major differences is propensities now include event propensities.
#contains appropriate getters and setters and initializers in python.
cdef class LineageModel(Model):
	cdef unsigned num_division_events
	cdef unsigned num_division_rules
	cdef unsigned num_death_events
	cdef unsigned num_death_rules
	cdef unsigned num_volume_events
	cdef unsigned num_volume_rules

	cdef list volume_events_list
	cdef list division_events_list
	cdef list death_events_list

	cdef vector[void*] c_lineage_propensities
	cdef list lineage_propensities
	cdef vector[void*] c_death_events
	cdef list death_events
	cdef vector[void*] c_division_events
	cdef list division_events
	cdef vector[void*] c_volume_events
	cdef list volume_events
	cdef vector[void*] c_other_events
	cdef list other_events
	cdef vector[void*] c_death_rules
	cdef list death_rules
	cdef vector[void*] c_division_rules
	cdef list division_rules
	cdef vector[void*] c_volume_rules
	cdef list volume_rules
	cdef list division_event_volume_splitters
	cdef list division_rule_volume_splitters

	cdef void _initialize_lineage(self)

	cdef (vector[void*])* get_c_division_rules(self)
	cdef (vector[void*])* get_c_volume_rules(self)
	cdef (vector[void*])* get_c_death_rules(self)
	cdef (vector[void*])* get_c_division_events(self)
	cdef (vector[void*])* get_c_volume_events(self)
	cdef (vector[void*])* get_c_death_events(self)
	cdef unsigned get_num_division_rules(self)
	cdef unsigned get_num_volume_rules(self)
	cdef unsigned get_num_death_rules(self)
	cdef unsigned get_num_division_events(self)
	cdef unsigned get_num_volume_events(self)
	cdef unsigned get_num_death_events(self)

cdef class LineageCSimInterface(ModelCSimInterface):
	cdef unsigned num_division_rules
	cdef unsigned num_volume_rules
	cdef unsigned num_death_rules
	cdef unsigned num_division_events
	cdef unsigned num_volume_events
	cdef unsigned num_death_events
	cdef unsigned num_propensities
	cdef list lineage_propensities
	cdef vector[void*] *c_lineage_propensities
	cdef vector[void*] *c_volume_rules
	cdef vector[void*] *c_division_rules
	cdef vector[void*] *c_death_rules
	cdef vector[void*] *c_volume_events
	cdef vector[void*] *c_division_events
	cdef vector[void*] *c_death_events
	cdef unsigned ind

	#similar to compute_stochastic_volume_propensities but events are new included as well
	cdef void compute_lineage_propensities(self, double *state, double *propensity_destination, double volume, double time)

	#Applies volume rules and returns a new volume (double). Volume rules are applied in the order they were added to the model
	cdef double apply_volume_rules(self, double *state, double volume, double time, double dt)

	#Applies death rules in the order they were added to the model. Returns the index of the first death rule that returns True. -1 otherwise.
	cdef int apply_death_rules(self, double *state, double volume, double time, LineageVolumeCellState v)

	#Applies divison rules in the order they were added to the model. Returns the index of the first division rule that returns True. -1 otherwise
	cdef int apply_division_rules(self, double *state, double volume, double time, LineageVolumeCellState v)

	cdef double apply_volume_event(unsigned event_index, double *state, double current_time, double current_volume)

	cdef np.ndarray divide_cell(unsigned vsplit_index, double *state, double current_time, double current_volume)


#A new wrapper for the VolumeCellState with new internal variables
cdef class LineageVolumeCellState(DelayVolumeCellState):
	cdef double initial_volume
	cdef double initial_time
	cdef void set_initial_vars(self, double volume, double time)
	cdef double get_initial_volume(self)
	cdef double get_initial_time(self)
	 
cdef class SingleCellSSAResult(VolumeSSAResult):
	cdef int dead
	cdef divided
	cdef void set_dead(self, int dead)
	cdef int get_dead(Self)
	cdef void set_divided(self, int divided)
	cdef int get_divided(self, divided)

#Simulator Classes for Cell Lineages
cdef class LineageSSASimulator():
	
	#SSA for a single cell. Simulates until it devides or dies using division / death rules and/or reactions.
	cdef SingleCellSSAResult SimulateSingleCell(self, LineageCSimInterface sim, LineageVolumeCellState v, np.ndarray timepoints)

	cdef Lineage SimulateCellLineage(self, LineageCSimInterface sim, Volume v, np.ndarray timepoints)

	cdef np.ndarray partition(self, unsigned vsplit_ind, LineageVolumeCellState parent)

	cdef np.ndarray delay_partition(self, unsigned vsplit_ind, LineageVolumeCellState parent)

