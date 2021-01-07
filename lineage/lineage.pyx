# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
from libc.math cimport round
cimport numpy as np
from bioscrape.simulator cimport DeterministicSimulator, VolumeCellState, \
								 VolumeSSAResult, DelayVolumeSSAResult, \
								 VolumeSplitter, CSimInterface, \
								 ModelCSimInterface, DelayVolumeSSASimulator, \
								 VolumeSSASimulator, DelayQueue, \
								 DelayVolumeCellState
from bioscrape.simulator import DeterministicSimulator, VolumeCellState, \
								VolumeSSAResult, DelayVolumeSSAResult, \
								VolumeSplitter, CSimInterface, \
								ModelCSimInterface, DelayVolumeSSASimulator, \
								VolumeSSASimulator, DelayQueue, \
								DelayVolumeCellState

from bioscrape.types cimport Model, Volume, Schnitz, Lineage, Propensity, Term, \
							 Rule

from bioscrape.types import sympy_species_and_parameters, parse_expression, \
							Volume
from bioscrape.random cimport normal_rv

cimport bioscrape.random as cyrandom
import bioscrape.random as cyrandom

from bioscrape.vector cimport vector
#from vector cimport vector

import warnings

import pandas

# Events are general objects that happen with some internal propensity but are 
# not chemical reactions
cdef class Event:

	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices):
		raise NotImplementedError("VolumeReactions Must be subclassed")

	def get_species_and_parameters(self, dict event_fields):
		raise NotImplementedError("VolumeReactions Must be subclassed")

	#cdef Propensity get_propensity(self):
	#	return <Propensity>self.propensity

	#Meant to be overwritten if an event is supposed to do something.
	cdef double evaluate_event(self, double* state, double *params, double volume, double time):
		return 0

#Volume Events are stochastic events which alter a single cell's volume.
cdef class VolumeEvent(Event):
	cdef double evaluate_event(self, double* state, double *params, double volume, double time):
		return self.get_volume(state, params,volume, time)

	cdef double get_volume(self, double* state, double *params, double volume, double time):
		raise NotImplementedError("VolumeEvent must be subclassed")

	def get_species_and_parameters(self, dict event_fields):
		raise NotImplementedError("VolumeReactions Must be subclassed")

cdef class LinearVolumeEvent(VolumeEvent):
	cdef unsigned growth_rate_ind

	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices):
		for (k, v) in event_params.items():
			if k == "growth_rate":
				self.growth_rate_ind = parameter_indices[v]
			else:
				warnings.warn("Useless paramter for LinearVolumeEvent: "+str(k))

	cdef double get_volume(self, double* state, double *params, double volume, double time):
		return volume + params[self.growth_rate_ind]

	def get_species_and_parameters(self, dict event_fields):
		return ([], [event_fields["growth_rate"]])

cdef class MultiplicativeVolumeEvent(VolumeEvent):
	cdef unsigned growth_rate_ind
	cdef Term volume_equation


	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices):
		for (k, v) in event_params.items():
			if k == "growth_rate":
				self.growth_rate_ind = parameter_indices[v]
			else:
				warnings.warn("Useless paramter for MultiplicativeVolumeEvent: "+str(k))

	cdef double get_volume(self, double* state, double *params, double volume, double time):
		return volume*(1+params[self.growth_rate_ind])

	def get_species_and_parameters(self, dict event_fields):
		return ([], [event_fields["growth_rate"]])

cdef class GeneralVolumeEvent(VolumeEvent):
	cdef Term volume_equation

	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices):
		for (k, v) in event_params.items():
			if k == "equation":
				self.volume_equation = parse_expression(v, species_indices, parameter_indices)
			else:
				warnings.warn("Useless paramter for GeneralVolumeEvent: "+str(k))

	cdef double get_volume(self, double* state, double *params, double volume, double time):
		return (<Term>self.volume_equation).volume_evaluate(state,params,volume,time)


	def get_species_and_parameters(self, dict event_fields):
		equation_string = event_fields['equation'].strip()
		species_r, parameters_r = sympy_species_and_parameters(equation_string)
		return (species_r, parameters_r)

#Division Events are stochastic events which cause a cell to divide using a particular VolumeSplitter
cdef class DivisionEvent(Event):
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices):
		#self.propensity = propensity
		#self.propensity.initialize(propensity_params, species_indices, parameter_indices)
		#self.vsplit = vsplit
		pass

	#cdef Propensity get_propensity(self):
	#	return <Propensity>self.propensity

	#Return 1 indicated the cell divided. 0 indicates it didn't. Could be used in a subclass.
	cdef double evaluate_event(self, double* state, double *params, double volume, double time):
		return 1

	def get_species_and_parameters(self, dict event_fields):
		return ([], [])


#Death Events are stochastic events which casue a cell to die
cdef class DeathEvent(Event):
	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices):
		pass

	def get_species_and_parameters(self, dict event_fields):
		return ([], [])

	#Could be subclassed to have more complex logic for if a cell dies when this event occurs
	#0 indicates not dead. 1 indicates dead.
	cdef double evaluate_event(self, double* state, double *params, double volume, double time):
		return 1

#Dummy class to help with inheritance, compilation, and code simplification. Does nothing
cdef class LineageRule:
	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		raise NotImplementedError("LineageRule must be subclassed")

	def get_species_and_parameters(self, dict fields):
		raise NotImplementedError("LineageRule must be subclassed")

#Volume Rules occur every dt (determined by the simulation timepoints) and update the volume
cdef class VolumeRule(LineageRule):
	cdef double get_volume(self, double* state, double *params, double volume, double time, double dt):
		raise NotImplementedError("get_volume must be implemented in VolumeRule Subclasses")

cdef class LinearVolumeRule(VolumeRule):
	cdef unsigned has_noise
	cdef unsigned growth_rate_ind
	cdef unsigned noise_ind
	cdef double get_volume(self, double* state, double *params, double volume, double time, double dt):
		if self.has_noise > 0:
			return volume + (params[self.growth_rate_ind]+normal_rv(0, params[self.noise_ind]))*dt
		else:
			return volume + params[self.growth_rate_ind]*dt

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		self.has_noise = 0
		for (k, v) in param_dictionary.items():
			if k == "growth_rate":
				self.growth_rate_ind = parameter_indices[v]
			elif k == "noise":
				self.has_noise = 1
				self.noise_ind = parameter_indices[v]
			else:
				warnings.warn("Useless paramter for LinearVolumeRule: "+str(k))

	def get_species_and_parameters(self, dict fields):
		if self.has_noise > 0:
			return ([], [fields["growth_rate"], fields["noise"]])
		else:
			return ([], [fields["growth_rate"]])

cdef class MultiplicativeVolumeRule(VolumeRule):
	cdef unsigned has_noise
	cdef unsigned growth_rate_ind
	cdef unsigned noise_ind
	cdef double get_volume(self, double* state, double *params, double volume, double time, double dt):
		if self.has_noise > 0:
			return volume+volume*(params[self.growth_rate_ind]+normal_rv(0, params[self.noise_ind]))*dt
		else:
			return volume+volume*params[self.growth_rate_ind]*dt

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		self.has_noise = 0
		for (k, v) in param_dictionary.items():
			if k == "growth_rate":
				self.growth_rate_ind = parameter_indices[v]
			elif k == "noise":
				self.has_noise = 1
				self.noise_ind = parameter_indices[v]
			else:
				warnings.warn("Useless paramter for MultiplicativeVolumeRule: "+str(k))

	def get_species_and_parameters(self, dict fields):
		if self.has_noise:
			return ([], [fields["growth_rate"], fields['noise']])
		else:
			return ([], [fields["growth_rate"]])

cdef class AssignmentVolumeRule(VolumeRule):
	cdef Term volume_equation
	cdef double get_volume(self, double* state, double *params, double volume, double time, double dt):
		return (<Term>self.volume_equation).volume_evaluate(state,params,volume,time)

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "equation":
				self.volume_equation = parse_expression(v, species_indices, parameter_indices)
			else:
				warnings.warn("Useless paramter for AssignmentVolumeRule: "+str(k))

	def get_species_and_parameters(self, dict fields):
		equation_string = fields['equation'].strip()
		species, parameters = sympy_species_and_parameters(equation_string)
		return (species, parameters)

cdef class ODEVolumeRule(VolumeRule):
	cdef Term volume_equation
	cdef double get_volume(self, double* state, double *params, double volume, double time, double dt):
		return volume+(<Term>self.volume_equation).volume_evaluate(state,params,volume,time)*dt

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "equation":
				self.volume_equation = parse_expression(v, species_indices, parameter_indices)
			else:
				warnings.warn("Useless paramter for ODEVolumeRule: "+str(k))

	def get_species_and_parameters(self, dict fields):
		equation_string = fields['equation'].strip()
		species, parameters = sympy_species_and_parameters(equation_string)
		return (species, parameters)

#Division rules are checked at the beginning of each simulation loop to see if a cell should divide
cdef class DivisionRule(LineageRule):

	cdef int check_divide(self, double* state, double *params, double time, double volume, double initial_time, double initial_volume):
		raise NotImplementedError("check_divide must be implemented in DivisionRule subclasses")

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		raise NotImplementedError("DivisionRule must be subclassed!")

	def get_species_and_parameters(self, dict fields):
		raise NotImplementedError("DivisionRule must be subclassed!")

#A division rule where division occurs after some amount of time (with an optional noise term)
cdef class TimeDivisionRule(DivisionRule):
	cdef unsigned has_noise
	cdef unsigned threshold_ind
	cdef unsigned threshold_noise_ind
	cdef int check_divide(self, double* state, double *params, double time, double volume, double initial_time, double initial_volume):
		if self.has_noise > 0:
			if time-initial_time >= params[self.threshold_ind]+normal_rv(0, params[self.threshold_noise_ind]):
				return 1
			else:
				return 0
		else:
			if time-initial_time >= params[self.threshold_ind] - 1E-9:
				return 1
			else:
				return 0

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "threshold":
				self.threshold_ind = parameter_indices[v]
			elif k == "noise":
				self.threshold_noise_ind = parameter_indices[v]
				self.has_noise = 1
			else:
				warnings.warn("Useless paramter for TimeDivisionRule: "+str(k))
		if "noise" not in param_dictionary.keys():
			self.has_noise = 0

	def get_species_and_parameters(self, dict fields):
		if "noise" in fields:
			return ([], [fields["threshold"], fields["noise"]])
		else:
			return ([], [fields["threshold"]])

#A division rule where division occurs at some volume threshold (with an optional noise term)
cdef class VolumeDivisionRule(DivisionRule):
	cdef unsigned has_noise
	cdef unsigned threshold_ind
	cdef unsigned threshold_noise_ind
	cdef int check_divide(self, double* state, double *params, double time, double volume, double initial_time, double initial_volume):
		if self.has_noise > 0:
			if volume >= params[self.threshold_ind]+normal_rv(0, params[self.threshold_noise_ind]):
				return 1
			else:
				return 0
		else:
			if volume >= params[self.threshold_ind] - 1E-9:
				return 1
			else:
				return 0

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "threshold":
				self.threshold_ind = parameter_indices[v]
			elif k == "noise":
				self.threshold_noise_ind = parameter_indices[v]
				self.has_noise = 1
			else:
				warnings.warn("Useless paramter for VolumeDivisionRule: "+str(k))
		if "noise" not in param_dictionary.keys():
			self.has_noise = 0

	def get_species_and_parameters(self, dict fields):
		if "noise" in fields:
			return ([], [fields["threshold"], fields["noise"]])
		else:
			return ([], [fields["threshold"]])

#A division rule where division occurs after the cell has grown by some amount delta (with an optional noise term)
cdef class DeltaVDivisionRule(DivisionRule):
	cdef unsigned has_noise
	cdef unsigned threshold_ind
	cdef unsigned threshold_noise_ind
	cdef int check_divide(self, double* state, double *params, double time, double volume, double initial_time, double initial_volume):
		if self.has_noise > 0:
			if volume - initial_volume >= params[self.threshold_ind]+normal_rv(0, params[self.threshold_noise_ind]):
				return 1
			else:
				return 0
		else:
			if volume - initial_volume >= params[self.threshold_ind] - 1E-9:
				return 1
			else:
				return 0

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "threshold":
				self.threshold_ind = parameter_indices[v]
			elif k == "noise":
				self.threshold_noise_ind = parameter_indices[v]
				self.has_noise = 1
			else:
				warnings.warn("Useless paramter for DeltaVDivisionRule: "+str(k))
		if "noise" not in param_dictionary.keys():
			self.has_noise = 0

	def get_species_and_parameters(self, dict fields):
		if "noise" in fields:
			return ([], [fields["threshold"], fields["noise"]])
		else:
			return ([], [fields["threshold"]])

#A general division rule
#returns 1 if equation(state, params, volume, time) > 0
cdef class GeneralDivisionRule(DivisionRule):
	cdef Term equation
	cdef int check_divide(self, double* state, double *params, double time, double volume, double initial_time, double intial_volume):
		if (<Term> self.equation).volume_evaluate(state,params,volume,time) > 0:
			return 1
		else:
			return 0

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "equation":
				self.equation = parse_expression(v, species_indices, parameter_indices)
			else:
				warnings.warn("Useless paramter for GeneralDivisionRule: "+str(k))

	def get_species_and_parameters(self, dict fields):
		equation_string = fields['equation'].strip()
		species, parameters = sympy_species_and_parameters(equation_string)
		return (species, parameters)

#Death rules are checked at the beginning of each simulation loop to see if a cell should die
cdef class DeathRule(LineageRule):
	cdef int check_dead(self, double* state, double *params, double time, double volume, double initial_time, double initial_volume):
		raise NotImplementedError("check_dead must be implemented in DeathRule subclasses.")

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		raise NotImplementedError("DeathRule must be Subclassed!")

	def get_species_and_parameters(self, dict fields):
		raise NotImplementedError("DeathRule must be Subclassed!")

#A death rule where death occurs when some species is greater than or less than a given threshold
cdef class SpeciesDeathRule(DeathRule):
	cdef unsigned has_noise
	cdef unsigned species_ind
	cdef unsigned threshold_ind
	cdef unsigned threshold_noise_ind
	cdef int comp
	cdef double threshold
	cdef int check_dead(self, double* state, double *params, double time, double volume, double initial_time, double initial_volume):
		self.threshold = params[self.threshold_ind]
		if self.has_noise > 0:
			self.threshold = self.threshold + normal_rv(0, params[self.threshold_noise_ind])

		if state[self.species_ind] > self.threshold - 1E-9 and state[self.species_ind] < self.threshold + 1E-9 and self.comp == 0:
			return 1
		elif state[self.species_ind] > self.threshold-1E-9 and self.comp == 1:
			return 1
		elif state[self.species_ind] < self.threshold+ 1E-9 and self.comp == -1:
			return 1
		else:
			return 0

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "specie":
				self.species_ind = species_indices[v]
			elif k == "threshold":
				self.threshold_ind = parameter_indices[v]
			elif k == "noise":
				self.threshold_noise_ind = parameter_indices[v]
				self.has_noise = 1
			elif k == "comp":
				if v == "=" or v =="equal":
					self.comp = 0
				elif v == "<" or v == "less":
					self.comp = -1
				elif v == ">" or v == "greater":
					self.comp = 1
			else:
				warnings.warn("Useless paramter for SpeciesDeathRule: "+str(k))
		if "noise" not in param_dictionary.keys():
			self.has_noise = 0
		if "comp" not in param_dictionary.keys():
			warnings.warn("No comparison time added for SpeciesDeathRule in param dictionary. Defaulting to >.")
			self.comp = 1

	def get_species_and_parameters(self, dict fields):
		if "noise" in fields:
			return ([fields["specie"]], [fields["threshold"], fields["noise"]])
		else:
			return ([fields["specie"]], [fields["threshold"]])

#A death rule where death occurs when some parameter is greater than or less than a given threshold
cdef class ParamDeathRule(DeathRule):
	cdef unsigned param_ind
	cdef unsigned threshold_ind
	cdef unsigned threshold_noise_ind
	cdef unsigned has_noise
	cdef int comp
	cdef double threshold

	cdef int check_dead(self, double* state, double *params, double time, double volume, double initial_time, double initial_volume):
		self.threshold = params[self.threshold_ind]
		if self.has_noise > 0:
			self.threshold = self.threshold + normal_rv(0, params[self.threshold_noise_ind])

		if params[self.param_ind] > self.threshold - 1E-9  and params[self.param_ind] < self.threshold + 1E-9 and self.comp == 0:
			return 1
		elif params[self.param_ind] > self.threshold - 1E-9 and self.comp == 1:
			return 1
		elif params[self.param_ind] < self.threshold + 1E-9 and self.comp == -1:
			return 1
		else:
			return 0

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "param":
				self.param_ind = parameter_indices[v]
			elif k == "threshold":
				self.threshold_ind = parameter_indices[v]
			elif k == "noise":
				self.threshold_noise_ind = parameter_indices[v]
				self.has_noise = 1
			elif k == "comp":
				if v == "=" or v =="equal":
					self.comp = 0
				elif v == "<" or v == "less":
					self.comp = -1
				elif v == ">" or v == "greater":
					self.comp = 1
			else:
				warnings.warn("Useless paramter for ParamDeathRule: "+str(k))
		if "noise" not in param_dictionary.keys():
			self.has_noise = 0
		if "comp" not in param_dictionary.keys():
			warnings.warn("No comparison time added for ParamDeathRule in param dictionary. Defaulting to >.")
			self.comp = 1

	def get_species_and_parameters(self, dict fields):
		if "noise" in fields:
			return ([], [fields["param"], fields["threshold"], fields["noise"]])
		else:
			return ([], [fields["param"], fields["threshold"]])

#A general death rule. Returns 1 if equation > 0. 0 otherwise
cdef class GeneralDeathRule(DeathRule):
	cdef Term equation
	cdef int check_dead(self, double* state, double *params, double time, double volume, double initial_time, double initial_volume):
		if (<Term> self.equation).volume_evaluate(state,params,volume,time) > 0:
			return 1
		else:
			return 0

	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		for (k, v) in param_dictionary.items():
			if k == "equation":
				self.equation = parse_expression(v, species_indices, parameter_indices)
			else:
				warnings.warn("Useless paramter for GeneralDivisionRule: "+str(k))

	def get_species_and_parameters(self, dict fields):
		equation_string = fields['equation'].strip()
		species, parameters = sympy_species_and_parameters(equation_string)
		return (species, parameters)

#A super class of Bioscrape.Types.Model which contains new cell lineage features
cdef class LineageModel(Model):
	cdef unsigned num_division_events
	cdef unsigned num_division_rules
	cdef unsigned num_death_events
	cdef unsigned num_death_rules
	cdef unsigned num_volume_events
	cdef unsigned num_volume_rules

	cdef list volume_events_list
	cdef list division_events_list
	cdef list division_rules_list
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
	cdef list rule_volume_splitters
	cdef list event_volume_splitters

	cdef list global_species
	cdef double global_volume

	def __init__(self, filename = None, species = [], reactions = [], parameters = [], rules = [], events = [], sbml_filename = None, initial_condition_dict = None, input_printout = False, initialize_model = True):


		self.volume_events_list = []
		self.volume_events = []
		self.division_events_list = []
		self.division_events = []
		self.death_events_list = []
		self.death_events = []
		self.volume_rules = []
		self.death_rules = []
		self.division_rules = []
		self.division_rules_list = []
		self.num_volume_rules = 0
		self.num_death_rules = 0
		self.num_division_rules = 0
		self.num_volume_events = 0
		self.num_death_events = 0
		self.num_division_events = 0

		#Filter out new rule types before calling super
		original_rules = []
		for rule in rules:
			if len(rule) == 2:
				rule_type, rule_attributes = rule
			elif len(rule) == 3:
				rule_type, rule_attributes, rule_frequency = rule

			if not ("Volume" in rule_type or "volume" in rule_type or
				"Death" in rule_type or "death" in rule_type or
				"Division" in rule_type or "division" in rule_type):
				original_rules.append(rule)



		#Call super constructor
		super().__init__(filename = filename, species = species, reactions = reactions, parameters = parameters, rules = original_rules, initial_condition_dict = initial_condition_dict, sbml_filename = sbml_filename,  input_printout = input_printout, initialize_model = False)

		#Add new types to the model
		for rule in rules:
			if len(rule) == 2:
				rule_type, rule_attributes = rule
			elif len(rule) == 3:
				rule_type, rule_attributes, rule_frequency = rule

			if "Volume" in rule_type or "volume" in rule_type:
				self.create_volume_rule(rule_type, rule_attributes)
			elif "Death" in rule_type or "death" in rule_type:
				self.create_death_rule(rule_type, rule_attributes)
			elif "Division" in rule_type or "division" in rule_type:
				self.create_division_rule(rule_type, rule_attributes)

		for event in events:
			if len(event) == 4:
				event_type, event_params, event_propensity, propensity_params = event
			else:
				raise ValueError("Events must be tuples: (event_type (str), event_params (dict), event_propensity (str), propensity_params (dict)).")
			if "Volume" in event_type or "volume" in event_type:
				self.create_volume_event(event_type, event_params, event_propensity, propensity_params)
			elif "Death" in event_type or "death" in event_type:
				self.create_death_event(event_type, event_params, event_propensity, propensity_params)
			elif "Division" in event_type or "division" in event_type:
				self.create_division_event(event_type, event_params, event_propensity, propensity_params)
			else:
				raise ValueError("Unknown Event Type:", event_type)

		if initialize_model:
			self._initialize()


	def _create_vectors(self):
		#Create c-vectors of different objects
		super()._create_vectors()

		for rule_object in self.repeat_rules:
			self.c_repeat_rules.push_back(<void*> rule_object)

		self.num_volume_rules = len(self.volume_rules)
		for i in range(self.num_volume_rules):
			rule = self.volume_rules[i]
			self.c_volume_rules.push_back(<void*> rule)

		self.num_death_rules = len(self.death_rules)
		for i in range(self.num_death_rules):
			rule = self.death_rules[i]
			self.c_death_rules.push_back(<void *> rule)


		self.rule_volume_splitters = []
		self.event_volume_splitters = []
		self.num_division_rules = len(self.division_rules_list)
		for i in range(self.num_division_rules):
			rule, volume_splitter = self.division_rules_list[i]
			self.division_rules.append(rule)
			self.c_division_rules.push_back(<void*>rule)
			self.rule_volume_splitters.append(volume_splitter)

		#Propensity Order:
		# Reactionns, Divison Events, Volume Events, Death Events
		self.lineage_propensities = []

		self.num_division_events = len(self.division_events_list)
		for i in range(self.num_division_events):
			event, prop_object, volume_splitter = self.division_events_list[i]
			self.division_events.append(event)
			self.c_division_events.push_back(<void*> event)
			self.lineage_propensities.append(prop_object)
			self.c_lineage_propensities.push_back(<void*>prop_object)
			self.event_volume_splitters.append(volume_splitter)

		self.num_volume_events = len(self.volume_events_list)
		for i in range(self.num_volume_events):
			event, prop_object = self.volume_events_list[i]
			self.lineage_propensities.append(prop_object)
			self.c_lineage_propensities.push_back(<void*>prop_object)
			self.volume_events.append(event)
			self.c_volume_events.push_back(<void*>event)

		self.num_death_events = len(self.death_events_list)
		for i in range(self.num_death_events):
			event, prop_object = self.death_events_list[i]
			self.lineage_propensities.append(prop_object)
			self.c_lineage_propensities.push_back(<void*>prop_object)
			self.death_events.append(event)
			self.c_death_events.push_back(<void*>event)


	def py_initialize(self):
		self._initialize()
		self.initialized = True

	def py_get_event_counts(self):
		return self.num_division_events, self.num_volume_events, self.num_death_events
	def py_get_rule_counts(self):
		return self.num_division_rules, self.num_volume_rules, self.num_death_rules

	def add_event(self, Event event_object, dict event_param_dict, 
				  Propensity prop_object, dict propensity_param_dict, 
				  str event_type = None, VolumeSplitter volume_splitter = None):
		self.initialized = False

		species_names_e, param_names_e = \
					event_object.get_species_and_parameters(event_param_dict)
		species_names_p, param_names_p = \
				prop_object.get_species_and_parameters(propensity_param_dict)

		for species_name in species_names_e+species_names_p:
			self._add_species(species_name)
		for param_name in param_names_e+param_names_p:
			self._add_param(param_name)

		event_object.initialize(event_param_dict, self.species2index, self.params2index)
		prop_object.initialize(propensity_param_dict, self.species2index, self.params2index)

		if event_type in ["division", "Division", "division event", "DivisionEvent", "Division Event"]:
			if volume_splitter == None:
				raise ValueError("DivisionRules require a VolumeSplitter Object to be passed into add_event")
			self.division_events_list.append((event_object, prop_object, volume_splitter))

		elif event_type in ["volume", "Volume", "volume event", "VolumeEvent", "Volume Event"]:
			self.volume_events_list.append((event_object, prop_object))
		elif event_type in ["death", "Death", "death event", "DeathEvent", "Death Event"]:
			self.death_events_list.append((event_object, prop_object))
		else:
			raise ValueError("Unknown Event Type: Misc Event Not Yet Implemented")
			self.other_events_list.append((event_object, prop_object))

	def create_death_event(self, str event_type, dict event_params, str event_propensity_type, dict propensity_params, print_out = False):
		event_params = dict(event_params)
		propensity_params = dict(propensity_params)
		prop_object = self.create_propensity(event_propensity_type, propensity_params, print_out = print_out)
		if event_type.lower() in ["", "death", "deathevent", "death event", "default"]:
			event_object = DeathEvent()
		else:
			raise ValueError("Unknwown DeathEvent type"+str(event_type))
		self.add_event(event_object, event_params, prop_object, propensity_params, event_type = "death")

	def create_division_event(self, str event_type, dict event_params, str event_propensity_type, dict propensity_params, VolumeSplitter volume_splitter, print_out = False):
		if print_out:
			print("Adding New DivisionEvent with event_type=", event_type, "params=", event_params, "propensity_type=",event_propensity_type, "propensity_params=", propensity_params, "and VolumeSplitter=", volume_splitter)
		event_params = dict(event_params)
		propensity_params = dict(propensity_params)
		prop_object = self.create_propensity(event_propensity_type, propensity_params, print_out = print_out)
		if event_type.lower() in ["", "division", "divisionevent", "division event", "default"]:
			event_object = DivisionEvent()
		else:
			raise ValueError("Unknown DivisionEvent type"+str(event_type))


		self.add_event(event_object, event_params, prop_object, propensity_params, event_type = "division", volume_splitter = volume_splitter)

	def create_volume_event(self, event_type, dict event_params, str event_propensity_type, dict propensity_params, print_out = False):
		event_params = dict(event_params)
		propensity_params = dict(propensity_params)
		if print_out:
			warnings.warn("Creating New Volume event\n\ttype="+event_type+"\n\tparams="+str(event_params)+"\n\tprop_type="+str(event_propensity_type)+"\n\tprop_params="+str(propensity_params))
		prop_object = self.create_propensity(event_propensity_type, propensity_params, print_out = print_out)

		if event_type.lower() in ["linear", "linear volume", "linearvolume" "linearvolumeevent", "linear volume event"]:
			self._param_dict_check(event_params, "growth_rate", "DummyVar_LinearVolumeEvent")
			event_object = LinearVolumeEvent()
		elif event_type.lower() in ["multiplicative", "multiplicative volume", "multiplicativevolume", "multiplicativevolumeevent", "multiplicative volume event"]:
			self._param_dict_check(event_params, "growth_rate", "DummyVar_MultiplicativeVolumeEvent")
			event_object = MultiplicativeVolumeEvent()
		elif event_type.lower() in ["general", "general volume", "generalvolume", "generalvolumeevent", "general volume event"]:
			event_object = GeneralVolumeEvent()
		else:
			raise ValueError("Unknown VolumeEvent Type: "+str(event_type))

		self.add_event(event_object, event_params, prop_object, propensity_params, event_type = "volume")

	def add_lineage_rule(self, LineageRule rule_object, dict rule_param_dict, str rule_type, VolumeSplitter volume_splitter = None):
		species_names, param_names = rule_object.get_species_and_parameters(rule_param_dict)

		for species_name in species_names:
			self._add_species(species_name)
		for param_name in param_names:
			self._add_param(param_name)
		if "division" in rule_type.lower():
			if volume_splitter == None:
				raise ValueError("DivisionRules must be added with a volume splitter object in add_lineage_rule.")
			rule_object.initialize(rule_param_dict, self.species2index,  self.params2index)
			self.division_rules_list.append((rule_object, volume_splitter))
		else:
			rule_object.initialize(rule_param_dict, self.species2index,  self.params2index)
			if "death" in rule_type.lower():
				self.death_rules.append(rule_object)
			elif "volume" in rule_type.lower():
				self.volume_rules.append(rule_object)
			else:
				raise ValueError("add_lineage_rule only takes rules of type 'DeathRule', 'DivisionRule', and 'VolumeRule'. For Other rule types, consider trying Model.add_rule.")

	def create_death_rule(self, str rule_type, dict rule_param_dict):
		if rule_type.lower() in ["species", "speciesdeathrule"]:
			self._param_dict_check(rule_param_dict, "specie", "DummyVar_SpeciesDeathRule")
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_SpeciesDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_SpeciesDeathRule")
			rule_object = SpeciesDeathRule()
		elif rule_type.lower() in ["param", "parameter", "paramdeathrule"]:
			self._param_dict_check(rule_param_dict, "param", "DummyVar_ParamDeathRule")
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_ParamDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_ParamDeathRule")
			rule_object = ParamDeathRule()
		elif rule_type.lower() in ["general" "generaldeathrule"]:
			rule_object = GeneralDeathRule()
		else:
			raise ValueError("Unknown DeathRule type: "+str(rule_type))

		self.add_lineage_rule(rule_object, rule_param_dict, rule_type = "death")

	def create_division_rule(self, str rule_type, dict rule_param_dict, VolumeSplitter volume_splitter):
		if rule_type.lower() in ["time", "timedivisionrule"]:
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_TimeDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_TimeDeathRule")
			rule_object = TimeDivisionRule()
		elif rule_type.lower() in ["volume", "volumedivisionrule"]:
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_VolumeDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_VolumeDeathRule")
			rule_object = VolumeDivisionRule()
		elif rule_type.lower() in ["delta", "deltav", "deltavdivisionrule"]:
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_DeltaVDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_DeltaVDeathRule")
			rule_object = DeltaVDivisionRule()
		elif rule_type.lower() in ["general", "generaldivisionrule"]:
			rule_object = GeneralDivisionRule()
		else:
			raise ValueError("Unknown DivisionRule type: "+str(rule_type))
		self.add_lineage_rule(rule_object, rule_param_dict, rule_type = 'division', volume_splitter = volume_splitter)

	def create_volume_rule(self, str rule_type, dict rule_param_dict):
		if rule_type.lower() in ["linear", "linearvolumerule"]:
			self._param_dict_check(rule_param_dict, "growth_rate", "DummyVar_LinearVolumeRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_LinearVolumeRule")
			rule_object = LinearVolumeRule()
		elif rule_type.lower() in ["multiplicative", "multiplicativevolume", "multiplicativevolumerule"]:
			self._param_dict_check(rule_param_dict, "growth_rate", "DummyVar_MultiplicativeVolumeRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_MultiplicativeVolumeRule")
			rule_object = MultiplicativeVolumeRule()
		elif rule_type.lower() in ["assignment", "assignmentvolumerule"]:
			rule_object = AssignmentVolumeRule()
		elif rule_type.lower() in ["ode", "odevolumerule"]:
			rule_object = ODEVolumeRule()
		else:
			raise ValueError("Unknown VolumeRule type: "+str(rule_type))
		self.add_lineage_rule(rule_object, rule_param_dict, rule_type = 'volume')

	cdef unsigned get_num_division_rules(self):
		return self.num_division_rules
	cdef unsigned get_num_volume_rules(self):
		return self.num_volume_rules
	cdef unsigned get_num_death_rules(self):
		return self.num_death_rules
	cdef unsigned get_num_division_events(self):
		return self.num_division_events
	cdef unsigned get_num_volume_events(self):
		return self.num_volume_events
	cdef unsigned get_num_death_events(self):
		return self.num_death_events
	cdef list get_lineage_propensities(self):
		return self.lineage_propensities
	def py_get_lineage_propensities(self):
		return self.lineage_propensities
	def py_get_num_lineage_propensities(self):
		return len(self.lineage_propensities)
	def py_get_num_division_rules(self):
		return self.get_num_division_rules()
	def py_get_num_volume_rules(self):
		return self.get_num_volume_rules()
	def py_get_num_death_rules(self):
		return self.get_num_death_rules()
	def py_get_num_division_events(self):
		return self.get_num_division_events()
	def py_get_num_volume_events(self):
		return self.get_num_volume_events()
	def py_get_num_death_events(self):
		return self.get_num_death_events()

	cdef (vector[void*])* get_c_lineage_propensities(self):
		return & self.c_lineage_propensities
	cdef (vector[void*])* get_c_division_rules(self):
		return & self.c_division_rules
	cdef (vector[void*])* get_c_volume_rules(self):
		return & self.c_volume_rules
	cdef (vector[void*])* get_c_death_rules(self):
		return & self.c_death_rules
	cdef (vector[void*])* get_c_division_events(self):
		return & self.c_division_events
	cdef (vector[void*])* get_c_volume_events(self):
		return & self.c_volume_events
	cdef (vector[void*])* get_c_death_events(self):
		return & self.c_death_events

	def py_get_volume_splitters(self):
		return self.rule_volume_splitters, self.event_volume_splitters

cdef class LineageCSimInterface(ModelCSimInterface):
	cdef unsigned num_division_events
	cdef unsigned num_division_rules
	cdef unsigned num_death_events
	cdef unsigned num_death_rules
	cdef unsigned num_volume_events
	cdef unsigned num_volume_rules
	cdef unsigned num_lineage_propensities

	cdef list volume_events_list
	cdef list division_events_list
	cdef list death_events_list

	cdef vector[void*] *c_lineage_propensities
	cdef vector[void*] *c_death_events
	cdef vector[void*] *c_division_events
	cdef vector[void*] *c_volume_events
	cdef vector[void*] *c_other_events
	cdef vector[void*] *c_death_rules
	cdef vector[void*] *c_division_rules
	cdef vector[void*] *c_volume_rules

	cdef list division_event_volume_splitters
	cdef list division_rule_volume_splitters

	def __init__(self, LineageModel M):
		super().__init__(M)
		self.num_division_rules = <unsigned>M.py_get_num_division_rules()
		self.num_volume_rules = <unsigned>M.py_get_num_volume_rules()
		self.num_death_rules = <unsigned>M.py_get_num_death_rules()
		self.num_division_events = <unsigned>M.py_get_num_division_events()
		self.num_volume_events = <unsigned>M.py_get_num_volume_events()
		self.num_death_events = <unsigned>M.py_get_num_death_events()
		self.num_lineage_propensities = <unsigned>M.py_get_num_lineage_propensities()
		self.c_lineage_propensities = M.get_c_lineage_propensities()
		self.c_volume_rules = M.get_c_volume_rules()
		self.c_division_rules = M.get_c_division_rules()
		self.c_death_rules = M.get_c_death_rules()
		self.c_volume_events = M.get_c_volume_events()
		self.c_division_events = M.get_c_division_events()
		self.c_death_events = M.get_c_death_events()
		self.division_rule_volume_splitters, self.division_event_volume_splitters = M.py_get_volume_splitters()


	#similar to compute_stochastic_volume_propensities but events are new included as well
	cdef void compute_lineage_propensities(self, double* state, double* propensity_destination, double volume, double time):
		cdef unsigned ind
		for ind in range(self.num_reactions):
			propensity_destination[ind] = (<Propensity>(self.c_propensities[0][ind])).get_stochastic_volume_propensity(&state[0], self.c_param_values, volume, time)

		for ind in range(self.num_lineage_propensities):
			propensity_destination[self.num_reactions+ind] = (<Propensity>(self.c_lineage_propensities[0][ind])).get_stochastic_volume_propensity(&state[0], self.c_param_values, volume, time)

	#Applies all rules that change the volume of a cell
	cdef double apply_volume_rules(self, double* state, double volume, double time, double dt):
		cdef int ind
		#print("num_volume_rules", self.num_volume_rules, "state", state[0], "volume", volume, "time", time, "dt", dt)
		for ind in range(self.num_volume_rules):
			#print("ind", ind)
			#print("(<VolumeRule>self.c_volume_rules[0][ind])", (<VolumeRule>self.c_volume_rules[0][ind]))
			volume = (<VolumeRule>self.c_volume_rules[0][ind]).get_volume(state, self.c_param_values, volume, time, dt)
		return volume

	#Applies death rules in the order they were added to the model. Returns the index of the first death rule that returns True. -1 otherwise.
	cdef int apply_death_rules(self, double* state, double volume, double time, double start_volume, double start_time):
		cdef int isdead = 0
		cdef int ind
		for ind in range(self.num_death_rules):
			isdead = (<DeathRule>self.c_death_rules[0][ind]).check_dead(state, self.c_param_values, time, volume, start_time, start_volume)
			if isdead > 0:
				return ind
		return -1

	#Applies divison rules in the order they were added to the model. Returns the index of the first division rule that returns True. -1 otherwise
	cdef int apply_division_rules(self, double* state, double volume, double time, double start_volume, double start_time):
		cdef int divided = 0
		cdef int ind
		#print("state", state[0], "volume", volume, "time", time, "start_volume", start_volume, "start_time", start_time)
		for ind in range(self.num_division_rules):
			#print("ind", ind)
			#print("division rule:", <DivisionRule>self.c_division_rules[0][ind])
			divided = (<DivisionRule>self.c_division_rules[0][ind]).check_divide(state, self.c_param_values, time, volume, start_time, start_volume)
			if divided > 0:
				return ind
		return -1

	#Applies a single volume event, determined by the index passed in
	cdef double apply_volume_event(self, int event_index, double* state, double current_time, double current_volume):
		#print("event_index", event_index, "state", state[0], "current_time", current_time, "current_vol", current_volume)
		#print("<VolumeEvent>self.c_volume_events[0][event_index]")
		current_volume = (<VolumeEvent>self.c_volume_events[0][event_index]).get_volume(state, self.c_param_values, current_volume, current_time)
		return current_volume

	#Divides a single cell using a VolumeSplitter determined by vsplit_index
	cdef np.ndarray partition(self, int vsplit_ind, LineageVolumeCellState parent):
		cdef VolumeSplitter vsplit
		if vsplit_ind >= self.num_division_rules and vsplit_ind<self.num_division_rules+self.num_division_events:
			vsplit_ind = vsplit_ind - self.num_division_rules
			vsplit = self.division_event_volume_splitters[vsplit_ind]
		elif vsplit_ind < self.num_division_rules and vsplit_ind >= 0:
			vsplit = self.division_rule_volume_splitters[vsplit_ind]
		else:
			raise ValueError('Invalid volume splitter index: vsplit_ind='+str(vsplit_ind))

		return vsplit.partition(parent)


	cdef np.ndarray delay_partition(self, unsigned vsplit_ind, LineageVolumeCellState parent):
		raise NotImplementedError("Implement me!")

	cdef unsigned get_num_lineage_propensities(self):
		return self.num_lineage_propensities

	cdef unsigned get_num_division_events(self):
		return self.num_division_events
	cdef unsigned get_num_volume_events(self):
		return self.num_volume_events
	cdef unsigned get_num_death_events(self):
		return self.num_death_events
	cdef unsigned get_num_volume_rules(self):
		return self.num_volume_rules
	cdef unsigned get_num_death_rules(self):
		return self.num_death_rules
	cdef unsigned get_num_division_rules(self):
		return self.num_division_rules

#A new wrapper for the VolumeCellState with new internal variables
cdef class LineageVolumeCellState(DelayVolumeCellState):
	cdef double initial_volume #Stores the birth Volume
	cdef double initial_time #Stores the time the Cell was "born"
	#divided = -1: Cell Not divided
	#divided E [0, num_division_rules): DivisionRule divided caused the cell to divide
	#divided E [num_division_rules, num_division_rules + num_division_events]: Division Event divided-num_division_rules caused the cell to divide
	cdef int divided
	#dead = -1: Cell Not dead
	#dead E [0, num_death_rules): DeathRule divided caused the cell to die
	#dead E [num_death_rules, num_death_rules + num_death_events]: DeathEvent dead-num_death_rules caused the cell to die
	cdef int dead
	cdef state_set

	def __init__(self, v0 = 0, t0 = 0, state = []):
		self.set_initial_vars(v0, t0)
		self.set_volume(v0)
		self.set_time(t0)
		self.volume_object = None
		self.divided = -1
		self.dead = -1
		if len(state) == 0:
			self.state_set = 0
		else:
			self.state_set = 1
			self.py_set_state(state)

	def get_state_set(self):
		return self.state_set

	def py_set_state(self, state):
		self.state_set = 1
		return super().py_set_state(np.asarray(state))

	cdef void set_initial_vars(self, double volume, double time):
		self.initial_volume = volume
		self.initial_time = time

	cdef double get_initial_volume(self):
		return self.initial_volume

	cdef void set_state_comp(self, double val, unsigned comp_ind):
		self.state[comp_ind] = val

	cdef double get_state_comp(self, unsigned comp_ind):
		return self.state[comp_ind]

	def py_get_initial_volume(self):
		return self.get_initial_volume()

	cdef double get_initial_time(self):
		return self.initial_time

	def py_get_initial_time(self):
		return self.get_initial_time()

	cdef void set_divided(self, divided):
		self.divided = divided

	cdef int get_divided(self):
		return self.divided

	cdef void set_dead(self, dead):
		self.dead = dead
	cdef int get_dead(self):
		return self.dead

cdef class SingleCellSSAResult(VolumeSSAResult):
	#divided = -1: Cell Not divided
	#divided E [0, num_division_rules): DivisionRule divided caused the cell to divide
	#divided E [num_division_rules, num_division_rules + num_division_events]: Division Event divided-num_division_rules caused the cell to divide
	cdef int divided
	#dead = -1: Cell Not dead
	#dead E [0, num_death_rules): DeathRule divided caused the cell to die
	#dead E [num_death_rules, num_death_rules + num_death_events]: DeathEvent dead-num_death_rules caused the cell to die
	cdef int dead

	#save the initial time and volume
	#typically not used - but occasionally useful for Propagating cells and not saving trajectories
	cdef double t0
	cdef double v0

	cdef void set_dead(self, int dead):
		self.dead = dead
	def py_set_dead(self, dead):
		self.set_dead(dead)

	cdef int get_dead(self):
		return self.dead
	def py_get_dead(self):
		return self.get_dead()

	cdef void set_divided(self, int divided):
		self.divided = divided
	def py_set_divided(self, divided):
		self.set_divided(divided)

	cdef int get_divided(self):
		return self.divided
	def py_get_divided(self):
		return self.get_divided()

	cdef void set_initial_time(self, double t):
		self.t0 = t
	def py_set_initial_time(self, t):
		self.set_initial_time(t)

	cdef void set_initial_volume(self, double v):
		self.v0 = v
	def py_set_initial_volume(self, v):
		self.set_initial_volume(v)

	cdef VolumeCellState get_final_cell_state(self):
		#cdef unsigned final_index = (<np.ndarray[np.double_t,ndim=1]> self.timepoints).shape[0]-1
		cdef unsigned final_index = self.timepoints.shape[0]-1


		if self.t0 is None:
			self.set_initial_time(self.timepoints[0])
		if self.v0 is None:
			self.set_initial_volume(self.volume[0])


		cdef LineageVolumeCellState cs  = LineageVolumeCellState(t0 = self.t0, v0 = self.v0, state = self.simulation_result[final_index,:])
		cs.set_time(self.timepoints[final_index])
		cs.set_volume(self.volume[final_index])
		cs.set_divided(self.divided)
		cs.set_dead(self.dead)
		return cs


cdef class LineageVolumeSplitter(VolumeSplitter):
	cdef unsigned how_to_split_v
	cdef unsigned how_to_split_default
	cdef vector[int] binomial_indices
	cdef vector[int] perfect_indices
	cdef vector[int] duplicate_indices
	cdef vector[int] custom_indices
	cdef double partition_noise
	cdef dict ind2customsplitter
	cdef dict custom_partition_functions

	def __init__(self, Model M, options = {}, custom_partition_functions = {}, partition_noise = .5):
		self.ind2customsplitter == {}
		self.custom_partition_functions = custom_partition_functions
		if self.partition_noise > 1:
			raise ValueError("Partition Noise must be between 0 and 1")
		self.partition_noise = partition_noise

		#Check if the default has been changed
		if "default" not in options or options["default"] == "binomial":
			default = "binomial"
		elif options["default"] == "duplicate":
			default = "duplicate"
		elif options["default"] == "perfect":
			default = "perfect"
		elif options["default"] in custom_partition_functions:
			default = "custom"
			self.ind2customsplitter["default"] = options["default"]
		else:
			raise ValueError("Custom partition function key, "+str(options["default"])+", for 'default' not in custom_partition_functions")

		#Figure out how volume will be split
		if ("volume" not in options and default == "binomial") or options["volume"] == "binomial":
			self.how_to_split_v = 0
		elif ("volume" not in options and default == "duplicate") or options["volume"] == "duplicate":
			self.how_to_split_v = 1
		elif ("volume" not in options and default == "perfect") or options["volume"] == "perfect":
			self.how_to_split_v = 2
		elif ("volume" not in options and default == "custom") or options["volume"] in custom_partition_functions:
			self.how_to_split_v = 3
			if "volume" in options:
				self.ind2customsplitter["volume"] = options["volume"]
			else:
				self.ind2customsplitter["volume"] = options["default"]
		else:
			raise ValueError("Custom partition function key, "+str(options["volume"])+", for 'volume' not in custom_partition_functions")

		#Figure out how other species are split
		for s in M.get_species2index():
			index = M.get_species_index(s)
			if s not in options:
				if default == "binomial":
					self.binomial_indices.push_back(index)
				elif default == "duplicate":
					self.duplicate_indices.push_back(index)
				elif default == "perfect":
					self.perfect_indices.push_back(index)
				elif default == "custom":
					self.custom_indices.push_back(index)
					self.ind2customsplitter[index] = options["default"]
			elif options[s] == "binomial":
				self.binomial_indices.push_back(index)
			elif options[s] == "duplicate":
				self.duplicate_indices.push_back(index)
			elif options[s] == "perfect":
				self.perfect_indices.push_back(index)
			elif options[s] in custom_partition_functions:
				self.custom_indices.push_back(index)
				self.ind2customsplitter[index] = options[s]
			else:
				raise ValueError("Custom partition function key, "+str(options["volume"])+", for "+s+" not in custom_partition_functions")

	cdef np.ndarray partition(self, VolumeCellState parent):
		cdef double v0d, v0e, t0, p, q, c1, c2
		# set times
		t0 = parent.get_time()

		# partition the states, copying already takes care of duplication replications.
		cdef np.ndarray[np.double_t,ndim=1] dstate = parent.get_state().copy()
		cdef np.ndarray[np.double_t,ndim=1] estate = parent.get_state().copy()
		cdef unsigned length = dstate.shape[0]

		cdef unsigned loop_index = 0
		cdef unsigned species_index = 0
		cdef unsigned amount = 0
		cdef unsigned amount2 = 0
		cdef double d_value = 0.0

		# simulate partitioning noise
		if self.how_to_split_v == 0: #Binomial
			p = 0.5 - cyrandom.uniform_rv()*self.partition_noise/2.
			q = 1 - p
			# calculate binomial volumes
			v0d = parent.get_volume() * p
			v0e = parent.get_volume() * q
		elif self.how_to_split_v == 1: #Duplicate
			p = 1
			q = 1
			v0d = parent.get_volume()
			v0e = parent.get_volume()
		elif self.how_to_split_v == 2: #perfect
			p = .5
			q = .5
			v0d = parent.get_volume()*.5
			v0e = parent.get_volume()*.5
		else:
			splitter = self.custom_partition_functions[self.ind2customsplitter["volume"]]
			v0d, v0e = splitter(parent, "volume")
			p = v0d/parent.get_volume()
			q = v0e/parent.get_volume()
			if v0d <= 0 or v0e <= 0:
				raise ValueError("splitter "+self.ind2customsplitter[species_index]+" returned negative quantities for volume")

		# take care of perfect splitting
		for loop_index in range(self.perfect_indices.size()):
			species_index = self.perfect_indices[loop_index]
			d_value = p * dstate[species_index]
			amount = <int> (d_value+0.5)
			if d_value-amount <= 1E-8 and amount>0:
				dstate[species_index] = <double> amount
			elif amount <0:
				raise ValueError('negative quantity in perfect partitioning')
			else:
				if cyrandom.uniform_rv() <= p:
					dstate[species_index] = <int> d_value + 1
				else:
					dstate[species_index] = <int> d_value
			estate[species_index] -= dstate[species_index]

		# take care of binomial splitting
		for loop_index in range(self.binomial_indices.size()):
			species_index = self.binomial_indices[loop_index]
			amount = cyrandom.binom_rnd_f(dstate[species_index],p)
			dstate[species_index] = <double> amount
			estate[species_index] -= dstate[species_index]

		for loop_index in range(self.custom_indices.size()):
			species_index = self.custom_indices[loop_index]
			splitter = self.custom_partition_functions[self.ind2customsplitter[species_index]]
			c1, c2 = splitter(species_index, parent)
			if c1 < 0 or c2 < 0:
				raise ValueError("splitter "+self.ind2customsplitter[species_index]+" returned negative quantities for species index "+str(species_index))
			dstate[species_index] = <double> c1
			estate[species_index] = <double> c2

		# create return structure
		cdef np.ndarray ans = np.empty(2, dtype=np.object)
		cdef LineageVolumeCellState d = LineageVolumeCellState(v0 = v0d, t0 = parent.get_time(), state = dstate)
		cdef LineageVolumeCellState e = LineageVolumeCellState(v0 = v0e, t0 = parent.get_time(), state = estate)
		ans[0] = d
		ans[1] = e
		return ans

cdef class LineageSSASimulator:

	#Memory Views are reused between each individual cell
	#Sometimes, to be compatable with the double* used in original bioscrape, these are cast to double*
	#these are used by SimulateSingleCell
	cdef double[:] c_timepoints, c_current_state, c_propensity, c_truncated_timepoints, c_volume_trace
	cdef double[:, :] c_stoich, c_results

	#An Interface is stored ina  linear model for fast helper functions
	cdef LineageCSimInterface interface

	#All the counters are also reused for each individual cell and set when interface is set
	cdef unsigned num_species, num_reactions, num_volume_events, num_death_events, num_division_events, num_volume_rules, num_death_rules, num_division_rules, num_propensities, num_timepoints

	#These are used by SimulateLineage and PropagateCells
	cdef Schnitz s, daughter_schnitz1, daughter_schnitz2
	cdef LineageVolumeCellState d1, d2, d1final, d2final, cs, dcs
	cdef SingleCellSSAResult r
	cdef list old_schnitzes, old_cell_states
	cdef Lineage lineage

	#Used for PropagateCells
	cdef list cell_states

	#Used to create a propensity buffer from an interface
	cdef create_propensity_buffer(self, LineageCSimInterface interface):
		cdef np.ndarray[np.double_t, ndim = 1] c_propensity = np.zeros(interface.get_num_lineage_propensities()+interface.get_num_reactions())
		return c_propensity

	#Python accessible version
	def py_create_propensity_buffer(self, LineageCSimInterface interface):
		return self.create_propensity_buffer(interface)

	#helper function to take an np.ndarrays and set them to internal memory views
	cdef void set_c_truncated_timepoints(self, np.ndarray timepoints):
		self.c_truncated_timepoints = timepoints
	#python accessible version
	def py_set_c_truncated_timepoints(self, np.ndarray timepoints):
		self.set_c_truncated_timepoints(timepoints)

	cdef void set_c_timepoints(self, np.ndarray timepoints):
		self.c_timepoints = timepoints
	#Python accessible version
	def py_set_c_timepoints(self, np.ndarray timepoints):
		self.set_c_timepoints(timepoints)

	#Sets the internal interface and associated internal variables
	#Main speed-up due to not having to set c_stoic (which could be very large) as often
	cdef void intialize_single_cell_interface(self, LineageCSimInterface interface):
		#Reset internal variables
		self.interface = interface

		#Memory View Setup
		#Stochiomettric Matrix
		self.c_stoich = self.interface.get_update_array() + self.interface.get_delay_update_array()

		self.num_species = self.interface.get_number_of_species()
		self.num_reactions = self.interface.get_number_of_reactions()

		self.num_volume_events = self.interface.get_num_volume_events()
		self.num_death_events = self.interface.get_num_death_events()
		self.num_division_events = self.interface.get_num_division_events()
		self.num_volume_rules = self.interface.get_num_volume_rules()
		self.num_death_rules = self.interface.get_num_death_rules()
		self.num_division_rules = self.interface.get_num_division_rules()
		self.num_propensities = self.num_reactions + self.num_volume_events + self.num_death_events + self.num_division_events

		#print("speces", self.num_species, "reactions",self.num_reactions)
		#print("volume events", self.num_volume_events, "death events", self.num_death_events, "division evets", self.num_division_events)
		#print("volume_rules", self.num_volume_rules, "death rules", self.num_death_rules, "division_rules", self.num_division_rules)
		#print("total propensities", self.num_propensities)
		#Prepare propensity buffer of the right size
		self.c_propensity = np.zeros(self.num_propensities)
		self.c_current_state = np.zeros(self.num_species)

	def py_initialize_single_cell_interface(self, LineageCSimInterface interface):
		self.initialize_single_cell_interface(interface)

	#Allocates buffer to store results from SimulateSingleCell
	cdef void initialize_single_cell_results_arrays(self, unsigned num_timepoints):
		cdef np.ndarray[np.double_t,ndim=2] results 
		cdef np.ndarray[np.double_t,ndim=1] volume_trace

		if num_timepoints == 0:
			num_timepoints = 1 
		results = np.zeros((num_timepoints,self.num_species))
		self.c_results = results
		volume_trace = np.zeros(num_timepoints,)
		self.c_volume_trace = volume_trace

	#python accessible version
	def py_initialize_single_cell_results_arrays(self, int num_timepoints):
		self.initialize_single_cell_results_arrays(num_timepoints)

	#SSA for a single cell. Simulates until it devides or dies using division / death rules and/or reactions.
	#Before calling this for a given interface, must call initialize_single_cell_interface to set up internal variables
	#Python wrapper below takes care of all the details, but will be slower if used repeatedly
	#unsigned mode: if 1 returns a full SingleCellSSAResult. If 0 returns a dummy SingleCellSSAResult of length 1
	cdef SingleCellSSAResult SimulateSingleCell(self, LineageVolumeCellState v, double[:] timepoints, unsigned mode):
		#print("SimulateSingleCell", "mode=", mode)

		#Memory views are reused from other objects for less allocation
		cdef unsigned num_timepoints = len(timepoints)

		cdef double initial_time = v.get_initial_time()
		cdef double current_time = v.get_time()

		cdef double final_time = timepoints[num_timepoints-1]
		cdef double proposed_time = 0.0
		cdef unsigned current_index = 0
		cdef unsigned reaction_choice = 4294967295 # https://en.wikipedia.org/wiki/4,294,967,295
		cdef unsigned species_index = 4294967295
		cdef double delta_t = timepoints[1]-timepoints[0]
		cdef double next_queue_time = timepoints[current_index+1]
		cdef double move_to_queued_time = 0
		cdef double initial_volume = v.get_initial_volume()
		cdef double current_volume = v.get_volume()
		cdef int cell_divided = -1
		cdef int cell_dead = -1
		cdef double Lambda = 0

		#print("single cell simulation starting at ", current_time, "till", final_time)

		#1 dimensional ndarrays used for mode=0 simulation only when fully results aren't saved
		cdef np.ndarray[np.double_t,ndim=1] dummy_t
		cdef np.ndarray[np.double_t,ndim=1] dummy_v
		cdef np.ndarray[np.double_t,ndim=2] dummy_r

		#These are kept as local numpy arrays because they are returned after every simulation
		cdef SingleCellSSAResult SCR

		#If Mode is 1 we will need to allocate new memory for this simulation
		if mode == 1:
			self.initialize_single_cell_results_arrays(num_timepoints)

		elif mode == 0:
			dummy_t = np.zeros(1)
			dummy_v = np.zeros(1)
			dummy_r = np.zeros((1, self.num_species))

		elif self.c_results is None or self.c_volume_trace is None:
			raise RuntimeError("Must call LineageSSASimulator.initialize_single_cell_results_arrays(num_timepoints) before calling SimulateSingleCell. Setting the mode=1 will do this automatically for you.")

		#Set Initial State
		if v.get_state_set() == 1:
			self.c_current_state = v.py_get_state().copy()
		else:
			warnings.warn("No initial state set (via LineageVolumeCellState v) in SingleCellSSAResuslt. Defaulting to the Model's initial state.")
			self.c_current_state = self.interface.get_initial_state().copy()
			v.py_set_state(self.c_current_state)

		#Warn user if delays are in the model (which will be converted to non-delay reactions)
		if (self.interface.py_get_delay_update_array() != np.zeros(self.interface.py_get_delay_update_array().shape)).any():
			warnings.warn("Delay reactions found in the model. SingleCellSSASimulator will simulate these reactions without delay. Delays are not yet supported for LineageModels but can be simulated as regular Models with the DelayVolumeSSASimulator.")

		# Do the SSA part now
		#print("SingleCell loop start", self.interface, v)
		while current_index < num_timepoints:
			#print("a1")
			# Compute rules in place
			self.interface.apply_repeated_volume_rules(&self.c_current_state[0], current_volume, current_time)
			#print("a2")
			#returns the index of the first DeathRule that returned True and -1 otherwise
			cell_dead = self.interface.apply_death_rules(&self.c_current_state[0], current_volume, current_time, initial_volume, initial_time)
			#print("a3", len(self.c_current_state))
			#print("self.c_current_state", self.c_current_state[0], self.c_current_state[1], self.c_current_state[2])
			#print("current_time", current_time, "initial_volume", initial_volume, "initial_time", initial_time)
			#returns the index of the first DivisionRule that returned True and -1 otherwise
			cell_divided = self.interface.apply_division_rules(&self.c_current_state[0], current_volume, current_time, initial_volume, initial_time)

			#print("A", "prop len", len(self.c_propensity))
			#Break the loop cell dead or divided
			if cell_dead >= 0 and cell_divided >= 0:
				warnings.warn("Cell Death and Division Occured Simultaneously - Death Takes Precedent")
				cell_divided = -1
				break
			elif cell_dead >= 0:
				break
			elif cell_divided >= 0:
				break
			#Compute Reaction and Event propensities in-place
			self.interface.compute_lineage_propensities(&self.c_current_state[0], &self.c_propensity[0], current_volume, current_time)
			#print("b1 self.num_propensities", self.num_propensities)
			Lambda = cyrandom.array_sum(&self.c_propensity[0], self.num_propensities)
			# Either we are going to move to the next queued time, or we move to the next reaction time.
			#print("B")
			if Lambda == 0:
				proposed_time = final_time+1
			else:
				proposed_time = current_time + cyrandom.exponential_rv(Lambda)
			if next_queue_time < proposed_time and next_queue_time < final_time:
				current_time = next_queue_time
				next_queue_time += delta_t
				move_to_queued_time = 1
			elif proposed_time > final_time-10e-8:
				current_time = final_time
				move_to_queued_time = 1
			else:
				current_time = proposed_time
				move_to_queued_time = 0

			v.set_time(current_time)

			# Update the results array with the state for the time period that we just jumped through.
			while current_index < num_timepoints and timepoints[current_index] <= current_time:
				for species_index in range(self.num_species):
					self.c_results[current_index, species_index] = self.c_current_state[species_index]
				self.c_volume_trace[current_index] = current_volume
				current_index += 1
			#print("C")
			# Now update the state accordingly.
			# IF the queue won, then update the volume and continue on or stop if the cell divided.
			if move_to_queued_time == 1:
				# Update the volume every dtyp
				current_volume = self.interface.apply_volume_rules(&self.c_current_state[0], current_volume, current_time, delta_t)
				v.set_volume(current_volume)

			# if an actual reaction happened, do the reaction and maybe update the queue as well.
			else:
				#print("D")
				# select a reaction
				reaction_choice = cyrandom.sample_discrete(self.num_propensities, &self.c_propensity[0], Lambda )
				#Propensities are Ordered:
				# Reactions, Divison Events, Volume Events, Death Events

				#Propensity is a reaction
				#print("d0", "reaction_choice", reaction_choice)
				if reaction_choice < self.num_reactions:

					#print("d1", "reaction_choice", reaction_choice, "self.num_reactions", self.num_reactions, "self.num_species", self.num_species)

					# Do the reaction's initial stoichiometry.
					for species_index in range(self.num_species):
						self.c_current_state[species_index] += self.c_stoich[species_index, reaction_choice]

				#Propensity is a VolumeEvent
				elif reaction_choice >= self.num_reactions and reaction_choice < self.num_reactions + self.num_volume_events:
					#print("d2", "num_reactions", self.num_reactions)
					#print("volume event! volume before", current_volume)
					current_volume = self.interface.apply_volume_event(reaction_choice - self.num_reactions, &self.c_current_state[0], current_time, current_volume)
					v.set_volume(current_volume)
					#print("volume after", v.get_volume())
				#Propensity is a DivisionEvent.
				elif reaction_choice >= self.num_reactions+self.num_volume_events and reaction_choice < self.num_reactions + self.num_volume_events+self.num_division_events:
					#print("d3")
					#Cell Divided = DivisionEvent Index + num_division_rules
					cell_divided = reaction_choice - self.num_reactions - self.num_volume_events + self.num_division_rules
					break
				#Propensity is a Death Event
				elif reaction_choice >= self.num_reactions + self.num_volume_events+self.num_division_events:
					#print("d4")
					#Cell Divided = DeathEvent Index + num_death_rules
					cell_dead = reaction_choice - self.num_reactions + self.num_volume_events+self.num_division_events+self.num_death_rules
					break
				else:
					raise ValueError("More reaction propensities than expected!")

		#print("Out of SingleCell loop")
		if cell_divided>=0 or cell_dead>=0:
			#Push current state to the nearest index
			if current_time < timepoints[current_index]:
				for species_index in range(self.num_species):
					self.c_results[current_index,species_index] = self.c_current_state[species_index]
				self.c_volume_trace[current_index] = current_volume
				current_index += 1

			if mode == 1:
				if current_index == 0: #Can't return an empty array!
					timepoints = timepoints[:current_index+1]
					self.c_volume_trace = self.c_volume_trace[:current_index+1]
					self.c_results = self.c_results[:current_index+1,:]
					timepoints = timepoints[:current_index+1]
				else:
					timepoints = timepoints[:current_index]
					self.c_volume_trace = self.c_volume_trace[:current_index]
					self.c_results = self.c_results[:current_index,:]
					timepoints = timepoints[:current_index]


		if current_time > final_time + 10e-8:
			warnings.warn(f"Current time ({current_time}) is greater than final time({final_time}) in SimulateSingleCell. You have possibly uncovered a bug!")

		#SCR (SingleCellSSAResult) contains the simulation results until cell death / division or simualtion termination.
		#cell_divided and cell_dead are returend via vsr so the events/rules/VolumeSplitters can be called by the lineage simualtion loop.
		if mode == 1: #default behavior to return the entire cell trajectory
			#print("SRC Instantiated in SimulateSingleCell")
			SCR = SingleCellSSAResult(np.asarray(timepoints), np.asarray(self.c_results), np.asarray(self.c_volume_trace), <unsigned>(cell_divided >= 0))
			SCR.set_divided(cell_divided)
			SCR.set_dead(cell_dead)
			SCR.set_volume_object(v.get_volume_object())
			SCR.set_initial_volume(self.c_volume_trace[0])
			SCR.set_initial_time(timepoints[0])
			#print("SRC creation complete")
		elif mode == 0: #memoryviemory saving behvior where a dummy SCR is returned with just one final state in it.
			#print("Single Cell Simulation Complete with current time, current_volume = ", current_time, current_volume)
			dummy_t[0] = current_time
			dummy_v[0] = current_volume
			for species_index in range(self.num_species):
				dummy_r[0, species_index] = self.c_current_state[species_index]
			SCR = SingleCellSSAResult(dummy_t, dummy_r, dummy_v, <unsigned>(cell_divided >= 0))
			SCR.set_divided(cell_divided)
			SCR.set_dead(cell_dead)
			SCR.set_initial_volume(self.c_volume_trace[0])
			SCR.set_initial_time(timepoints[0])


		#print("Simulate Single Cell Complete")
		return SCR

	def py_SimulateSingleCell(self, np.ndarray timepoints, LineageModel Model = None, LineageCSimInterface interface = None, LineageVolumeCellState v = None):

		if Model == None and interface == None:
			raise ValueError('py_SimulateSingleCell requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
		elif interface == None:
			interface = LineageCSimInterface(Model)
			interface.py_set_initial_time(timepoints[0])
		if v == None:
			v = LineageVolumeCellState(v0 = 1, t0 = 0, state = Model.get_species_array())

		self.set_c_truncated_timepoints(timepoints)
		self.intialize_single_cell_interface(interface)

		return self.SimulateSingleCell(v, self.c_truncated_timepoints, 1)#The 1 indicates that trajectory results will be saved to memory



	#Functions to simulate linages of cells
	#returns a truncated version of the memoryview array[x:]
	#starting at the first index x such that array[x] >= value
	cdef double[:] truncate_timepoints_less_than(self, double[:] array, double value):
		cdef unsigned j = 0
		for j in range(array.shape[0]):
			if array[j] >= value:
				return array[j:]
		return array #If nothing is less than the value, return the entire array

	#Functions to simulate linages of cells
	#returns a truncated version of the memoryview array[:x]
	#starting at the first index x such that array[x] > value
	cdef double[:] truncate_timepoints_greater_than(self, double[:] array, double value):
		cdef unsigned j = 0
		cdef int return_ind = -1
		for j in range(array.shape[0]):
			if array[j] > value:
				return_ind = j
				break

		if return_ind == -1: #Do not truncate anything
			return array
		else:
			return array[:return_ind]

	#Starts the simulation
	#add_to_lineage toggled for lineage versus propogation simulation
	cdef void simulate_cell_list(self, list initial_cell_states, double[:] timepoints, unsigned add_to_lineage, unsigned create_schnitzes):
		#print("Simulating cells N=", len(initial_cell_states))
		cdef unsigned i = 0

		cdef unsigned mode = 1
		if add_to_lineage == 0 and create_schnitzes == 0:
			mode = 0 #If this is 0, will run with less memory and not save trajectories

		for i in range(len(initial_cell_states)):
			self.r = self.SimulateSingleCell(initial_cell_states[i], timepoints, mode)

			self.old_cell_states.append(self.r.get_final_cell_state())

			if create_schnitzes or add_to_lineage:
				self.s = self.r.get_schnitz()
				self.s.set_parent(None)

				if create_schnitzes:
					self.old_schnitzes.append(self.s)

				if add_to_lineage:
					self.lineage.add_schnitz(self.s)



	#Simulate inside the Simulation Queue, 1 cell at a time
	#add_to_lineage toggled for lineage versus propogation simulation
	cdef void simulate_daughter_cells(self, double[:] timepoints, unsigned add_to_lineage, unsigned create_schnitzes):
		#print("simulate_daughter_cells")
		cdef double final_time = timepoints[timepoints.shape[0]-1]
		cdef unsigned mode = 1
		if add_to_lineage == 0 and create_schnitzes == 0:
			mode = 0 #If this is 0, will run with less memory and not save trajectories

		self.r = self.SimulateSingleCell(self.d1, timepoints, mode)
		self.d1final = self.r.get_final_cell_state()

		if add_to_lineage or create_schnitzes:
			self.daughter_schnitz1 = self.r.get_schnitz()

		# Add on the new daughter if final time wasn't reached.
		if self.d1final.get_time() < final_time + 1E-9:
			self.old_cell_states.append(self.d1final)

			if create_schnitzes:
				self.old_schnitzes.append(self.daughter_schnitz1)
		else:
			warnings.warn(f"Daughter cell simulation went over the total time. Simulation has been discarded. Check for model errors. daughter time = {self.d1final.get_time()} final time = {final_time}")

		self.r = self.SimulateSingleCell(self.d2, timepoints, mode)
		self.d2final = self.r.get_final_cell_state()

		if add_to_lineage or create_schnitzes:
			self.daughter_schnitz2 = self.r.get_schnitz()

		if self.d2final.get_time() < final_time + 1E-9:
			self.old_cell_states.append(self.d2final)

			if create_schnitzes:
				self.old_schnitzes.append(self.daughter_schnitz2)
		else:
			warnings.warn(f"Daughter cell simulation went over the total time. Simulation has been discarded. Check for model errors. daughter time = {self.d2final.get_time()} final time = {final_time}")

		if add_to_lineage or create_schnitzes:
			# Set up daughters and parent appropriately.
			self.daughter_schnitz1.set_parent(self.s)
			self.daughter_schnitz2.set_parent(self.s)
			self.s.set_daughters(self.daughter_schnitz1,self.daughter_schnitz2)

			if add_to_lineage:
				# Add daughters to the lineage
				self.lineage.add_schnitz(self.daughter_schnitz1)
				self.lineage.add_schnitz(self.daughter_schnitz2)


	#Simulates a lineage of cells keeptring track of mother-daughter relations
	cdef Lineage SimulateCellLineage(self, list initial_cell_states, double[:] timepoints):
		cdef unsigned i, j, list_index
		cdef int cell_divided, cell_dead
		cdef double final_time

		#Check instantation of core data structures
		#Note: these are not automatically instantiated here so they can be reused in more complex simulation types
		if self.lineage is None:
			raise RuntimeError("LineageSSASimulator.lineage must be instantiated to a Lineage before calling SimulateCellLineage. py_SimulateCellLineage automatically does this for you but is slower.")
		if self.old_cell_states is None:
			raise RuntimeError("LineageSSASimulator.old_cell_states must be instantiated to a list before calling SimulateCellLineage. py_SimulateCellLineage automatically does this for you but is slower.")
		if self.old_schnitzes is None:
			raise RuntimeError("LineageSSASimulator.old_schnitzes must be instantiated to a list before calling SimulateCellLineage. py_SimulateCellLineage automatically does this for you but is slower.")


		#These will be used to store outputs
		cdef np.ndarray daughter_cells

		#initialize variables
		final_time = timepoints[timepoints.shape[0]-1]
		i = 0
		list_index = 0
		cell_divided = -1
		cell_dead = -1

		# Simulate the first cell until death division or max time

		self.simulate_cell_list(initial_cell_states, timepoints, 1, 1) #toggle add_to_lineage = 1 create snitches 1

		while list_index < len(self.old_cell_states):
			self.cs = self.old_cell_states[list_index]
			self.s = self.old_schnitzes[list_index]

			list_index += 1

			#If the cell has already simulated all its time, do nothing
			if self.cs.get_time() >= final_time- 1E-9:
				pass
			#If the cell is dead, do nothing
			elif self.cs.get_dead() >= 0:
				pass
			#If the cell has divided, apply the appropriate division rule
			elif self.cs.get_divided() >= 0:

				#Check if dt is too small for accurate lineage simulation
				if self.cs.get_initial_time() == self.cs.get_time():
					raise ValueError("Cells are dividing too faster for the timepoints passed into SimulateCellLineage. Try decreasing the spacing between timepoints or limiting cell growth.")

				daughter_cells = self.interface.partition(self.cs.get_divided(), self.cs)

				self.d1 = <LineageVolumeCellState>(daughter_cells[0])
				self.d2 = <LineageVolumeCellState>(daughter_cells[1])

				#Create a new timepoint array and simulate the first daughter and queue if it doesn't reach final time.
				#self.c_truncated_timepoints = timepoints[timepoints >= self.cs.get_time()]
				self.c_truncated_timepoints = self.truncate_timepoints_less_than(timepoints, self.cs.get_time())
				self.simulate_daughter_cells(self.c_truncated_timepoints, 1, 1) #toggle add_to_lineage = 1 and add_schnitz = 1
			else:
				warnings.warn("SimulateLineage Should not reach this point. If it does, you have uncovered a bug!")

		return self.lineage

	#Python wrapper of the above
	def py_SimulateCellLineage(self, np.ndarray timepoints, initial_cell_states, LineageCSimInterface interface):
		#Instantiate variables
		self.lineage = Lineage()
		self.old_cell_states = []
		self.old_schnitzes = []
		self.set_c_timepoints(timepoints)
		self.intialize_single_cell_interface(interface)
		return self.SimulateCellLineage(initial_cell_states, self.c_timepoints)


	#Simulates an ensemble of cells over some amount of time.
	#Returns the existing cell states at all times in the array sample_times
	#  dead cells are included based upon the include_dead_cells parameter (1 = included, otherwise = excluded)
	cdef list PropagateCells(self, list initial_cell_states, double[:] timepoints, double[:] sample_times, unsigned include_dead_cells):
		#print("Propagating Cells")
		cdef unsigned list_index
		cdef double final_time
		cdef unsigned sample_ind
		cdef double dt = timepoints[1]-timepoints[0]

		samples = []
		if self.c_results is None or self.c_volume_trace is None:
			raise RuntimeError("Must call LineageSSASimulator.initialize_single_cell_results_arrays(num_timepoints) before calling PropagateCells. py_PropagateCells does this for you automatically but may be slower.")
		if self.old_cell_states is None:
			raise RuntimeError("LineageSSASimulator.old_cell_states must be instantiated to a list before calling PropagateCells. py_SimulateCellLineage automatically does this for you but is slower.")
		if self.cell_states is None:
			raise RuntimeError("LineageSSASimulator.cell_states must be instantiated to a list before calling PropagateCells. py_SimulateCellLineage automatically does this for you but is slower.")

		if len(timepoints) < sample_times.shape[0] or timepoints[timepoints.shape[0]-1] < sample_times[sample_times.shape[0]-1] or timepoints[0] > sample_times[0]:
			raise ValueError("sample_times must be a subset of the timepoints array.")

		#If the first sample is just the start of the simulation, append the initial condition
		if sample_times[0] == timepoints[0]:
			samples.append(initial_cell_states)
			sample_times = sample_times[1:]

		for sample_ind in range(sample_times.shape[0]):
			final_time = sample_times[sample_ind]
			list_index = 0

			self.c_timepoints = self.truncate_timepoints_greater_than(timepoints, final_time+1E-10)
			#Simulate the initial cells for the first sample
			if sample_ind == 0:
				#print("initial simulation")
				self.simulate_cell_list(initial_cell_states, self.c_timepoints, 0, 0) #Toggle add to lineage 0 create schnitzes 0

			#print("Entering while loop")
			#Enter Simulation Queue
			while list_index < len(self.old_cell_states):
				self.cs = self.old_cell_states[list_index]
				list_index += 1
				#print("list_index", list_index, "cs.get_time()", self.cs.get_time(), "of final_time = ", final_time)

				#If the cell is dead, do not simulate. Save if include_dead_cells toggled
				if self.cs.get_dead() >= 0 and include_dead_cells == 1:
					#print("cell dead")
					self.cell_states.append(self.cs)

				#If the cell has already simulated all its time, do not simulate and save
				elif self.cs.get_time() > final_time - dt:
					#print("cell hit final time", final_time, self.cs.get_time())
					self.cell_states.append(self.cs)

				#If the cell has divided, apply the appropriate division rule
				elif self.cs.get_divided() >= 0:
					#print("cell divided")
					daughter_cells = self.interface.partition(self.cs.get_divided(), self.cs)

					self.d1 = <LineageVolumeCellState>(daughter_cells[0])
					self.d2 = <LineageVolumeCellState>(daughter_cells[1])

					#Create a new timepoint array and simulate the first daughter and queue if it doesn't reach final time.
					self.c_truncated_timepoints = self.truncate_timepoints_less_than(self.c_timepoints, self.cs.get_time())
					self.simulate_daughter_cells(self.c_truncated_timepoints, 0, 0) #toggle add_to_lineage = 0 create_schnitzes = 0

				#Otherwise simulate some more
				else:
					#print("continuting simulation self.cs.get_time()", self.cs.get_time())
					self.c_truncated_timepoints = self.truncate_timepoints_less_than(self.c_timepoints, self.cs.get_time())
					self.r = self.SimulateSingleCell(self.cs, self.c_truncated_timepoints, 0)#Mode set to 0
					self.cs = self.r.get_final_cell_state()
					self.old_cell_states.append(self.cs)

			#Add samples to list to return
			samples.append(self.cell_states)
			#reset lists
			self.old_cell_states = self.cell_states
			self.cell_states = []

		return samples #Return list of samples

	#Python wrapper of the above
	def py_PropagateCells(self, np.ndarray timepoints, list initial_cell_states, LineageCSimInterface interface, np.ndarray sample_times, unsigned include_dead_cells):
		self.cell_states = []
		self.old_cell_states = []
		self.set_c_timepoints(timepoints)
		self.intialize_single_cell_interface(interface)
		self.initialize_single_cell_results_arrays(timepoints.shape[0])
		return self.PropagateCells(initial_cell_states, timepoints, sample_times, include_dead_cells)


	#Simulates a single cell trajectory, ignoring half the daughters every division.
	cdef list SingleCellLineage(self, LineageVolumeCellState initial_cell, double[:] timepoints):
		cdef double final_time = timepoints[timepoints.shape[0]-1]
		cdef unsigned list_index = 0
		cdef unsigned i

		self.cell_states = []
		self.old_cell_states = []

		#These will be used to store outputs
		cdef np.ndarray daughter_cells

		self.r = self.SimulateSingleCell(initial_cell, timepoints, 1) #Mode set to 1
		self.old_cell_states.append(self.r.get_final_cell_state())
		self.cell_states.append(self.r)

		while list_index < len(self.old_cell_states):
			self.cs = self.old_cell_states[list_index]
			list_index += 1

			#If the cell has already simulated all its time, do nothing
			if self.cs.get_time() >= final_time- 1E-9:
				pass
			#If the cell is dead, return the list of cells
			elif self.cs.get_dead() >= 0:
				return self.cell_states
			#If the cell has divided, apply the appropriate division rule
			elif self.cs.get_divided() >= 0:
				daughter_cells = self.interface.partition(self.cs.get_divided(), self.cs)

				#Create a new timepoint array and simulate a random daughter and queue if it doesn't reach final time.
				i = <unsigned>cyrandom.uniform_rv()>.5
				self.d1 = <LineageVolumeCellState>(daughter_cells[i])
				self.c_truncated_timepoints = self.truncate_timepoints_less_than(timepoints, self.cs.get_time())
				self.r = self.SimulateSingleCell(self.d1, self.c_truncated_timepoints, 1)
				self.cell_states.append(self.r)
				self.d1final = self.r.get_final_cell_state()

				# Add on the new daughter if final time wasn't reached.
				if self.d1final.get_time() < final_time + 1E-9:
					self.old_cell_states.append(self.d1final)
				else:
					warnings.warn("Daughter cell simulation went over the total time. Simulation has been discarded. Check for model errors.")

		return self.cell_states
	#Python wrapper of the above
	def py_SingleCellLineage(self, np.ndarray timepoints, LineageVolumeCellState initial_cell, LineageCSimInterface interface):
		self.set_c_timepoints(timepoints)
		self.intialize_single_cell_interface(interface)
		return self.SingleCellLineage(initial_cell, timepoints)


#Auxilary wrapper functions for quick access to Lineage Simulations

#SingleCellLineage simulates the trajectory of a single cell, randomly discarding one of its daughters every division.
def py_SingleCellLineage(timepoints, initial_cell_state = None, LineageModel Model = None, LineageCSimInterface interface = None, LineageSSASimulator simulator = None, return_dataframes = True):
	if Model == None and interface == None:
		raise ValueError('py_PropagateCells requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
	elif interface == None:
		interface = LineageCSimInterface(Model)
		interface.py_set_initial_time(timepoints[0])

	if initial_cell_state is None:
		initial_cell_state = LineageVolumeCellState(v0 = 1, t0 = 0, state = interface.py_get_initial_state())
	elif not isinstance(initial_cell_state, LineageVolumeCellState):
		raise ValueError("initial_cell_state must be of type LineageVolumeCellState or None (in which case it will default to the Model's initial state).")

	if simulator == None:
		simulator = LineageSSASimulator()
	result = simulator.py_SingleCellLineage(timepoints, initial_cell_state, interface)
	cell_lineage = []
	if return_dataframes:#Converts list of cell states into a Pandas dataframe
		try:
			import pandas
			df_list = [r.py_get_dataframe(Model = Model) for r in result]
			return pandas.concat(df_list)
		except ModuleNotFoundError:
			warnings.warn("return_dataframes=True requires that pandas be installed. Instead a numpy array is being returned (each column is a species, the last column is volume, and rows are cell states)")
	else:
		return result

#py_PropagateCells simulates an ensemble of growing dividing cells, returning only the cell states at the end of timepoints
#include_dead_cells toggles whether all dead cells accumulated along the way will also be returned.
#return data_frames returns all the results as a pandas dataframe. Otherwise results are returned as a list of LineageVolumeCellStates
def  py_PropagateCells(timepoints, initial_cell_states = [], LineageModel Model = None, LineageCSimInterface interface = None, LineageSSASimulator simulator = None, sample_times = 1, include_dead_cells = False, return_dataframes = True, return_sample_times = True):

	if Model == None and interface == None:
		raise ValueError('py_PropagateCells requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
	elif interface == None:
		interface = LineageCSimInterface(Model)
		interface.py_set_initial_time(timepoints[0])

	if isinstance(initial_cell_states, int):
		initial_cell_states = [LineageVolumeCellState(v0 = 1, t0 = 0, state = interface.py_get_initial_state())]*initial_cell_states
	elif (isinstance(initial_cell_states, list) and len(initial_cell_states) == 0):
		initial_cell_states = [LineageVolumeCellState(v0 = 1, t0 = 0, state = interface.py_get_initial_state())]
	elif not isinstance(initial_cell_states, list):
		raise ValueError("Initial Cell States must be a list of LineageVolumeCell states or and positive integer")
	if simulator == None:
		simulator = LineageSSASimulator()

	if isinstance(sample_times, int): #Return N=sample_times evenly spaced samples starting at the end of the simulation
		sample_times = timepoints[::-int(len(timepoints)/sample_times)]
		sample_times = np.flip(sample_times) #reverse the order
	else:
		sample_times = np.array(sample_times, dtype = np.double) #convert sample_times into doubles

	final_cell_state_samples = simulator.py_PropagateCells(timepoints, initial_cell_states, interface, sample_times, include_dead_cells)

	return_data = None
	if return_dataframes:#Converts list of cell states into a Pandas dataframe
		try:
			import pandas
			return_data = []
			for L in final_cell_state_samples:
				if len(L) > 0: darray = np.array([np.append(cs.py_get_state(), cs.py_get_volume()) for cs in L])
				if Model == None:
					warnings.warn("Without passing in a model, the data frame will not be indexable by species name.")
					df = pandas.DataFrame(darray)
				else:
					columns = Model.get_species_list()+["volume"]
					df = pandas.DataFrame(darray, columns = columns)
				return_data.append(df)
		except ModuleNotFoundError:
			warnings.warn("return_dataframes=True requires that pandas be installed. Instead a numpy array is being returned (each column is a species, the last column is volume, and rows are cell states)")
	else:
		return_data = final_cell_state_samples

	if return_sample_times:
		return return_data, sample_times
	else:
		return return_data

#SimulateCellLineage simulates a lineage of growing, dividing, and dieing cells over timepoints.
#The entire time trajectory of the simulation is returned as a Lineage which contains a binary tree of Schnitzes each containing a LineageSingleCellSSAResult.
def py_SimulateCellLineage(timepoints, initial_cell_states = [], initial_cell_count = 1, interface = None, Model = None):

	simulator = LineageSSASimulator()
	if Model == None and interface == None:
		raise ValueError('py_SimulateCellLineage requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
	elif interface == None:
		interface = LineageCSimInterface(Model)
		interface.py_set_initial_time(timepoints[0])

	if isinstance(initial_cell_states, int):
		initial_cell_states = [LineageVolumeCellState(v0 = 1, t0 = 0, state = interface.py_get_initial_state())]*initial_cell_states
	elif (isinstance(initial_cell_states, list) and len(initial_cell_states) == 0):
		initial_cell_states = [LineageVolumeCellState(v0 = 1, t0 = 0, state = interface.py_get_initial_state())]
	elif not isinstance(initial_cell_states, list):
		raise ValueError("Initial Cell States must be a list of LineageVolumeCell states or and positive integer")

	return simulator.py_SimulateCellLineage(timepoints, interface = interface, initial_cell_states = initial_cell_states)


#SimulateSingleCell performs an SSA simulation on a single cell until it divides, dies, or the final timepoint arrives.
def py_SimulateSingleCell(timepoints, Model = None, interface = None, initial_cell_state = None, return_dataframes = True):
	simulator = LineageSSASimulator()
	if Model == None and interface == None:
		raise ValueError('py_SimulateSingleCell requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
	elif interface == None:
		interface = LineageCSimInterface(Model)
		interface.py_set_initial_time(timepoints[0])

	if initial_cell_state == None:
		v = LineageVolumeCellState(v0 = 1, t0 = 0, state = interface.py_get_initial_state())

	result = simulator.py_SimulateSingleCell(timepoints, Model = Model, interface = interface, v = v)

	if return_dataframes:
		return result.py_get_dataframe(Model = Model)
	else:
		return result


#A simulator class for interacting cell lineages
cdef class InteractingLineageSSASimulator(LineageSSASimulator):

	#Used for Simulating Interacting lineages
	cdef int spec_ind, global_crn_initialized
	cdef unsigned num_global_species, num_interfaces, total_cell_count, num_global_crn_species
	cdef double[:] c_global_species, c_period_timepoints, active_lineages, sample_times
	cdef int[:, :] global_species_inds #stores global_species_inds[i, j] --> species index of interface j for global species i

	cdef double total_cell_volume, global_volume, leftover_global_volume, temp_volume, global_volume_param, average_dist_threshold, global_sync_period, dt

	cdef list interface_list, new_schnitzes, new_cell_states, lineage_list, samples, old_cell_state_list, new_cell_state_list, old_schnitz_list, new_schnitz_list
	cdef SingleCellSSAResult new_r, merge_r
	cdef LineageVolumeCellState new_cs
	cdef Schnitz new_s
	cdef double[:, :] c_global_species_array #output results for global species

	#variables for the simulators used for simulating CRNs in the global volume
	cdef VolumeSSASimulator global_ssa_simulator
	cdef DeterministicSimulator global_deterministic_simulator
	cdef ModelCSimInterface global_interface
	cdef int[:] global_species_global_crn_inds #stores global_species_global_crn_inds[i] --> species index of global_crn interface for global species i
	cdef double[:] global_crn_state
	cdef Volume global_volume_object
	cdef VolumeSSAResult global_crn_result, period_global_crn_result

	def __init__(self):
		self.global_crn_initialized = 0

	def get_global_species_array(self):
		#print("get_global_species_array")
		if self.c_global_species_array is None:
			warnings.warn("Global Species Array has not been created. Perhaps a simulation hasn't been run yet? Returning empty array.")
			return np.ndarray()
		else:
			return np.asarray(self.c_global_species_array)

	def get_global_crn_results(self):
		#print("get_global_crn_results")
		if self.global_crn_initialized == 0 or self.global_crn_result == None:
			warnings.warn("No Global simulation was performed. Will return None.")
			return None
		else:
			return self.global_crn_result

	def get_lineage_list(self):
		if self.lineage_list is not None:
			return self.lineage_list
		else:
			warnings.warn("There is no lineage_list - for Propogation results use get_samples(). Returning None.")
			return None

	def get_samples(self):
		if self.samples is not None:
			return self.samples, self.sample_times
		else:
			warnings.warn("There are no samples - for Lineage results use get_lineage_list(). Returning None.")
			return None

	#Functions to set up Global simulation
	def setup_global_volume_simulation(self, simulator, interface, global_species_global_crn_inds):
		cdef np.ndarray results
		cdef np.ndarray time = np.zeros(1)
		cdef np.ndarray volume = np.zeros(1)
		if isinstance(simulator, VolumeSSASimulator):
			self.global_ssa_simulator = simulator
			self.global_deterministic_simulator = None
			self.global_volume_object = Volume()
		elif isinstance(simulator, DeterministicSimulator):
			self.global_deterministic_simulator = simulator
			self.global_ssa_simulator = None
		else:
			raise ValueError("set_global_simulator requires a VolumeSSASimulator or DeterministicSimulator")

		if isinstance(interface, ModelCSimInterface):
			self.global_interface = interface
			self.global_crn_state = self.global_interface.get_initial_state()
			#print("global_crn_state", self.global_crn_state)

			if self.global_deterministic_simulator is not None:
				self.global_interface.py_prep_deterministic_simulation()
		else:
			raise ValueError("set_global_interface requires as ModelCSimInterface")

		global_species_global_crn_inds.astype(np.int, copy = False)
		self.global_species_global_crn_inds = global_species_global_crn_inds

		#instantiate VSR to store results
		self.num_global_crn_species = self.global_interface.get_number_of_species()
		results = np.zeros((1, self.num_global_crn_species))
		results[0, :] = self.global_crn_state[:]
		volume[0] = 1
		time[0] = 0

		self.global_crn_result = VolumeSSAResult(time, results, volume, 0)
		self.global_crn_initialized = 1



	#Initializes a new interface and sets all internal variables such as lineage and cell_lists appropriately
	cdef void switch_interface(self, unsigned interface_ind, create_schnitzes):
		#print("switch_interface", interface_ind)
		self.intialize_single_cell_interface(self.interface_list[interface_ind])

		self.old_cell_states = self.old_cell_state_list[interface_ind]
		self.new_cell_states = self.new_cell_state_list[interface_ind]

		if create_schnitzes == 1: #Note needed when propagating cells
			self.lineage = self.lineage_list[interface_ind]
			self.old_schnitzes = self.old_schnitz_list[interface_ind]
			self.new_schnitzes = self.new_schnitz_list[interface_ind]

			#print("self.lineage_list", self.lineage_list, "self.old_schnitzes", len(self.old_schnitzes))
		#print("self.old_cell_states", len(self.old_cell_states))

	#Calculates volume variables from the entire population of cells
	cdef void calculate_global_volumes(self):
		#print("calculate_global_volumes")
		cdef unsigned interface_ind = 0
		cdef unsigned i = 0
		cdef dead_cell_count = 0
		self.total_cell_count = 0
		self.total_cell_volume = 0 #reset total cell volume

		for interface_ind in range(self.num_interfaces):
			#print('interface ind', interface_ind, "self.old_cell_state_list[interface_ind]", self.old_cell_state_list[interface_ind])

			self.old_cell_states = self.old_cell_state_list[interface_ind]

			for i in range(len(self.old_cell_states)):
				self.cs = self.old_cell_states[i]
				self.total_cell_volume += self.cs.get_volume()
				self.total_cell_count += 1
				dead_cell_count += self.cs.get_dead()

			if dead_cell_count == len(self.old_cell_states):
				self.active_lineages[interface_ind] = 0 #toggle to say all cells in this lineage are dead

		if self.global_volume_param == 0: #If global volume is 0, assume global_volume = total_cell_volume for the entire simulation
			self.leftover_global_volume = 0
			self.global_volume = self.total_cell_volume

		elif self.total_cell_volume > self.global_volume_param:
			warnings.warn("Total cell volume exceeded global volume. All cells set to dead and simulation terminated.")
			#Set all cells to dead
			#Death by crowding
			for interface_ind in range(self.num_interfaces):
				self.old_cell_states = self.old_cell_state_list[interface_ind]
				for i in range(len(self.old_cell_states)):
					self.cs = self.old_cell_states[i]
					self.cs.set_dead(1)
				self.active_lineages[interface_ind] = 0 #toggle to say all cells in this lineage are dead
		else:
			self.global_volume = self.global_volume_param
			self.leftover_global_volume = self.global_volume - self.total_cell_volume


	#Calculates the total number of global species across cell states
	cdef void calculate_global_species_totals(self):
		#print("calculate_global_species_totals")
		cdef unsigned interface_ind = 0
		cdef unsigned list_index = 0
		cdef unsigned i = 0
		cdef int spec_ind = 0

		#Cycle through all the cell_states by interface
		for interface_ind in range(self.num_interfaces):
			self.old_cell_states = self.old_cell_state_list[interface_ind]
			#cycle through each cell
			for list_index in range(len(self.old_cell_states)):
				self.cs = self.old_cell_states[list_index]
				self.c_current_state = self.cs.get_state().copy()
				#Add the global species if they are in the model
				for i in range(self.num_global_species):
					spec_ind = self.global_species_inds[i, interface_ind]
					if spec_ind >= 0: #If the cell contains the global species, set its internal s_i to 0 and add that to the global total
						self.c_global_species[i] += self.c_current_state[spec_ind]
						self.cs.set_state_comp(0, spec_ind)

	#Synchronizes global species by redistributing them between different volumes, including the global volume
	cdef void synchronize_global_species(self):
		#print("synchronize_global_species")
		cdef unsigned i = 0
		cdef unsigned living_cells = 0
		cdef unsigned interface_ind

		#Calculate the global volumes
		self.calculate_global_volumes()
		#Calculate the global species totals
		self.calculate_global_species_totals()

		#check if all lineages are dead
		for interface_ind in range(self.num_interfaces):
			if self.active_lineages[interface_ind] == 1:
				living_cells = 1

		if living_cells:
			for i in range(self.num_global_species):
				#If the amount of a global species is above the threshold for stochastic distribute, distribute the average
				if self.total_cell_volume/self.global_volume*self.c_global_species[i]/self.total_cell_count > self.average_dist_threshold:
					self.c_global_species[i] = self.distribute_global_species_average(self.c_global_species[i], i)

				#Otherwise distribute stochastically
				else:
					self.c_global_species[i] = self.distribute_global_species_multinomial(self.c_global_species[i], i)

	#Distribute global species to their expected values
	#global_species is the number of species to distribue between all old_cell_states
	#will distribute these species and return the number of species left in the global volume
	cdef double distribute_global_species_average(self, double global_count, unsigned global_species_ind):
		#print("distribute_global_species_average")
		cdef unsigned i = 0
		cdef unsigned temp_count = 0
		cdef double new_global_species = global_count  #stores the number of global species passed to the global volume
		cdef unsigned list_index = 0
		cdef unsigned interface_ind = 0
		cdef int spec_ind = 0

		for interface_ind in range(self.num_interfaces):
			if self.active_lineages[interface_ind] == 1:
				self.old_cell_states = self.old_cell_state_list[interface_ind]
				for list_index in range(len(self.old_cell_states)):
					self.cs = self.old_cell_states[list_index]
					temp_count = int(self.cs.get_volume()/self.global_volume*global_count)

					spec_ind = self.global_species_inds[global_species_ind, interface_ind]
					#Add the global species if they are in the model and not dead
					if spec_ind >= 0 and self.cs.get_dead() < 0: #Check if the cell contains that species and do not distribute to dead cells
						new_global_species -= temp_count
						self.cs.set_state_comp(temp_count, spec_ind)

		return new_global_species

	#Distribute global species stochastically
	#global_count is the number of species to distribue between all old_cell_states
	#will distribute these species and return the number of species left in the global volume
	cdef double distribute_global_species_multinomial(self, double global_count, unsigned global_species_ind):
		#print("distribute_global_species_multinomial", global_count, global_species_ind)
		cdef unsigned i = 0
		cdef double rand
		cdef double temp_volume = 0
		cdef double new_global_species = 0  #stores the number of global species passed to the global volume
		cdef unsigned list_index = 0
		cdef unsigned interface_ind = 0
		cdef int spec_ind = 0

		while global_count >= 1: #Count may go less than 1 but greater than 0 if deterministic global crn simulation is used
			rand = cyrandom.uniform_rv()
			temp_volume = self.leftover_global_volume
			#randomly add to the global volume first because it is more probable
			if rand <= temp_volume/self.global_volume and self.global_volume_param > 0:
				#print("stuck here?",rand,  temp_volume, self.global_volume)
				new_global_species += 1
				global_count -= 1
				continue

			#cycle through cells
			for interface_ind in range(self.num_interfaces):
				#print("or here????",rand, temp_volume, self.global_volume)
				if self.active_lineages[interface_ind] == 1:
					self.old_cell_states = self.old_cell_state_list[interface_ind]
					for list_index in range(len(self.old_cell_states)):
						#print("or this place?")
						self.cs = self.old_cell_states[list_index]
						temp_volume += self.cs.get_volume()
						if rand <= temp_volume/self.global_volume:

							spec_ind = self.global_species_inds[global_species_ind, interface_ind]
							if spec_ind >= 0 and self.cs.get_dead() < 0: #Check if the cell contains that species and do not distribute to dead cells
								#print("is this here the place?", spec_ind)
								self.c_current_state = self.cs.get_state()
								self.cs.set_state_comp(self.c_current_state[spec_ind]+1, spec_ind) #add one to the cell state
								global_count -= 1 #decrement the global species
							else: #if the cell doesn't contain the species, add it to the new global species vector
								new_global_species += 1
								global_count -= 1
							break

		return new_global_species





	#Helper function to simulate the global CRN, if there is one
	#This CRN can be Stochastic or Deterministic
	cdef void simulate_global_volume_crn(self, np.ndarray timepoints):
		#print("Simulating Global CRN")
		cdef unsigned i = 0
		cdef unsigned spec_ind = 0
		cdef unsigned num_timepoints = timepoints.shape[0]
		cdef double final_time = timepoints[num_timepoints - 1]
		cdef unsigned num_global_crn_Species = self.global_interface.get_number_of_species()

		self.global_interface.set_initial_time(timepoints[0])
		self.global_interface.set_dt(timepoints[1]-timepoints[0])

		#Deterministic Simulation
		if self.global_deterministic_simulator is not None:
			#Set the state to the concentration=count/volume
			for i in range(self.num_global_species):
				spec_ind = self.global_species_global_crn_inds[i]
				self.global_crn_state[spec_ind] = self.c_global_species[i]/self.leftover_global_volume

			#Simulate
			self.global_interface.set_initial_state(np.asarray(self.global_crn_state))
			self.period_global_crn_result = self.global_deterministic_simulator.simulate(self.global_interface, timepoints)
			for i in range(self.num_global_crn_species):
				self.global_crn_state[i] = self.period_global_crn_result.get_result()[num_timepoints-1, i] #Set global crn state

			#reset self.global_species after simulation count = volume*concentration
			for i in range(self.num_global_species):
				spec_ind = self.global_species_global_crn_inds[i]
				self.c_global_species[i] = self.global_crn_state[spec_ind]*self.leftover_global_volume


		#Stochastic Simulation
		elif self.global_ssa_simulator is not None:
			#print("Stochastic Simulation")
			#Set the global species
			#print("Setting global species")
			for i in range(self.num_global_species):
				spec_ind = self.global_species_global_crn_inds[i]
				self.global_crn_state[spec_ind] = self.c_global_species[i]

			#Simulate
			#print("Simulate")
			self.global_interface.set_initial_state(np.asarray(self.global_crn_state))
			#print("here?", np.asarray(self.global_crn_state))
			self.global_volume_object.set_volume(self.leftover_global_volume) #set the global volume
			#print("here??", timepoints)
			#print("self.global_interface", self.global_interface, "self.global_volume_object", self.global_volume_object)
			self.period_global_crn_result = self.global_ssa_simulator.volume_simulate(self.global_interface, self.global_volume_object, timepoints)
			#print("updating crn_state")
			for i in range(self.num_global_crn_species):
				self.global_crn_state[i] = self.period_global_crn_result.get_result()[num_timepoints-1, i] #Set global crn state

			#print("updating global species")
			#reset self.global_species after simulation
			for i in range(self.num_global_species):
				spec_ind = self.global_species_global_crn_inds[i]
				self.c_global_species[i] = self.global_crn_state[spec_ind]


	#Helper function to simulate one sync-period of an interacting lineage
	cdef void simulate_interacting_lineage_period(self, double[:] timepoints, unsigned create_schnitzes):
		#print("SimulateInteractingLineagePeriod")
		cdef np.ndarray daughter_cells
		cdef unsigned interface_ind = 0
		cdef unsigned list_index = 0
		cdef unsigned num_timepoints = timepoints.shape[0]
		cdef double final_time = timepoints[timepoints.shape[0]-1]

		#Cycle through interfaces
		for interface_ind in range(self.num_interfaces):
			list_index = 0
			#print("start of simulate lineage period")
			self.switch_interface(interface_ind, create_schnitzes)#set the correct interface and internal variables
			#print("simulate interacting lineage period: len(old_cell_states) = ", len(self.old_cell_states))
			#If propogating cells and not saving trajectories, create a global memory buffer here
			if create_schnitzes == 0:
				self.initialize_single_cell_results_arrays(num_timepoints)

			#cycle through cells
			#print("entering while loop")
			while list_index < len(self.old_cell_states):
				self.cs = self.old_cell_states[list_index]
				if create_schnitzes:
					self.s = self.old_schnitzes[list_index]

				list_index += 1

				#If the cell is dead add it to the lineage
				if self.cs.get_dead() >= 0:
					#print("cell dead", "time", self.cs.get_time(), "volume", self.cs.get_volume())
					# Add daughters to the lineage
					if create_schnitzes: self.lineage.add_schnitz(self.s)

				#If the cell has divided right at teh end of the period, add it to the next period for division
				#Do not add to lineage because that will happen next period
				elif self.cs.get_divided() >= 0 and self.cs.get_time() > final_time - self.dt:
					#print("cell divided and time up", "time", self.cs.get_time(), "volume", self.cs.get_volume())
					self.new_cell_states.append(self.cs)
					if create_schnitzes:
						self.new_schnitzes.append(self.s)

				#If a cell has divided and still has time left in the period, simulate the daughters.
				elif self.cs.get_divided() >= 0:
					#print("cell divided - split 'er in two!", "time", self.cs.get_time(), "volume", self.cs.get_volume())
					if create_schnitzes:
						self.lineage.add_schnitz(self.s)
					#print("cd 1")
					daughter_cells = self.interface.partition(self.cs.get_divided(), self.cs)
					self.d1 = <LineageVolumeCellState>(daughter_cells[0])
					self.d2 = <LineageVolumeCellState>(daughter_cells[1])
					#print("cd 2")
					self.c_truncated_timepoints = self.truncate_timepoints_less_than(timepoints, self.cs.get_time())
					self.simulate_daughter_cells(self.c_truncated_timepoints, 0, create_schnitzes) #Toggle add to lineage False and create schnitzes appropriately

				#If the cell has reached its period time, add it to new_schnitzes and cell_states
				elif self.cs.get_time() > final_time - self.dt:
					#print("cell has reached period time", "time", self.cs.get_time(), "volume", self.cs.get_volume())
					self.new_cell_states.append(self.cs)
					if create_schnitzes:
						self.new_schnitzes.append(self.s)

				#If the cell isn't dead or divided or at period time simulate it more
				else:
					#If there is only one or two timepoints left, push to the next period.
					self.c_truncated_timepoints = self.truncate_timepoints_less_than(timepoints, self.cs.get_time())

					if len(self.c_truncated_timepoints) <= 2:
						self.new_cell_states.append(self.cs)
						if create_schnitzes:
							self.new_schnitzes.append(self.s)
					else:
						#print("Calling from Continue", "cs.get_time", self.cs.get_time(),"volume", self.cs.get_volume(), "final time", final_time, "mode", create_schnitzes)
						self.new_r = self.SimulateSingleCell(self.cs, self.c_truncated_timepoints, create_schnitzes) #create_schnitzes toggles the mode of the SimulateSingleCell to save/not save the entire trajectory
						self.new_cs = self.new_r.get_final_cell_state()
						self.new_cs.set_initial_vars(self.cs.get_initial_volume(), self.cs.get_initial_time())
						self.new_cs.set_time(self.cs.get_time())
						self.new_cs.set_volume(self.cs.get_volume())

						#After simulation, merge the two SSA results
						if create_schnitzes:
							#print("A")
							self.merge_r = SingleCellSSAResult(np.concatenate((self.s.get_time(), self.new_r.get_timepoints())),
								np.concatenate((self.s.get_data(), self.new_r.get_result())),
								np.concatenate((self.s.get_volume(), self.new_r.get_volume())),
								self.new_r.get_divided() >= 0)
							#print("after merged")
							self.merge_r.set_divided(self.new_r.get_divided())
							self.merge_r.set_dead(self.new_r.get_dead())
							self.new_s = self.merge_r.get_schnitz()
							self.new_s.set_parent(self.s.get_parent())
							#print("B")


						#Add schnitzes to the next period if they are done simulating
						if self.new_cs.get_time() > final_time - self.dt:
							#print("C")
							self.new_cell_states.append(self.new_cs)
							if create_schnitzes:
								self.new_schnitzes.append(self.new_s)

						#stay in the same period (perhaps they have died or divided)
						else:
							#print("D")
							self.old_cell_states.append(self.new_cs)
							if create_schnitzes:
								self.old_schnitzes.append(self.new_s)
							#print("D2")

			#print("end of period loop bfefor updating the lists")
			#After going through all the old_cell_states, make new_cell_states old.
			#Then reset the lists to be empty
			#print("end of simulate interacting lineage period: len(old_cell_states) = ", len(self.old_cell_states), "len new cell states", len(self.new_cell_states))
			#print("old_cell_state_list", self.old_cell_state_list)
			#print("new_cell_state_list", self.new_cell_state_list)

			self.old_cell_state_list[interface_ind] = list(self.new_cell_states)
			self.new_cell_state_list[interface_ind] = list()

			if create_schnitzes:
				#print("creating schnitzes?")
				self.old_schnitz_list[interface_ind] = list(self.new_schnitzes)
				self.new_schnitz_list[interface_ind] = list()
		#print("Period simulation complete.")



	cdef list SimulateInteractingCellLineage(self, list interface_list, 
											 list initial_cell_states, 
											 np.ndarray timepoints, 
											 double global_sync_period, 
											 np.ndarray global_species_inds, 
											 double global_volume_param, 
											 double average_dist_threshold):
		#print("Starting Interacting Lineage Simulation")

		cdef unsigned i = 0
		cdef unsigned j = 0
		cdef unsigned spec_ind = 0
		cdef unsigned list_index = 0
		cdef unsigned interface_ind = 0
		cdef unsigned current_time_index = 0

		cdef double final_time = timepoints[timepoints.shape[0]-1] #when the entire simulation ends (adding dt onto the end for rounding reasons)
		cdef double current_time = timepoints[0] #current time

		self.dt = timepoints[1] - timepoints[0]
		self.global_sync_period = global_sync_period #How often global species are synchronized
		cdef double period_time = timepoints[0]+self.global_sync_period #When the next sync period happens

		#Store seperate cell lists and lineages for each interface
		self.lineage_list = [] #stores one lineage for each interface
		self.old_cell_state_list = [] #stores one list of cell states for each interface
		self.new_cell_state_list = [] #stores one list of new cell states for each interface
		self.old_schnitz_list = [] #stores one list of schnitzes for each interface
		self.new_schnitz_list = []


		self.c_global_species_array = np.zeros((timepoints.shape[0], self.num_global_species)) #stores the global species at each timepoint
		self.global_species_inds = global_species_inds #stores global_species_inds[i, j] --> species index of interface j for global species i
		self.c_global_species = np.zeros(self.num_global_species) #stores the global species vector


		#These parameters are global because they are used by helper functions
		self.total_cell_volume = 0 #stores sum_i volume(cell_i)
		self.num_global_species = global_species_inds.shape[0]
		self.leftover_global_volume = 0 #stores global_volume - total_cell_volume
		self.global_volume = 0 #stores global_volume_param OR total_cell_volume if global_volume_param == 0.
		self.global_volume_param = global_volume_param
		self.average_dist_threshold = average_dist_threshold

		#Check that len sims == len initial_cells. As cells divide, they inherit their mother's CSimInterface
		#Now done in python wrappers
		#if len(self.interface_list) != len(initial_cell_states):
		#	raise ValueError(f"interface list (length {len(self.interface_list)}) [a list of LineageCSimInterfaces] must be the same length as initial cells (length {len(initial_cell_states)}) [[a list of LineageVolumeCellStates] for each interface].")
		#Check that global period is greater than dt
		if self.global_sync_period < self.dt:
			raise ValueError("global sync period must be larger than the timestep, dt, in the timepoints passed in for simulation.")
		#Check global sync period is smaller than the total time simulated
		if self.global_sync_period >= final_time:
			raise ValueError("global sync period must be smaller than the entire length of the timepoints passed in for simulation.")

		#Calculate global volume
		self.calculate_global_volumes()

		self.interface_list = interface_list
		self.num_interfaces = len(self.interface_list)
		for interface_ind in range(self.num_interfaces):
			self.lineage_list.append(Lineage())
			self.new_cell_state_list.append([])
			self.old_cell_state_list.append([])
			self.old_schnitz_list.append([])
			self.new_schnitz_list.append([])
		#print("Lineage List before simulation", self.lineage_list)
		self.active_lineages = np.ones(self.num_interfaces)

		#Timepoints for the first set of simulations
		self.c_timepoints = timepoints.copy()
		self.c_period_timepoints = self.truncate_timepoints_greater_than(self.c_timepoints, period_time)


		#Simulate the initial cell states
		#All cells simulated till either they die, divide or they hit period time
		#Cell states placed in self.old_cell_state_list[interface_ind] organized by interface_ind
		for interface_ind in range(self.num_interfaces):
			self.switch_interface(interface_ind, 1) #Very important to call this before the simulation. 1 to toggle to set up the schnitz lists and lineage

			#print("Before Simulating Initial Cells")
			#print("interface_ind", interface_ind, "len new_cell_states", len(self.new_cell_states), "len old cell states", len(self.old_cell_states))
			#print("interface_ind", interface_ind, "len new_schnitz_list", len(self.new_schnitzes), "len old_schnitz_list", len(self.old_schnitzes))

			self.simulate_cell_list(initial_cell_states[interface_ind], self.c_period_timepoints, 0, 1) #toggle add_to_lineage = 0 create_schnitzes = 1

			#print("After Simulating Initial Cells")
			#print("interface_ind", interface_ind, "len new_cell_states", len(self.new_cell_states), "len old cell states", len(self.old_cell_states))
			#print("interface_ind", interface_ind, "len new_schnitz_list", len(self.new_schnitzes), "len old_schnitz_list", len(self.old_schnitzes))

		#Main Simulation loop
		#print("Entering Main loop")
		while period_time <= final_time and current_time <= final_time:
			#print("beginning of the loop")
			#Update timepoints
			current_time = period_time

			#Do nothing and the loop will end
			if period_time > final_time:
				pass
			#If the total sim time isn't a multiple of the period time, update carefully
			elif period_time + self.global_sync_period >= final_time and period_time < final_time:
				period_time = final_time+self.dt
			#Otherwise update normally
			elif period_time + self.global_sync_period < final_time+self.dt:
				period_time = period_time + self.global_sync_period

			self.c_period_timepoints = self.truncate_timepoints_greater_than(self.c_timepoints, period_time)

			#Calculate global volume and synchronize global species across all cells in self.old_cell_state_list
			self.synchronize_global_species() #calculates global values and redistributes the species
			self.simulate_interacting_lineage_period(self.c_period_timepoints, 1) #create_schnitzes toggled to 1

			#If a simulator is set, do the global CRN simulation
			if self.global_ssa_simulator is not None or self.global_deterministic_simulator is not None:
				#print("global simulation")
				self.simulate_global_volume_crn(timepoints[(current_time<=timepoints)*(timepoints<=period_time)]) #results are stored in self.r
				#Save final results via concatenation with previous results

				#Only concatenate one thing at a time to debug old code below

				self.global_crn_result = VolumeSSAResult(
					np.concatenate((self.global_crn_result.get_timepoints(), self.period_global_crn_result.get_timepoints()[1:])),
					np.concatenate((self.global_crn_result.get_result(), self.period_global_crn_result.get_result()[1:, :])),
					np.concatenate((self.global_crn_result.get_volume(), self.period_global_crn_result.get_volume()[1:])), 0)
				#print("set.")

			#print("setting global species array")
			while current_time_index < timepoints.shape[0] and timepoints[current_time_index]<=period_time:
				for spec_ind in range(self.num_global_species):
					self.c_global_species_array[current_time_index, spec_ind] = self.c_global_species[spec_ind]
				current_time_index += 1
			#print("complete.")


		#print("loop complete")
		#Add the final schnitzes to their lineages
		for interface_ind in range(self.num_interfaces):
			self.switch_interface(interface_ind, 1)

			for list_index in range(len(self.old_schnitzes)):
				self.s = self.old_schnitzes[list_index]
				self.lineage.add_schnitz(self.s)

		#print("about to return lineage list", len(self.lineage_list), (<Lineage>self.lineage_list[0]).py_size)
		return self.lineage_list

	#Python accessor

	def py_SimulateInteractingCellLineage(self, np.ndarray timepoints, 
										  list interface_list, 
										  list initial_cell_states, 
										  double global_sync_period, 
										  np.ndarray global_species_inds, 
										  double global_volume_param, 
										  double average_dist_threshold):
		self.set_c_timepoints(timepoints)
		lineage_list = self.SimulateInteractingCellLineage(interface_list, 
														   initial_cell_states, 
														   timepoints, 
														   global_sync_period, 
														   global_species_inds, 
														   global_volume_param, 
														   average_dist_threshold)

		#print("lineage_list in py_Simulate...", len(lineage_list))
		return lineage_list


	#Simulates an ensemble of interacting cells over some amount of time.
	#Returns the existing cell states at all times in the array sample_times
	#  dead cells are included based upon the include_dead_cells parameter (1 = included, otherwise = excluded)
	cdef list PropagateInteractingCells(self, list interface_list,
										list initial_cell_states,
										np.ndarray timepoints,
										double[:] sample_times,
										double global_sync_period,
										np.ndarray global_species_inds,
										double global_volume_param,
										double average_dist_threshold):
		#print("Starting PropagateInteractingCells")
		#Final result will be returned here as a list of lists
		self.samples = []
		self.sample_times = sample_times

		cdef unsigned sample_ind = 0
		cdef unsigned i = 0
		cdef unsigned j = 0
		cdef unsigned spec_ind = 0
		cdef unsigned list_index = 0
		cdef unsigned interface_ind = 0

		cdef unsigned num_timepoints = timepoints.shape[0]
		cdef double final_time = timepoints[timepoints.shape[0]-1] #when the entire simulation ends (adding dt onto the end for rounding reasons)
		cdef double current_time = timepoints[0] #current time
		cdef double next_sample_time = self.sample_times[0] #when the next sample is taken

		cdef unsigned num_samples = self.sample_times.shape[0]
		cdef unsigned first_simulation = 0

		self.dt = timepoints[1] - timepoints[0]
		self.global_sync_period = global_sync_period #How often global species are synchronized
		cdef double period_time = timepoints[0]+self.global_sync_period#When the next sync period happens

		#Store seperate cell lists and lineages for each interface
		self.old_cell_state_list = [] #stores one list of cell states for each interface
		self.new_cell_state_list = [] #stores one list of new cell states for each interface
		self.interface_list = interface_list
		self.global_species_inds = global_species_inds #stores global_species_inds[i, j] --> species index of interface j for global species i
		self.active_lineages = np.ones(self.num_interfaces) #Just in case everything dies!

		#These parameters are global because they are used by helper functions
		self.total_cell_volume = 0 #stores sum_i volume(cell_i)
		self.num_global_species = global_species_inds.shape[0]
		self.c_global_species = np.zeros(self.num_global_species) #stores the global species vector
		self.c_global_species_array = np.zeros((num_samples, self.num_global_species))#stores the global species at each sample time
		self.leftover_global_volume = 0 #stores global_volume - total_cell_volume
		self.global_volume = 0 #stores global_volume_param OR total_cell_volume if global_volume_param == 0.
		self.global_volume_param = global_volume_param
		self.average_dist_threshold = average_dist_threshold

		#Check that len sims == len initial_cells. As cells divide, they inherit their mother's CSimInterface
		#Now done in python wrappers
		#if len(self.interface_list) != len(initial_cell_states):
		#	raise ValueError(f"interface list (length {len(self.interface_list)}) [a list of LineageCSimInterfaces] must be the same length as initial cells (length {len(initial_cell_states)}) [[a list of LineageVolumeCellStates] for each interface].")
		#Check that global period is greater than dt
		if self.global_sync_period < self.dt:
			raise ValueError("global sync period must be larger than the timestep, dt, in the timepoints passed in for simulation.")
		#Check global sync period is smaller than the total time simulated
		if self.global_sync_period >= final_time:
			raise ValueError("global sync period must be smaller than the entire length of the timepoints passed in for simulation.")
		#Check that sample_times and timepoints are compatable
		if len(timepoints) < self.sample_times.shape[0] or timepoints[timepoints.shape[0]-1] < self.sample_times[self.sample_times.shape[0]-1] or timepoints[0] > self.sample_times[0]:
			raise ValueError("sample_times must be a subset of the timepoints array.")

		#Calculate global volume
		self.calculate_global_volumes()

		self.interface_list = interface_list
		self.num_interfaces = len(self.interface_list)
		for interface_ind in range(self.num_interfaces):
			self.new_cell_state_list.append([])
			self.old_cell_state_list.append([])

		#If the first sample is just the start of the simulation, append the initial condition
		if self.sample_times[0] <= timepoints[0]+self.dt:
			self.samples.append(list(initial_cell_states))
			self.c_global_species_array[0, :] = self.c_global_species
			sample_ind += 1
			next_sample_time = self.sample_times[sample_ind]

		if next_sample_time < period_time:
			final_time = next_sample_time
		elif period_time < next_sample_time:
			final_time = period_time

		#Timepoints for the first set of simulations
		#print("start of popagating: final time", final_time, "period_time", period_time, "next_sample_time", next_sample_time)

		self.c_period_timepoints = self.truncate_timepoints_greater_than(timepoints, final_time)

		#Simulate the initial cell states
		#All cells simulated till either they die, divide or they hit period time
		#Cell states placed in self.old_cell_state_list[interface_ind] organized by interface_ind
		for interface_ind in range(self.num_interfaces):
			self.switch_interface(interface_ind, 0) #Very important to call this before the simulation. 0 toggles no schnitzes or lineages
			self.initialize_single_cell_results_arrays(num_timepoints)
			self.simulate_cell_list(initial_cell_states[interface_ind], self.c_period_timepoints, 0, 0)



		#Main Simulation loop
		while sample_ind < num_samples:
			list_index = 0

			#Synchronize global species and calculate volumes
			if final_time == period_time:
				period_time += global_sync_period
				#print("synchronizing global species")
				self.synchronize_global_species()

			#Add samples to list to return if at sample_time
			if next_sample_time == final_time:
				sample_ind += 1
				next_sample_time = self.sample_times[sample_ind]

				self.samples.append(list(self.old_cell_state_list))
				self.c_global_species_array[sample_ind, :] = self.c_global_species

				#If global CRN simulation, update the results
				if self.global_crn_initialized == 1:
					#print("updating global crn result")
					#print("self.global_crn_result.get_timepoints()", self.global_crn_result.get_timepoints())
					#print("np.ndarray(final_time)", np.array(final_time))
					self.global_crn_result = VolumeSSAResult(
						np.append(self.global_crn_result.get_timepoints(), np.array(final_time)),
						np.concatenate((self.global_crn_result.get_result(), np.array(self.global_crn_state)[None, ...])),
						np.append(self.global_crn_result.get_volume(), np.array(self.leftover_global_volume)), 0)


			#If a simulator is set, do the global CRN simulation
			if self.global_crn_initialized == 1:
				self.simulate_global_volume_crn(timepoints[(current_time<=timepoints)*(timepoints<=final_time)]) #results are stored in self.r
				#Save final results via concatenation with previous results
				#print("self.global_crn_result.get_result().shape", self.global_crn_result.get_result().shape[0], self.global_crn_result.get_result().shape[1])
				#print("self.r.get_result().shape", self.period_global_crn_result.get_result().shape[0], self.period_global_crn_result.get_result().shape[1])
				#self.global_crn_result = VolumeSSAResult(
				#	np.concatenate((self.global_crn_result.get_timepoints(), self.period_global_crn_result.get_timepoints())),
				#	np.concatenate((self.global_crn_result.get_result(), self.period_global_crn_result.get_result())),
				#	np.concatenate((self.global_crn_result.get_volume(),self.period_global_crn_result.get_volume())), 0)
			#Choose timepoints for next simulation

			#Go to the first of next_sample_time and period_time and update appropriately
			current_time = final_time
			if next_sample_time <= period_time:
				final_time = next_sample_time
			if period_time <= next_sample_time:
				final_time = period_time

			#print("times updated in loop:", sample_ind, "final_time", final_time, "period_time", period_time, "next_sample_time", next_sample_time)

			self.c_timepoints = self.truncate_timepoints_greater_than(timepoints, final_time+1E-10)

			#Enter Simulation Queue
			self.simulate_interacting_lineage_period(self.c_timepoints, 0) #create_schnitzes toggled to off

#

		#print("about to return samples", "len samples, len samples[0], len samples[1]", len(self.samples), len(self.samples[0]), len(self.samples[1]))
		return self.samples
	#Python Accessor
	def py_PropagateInteractingCells(self, np.ndarray timepoints,
									 list interface_list,
									 list initial_cell_states, sample_times,
									 double global_sync_period,
									 np.ndarray global_species_inds,
									 double global_volume_param,
									 double average_dist_threshold):
		self.set_c_timepoints(timepoints)
		cell_sample_list =  self.PropagateInteractingCells(interface_list,
										initial_cell_states,
										timepoints,sample_times,
										global_sync_period, global_species_inds,
										global_volume_param,
										average_dist_threshold)
		return cell_sample_list



def py_set_up_InteractingLineage(global_species = [], interface_list = [],
								model_list = [], initial_cell_states = [],
								t0 = 0, simulator = None,
								global_species_inds = None,
								global_volume_simulator = "stochastic",
								global_volume_model = None):
	#Set up main simulator if needed
	if simulator == None:
		simulator = InteractingLineageSSASimulator()

	#set up global simulator
	if global_volume_model is not None:

		if len(global_species) == 0:
			warnings.warn(("Attemping InteractingCellL Simulation without ",
						   "global species defined. Consider a different ",
						   "simulator for increased efficiency."))

		global_species_global_crn_inds = np.zeros(len(global_species),
												  dtype = np.int32)

		for i in range(len(global_species)):
			s = global_species[i]
			ind = global_volume_model.get_species_index(s)
			#print("s", ind)
			if ind in [None, -1]:
				warnings.warn(f"Global Species {s} not in global_volume_model. Species is being added to Model with initial condition 0.")
				global_volume_model._add_species(s)
				global_volume_model.set_species({s: 0})
				global_volume_model.py_initialize()
				ind = global_volume_model.get_species_index(s)
			global_species_global_crn_inds[int(i)] = int(ind)

		global_interface = ModelCSimInterface(global_volume_model)

		if global_volume_simulator == "stochastic":
			global_volume_simulator = VolumeSSASimulator()
		elif global_volume_simulator == "deterministic":
			global_volume_simulator = DeterministicSimulator

		simulator.setup_global_volume_simulation(global_volume_simulator,
							global_interface, global_species_global_crn_inds)

	#Set up interface list from models and ensure compatability with global species
	if len(model_list) == 0 and len(interface_list) == 0:
		raise ValueError("Missing Required Keyword Arguments:models = [LineageModel] or interface_list = [LineageCSimInterface]")
	elif len(interface_list) == 0:
		interface_list = [LineageCSimInterface(m) for m in model_list]
		for m in model_list: #Initialize models
			m.py_initialize()
		global_species_inds = np.zeros((len(global_species), len(model_list)),
									    dtype = np.int32)
		for i in range(len(global_species)):
			for j in range(len(model_list)):
				m = model_list[j]
				s = global_species[i]
				ind = m.get_species_index(s)
				if ind != None:
					global_species_inds[i, j] = ind
				else:
					global_species_inds[i, j] = -1
	elif len(model_list) != 0 and len(interface_list) != 0:
		raise ValueError("Must call py_SimulateInteractingCellLineage with either the keyword argument model_list or the keyword argument interface_list, not both.")
	elif len(model_list) == 0 and (global_species_inds == None or global_species_inds.shape[1] != len(interface_list)):
		raise ValueError("When calling py_SimulateInteractingCellLineage with the keyword argument interface_list, the argument global_species_inds is required where global_species_inds[i, j] corresponds the species index of the ith global species in the jth interface.")
	elif len(global_species) == 0 and global_species_inds == None:
		warnings.warn('Calling SimulateInteractintCellLineage without any global species defined. Use the global_species or global_species_inds keywords.')
		global_species_inds = np.array(dtype = np.int32)

	#Set up initial cell states
	if len(initial_cell_states) == len(interface_list) and  isinstance(initial_cell_states[0], int):
		initial_cell_counts = initial_cell_states
		initial_cell_states = []

		for i in range(len(interface_list)):
			initial_cell_states.append([])
			interface = interface_list[i]
			initial_state = interface.py_get_initial_state()
			for j in range(initial_cell_counts[i]):
				lvcs = LineageVolumeCellState(v0 = 1.0, t0 = t0, state = initial_state.copy())
				initial_cell_states[i].append(lvcs)
	elif len(initial_cell_states) == 0 and len(interface_list) > 0:
		warnings.warn("Calling py_SimulateInteractintCellLineage without any initial_cell_states. Defaulting to creating one initial cell for each interface.")
		initial_cell_states = [LineageVolumeCellState(v0 = 1, t0 = 0, state = i.py_get_initial_state()) for i in interface_list]
	elif len(initial_cell_states) > len(interface_list):
		raise ValueError("When passing in more initial_cell_states than models, the keyword argument interface_inds is also required where interface_inds[i] corresponds to the index of the interface/model beloning to initial_cell_state[i].")


	return interface_list, simulator, initial_cell_states, global_species_inds

#Auxilary Python Function
def py_PropagateInteractingCells(timepoints, global_sync_period,
								 sample_times = 1, simulator = None,
								 global_volume_simulator = "stochastic",
								 global_volume_model = None, global_volume = 0,
								 global_species = None, global_species_inds = None,
								 average_dist_threshold = 1.0,
								 interface_list = None, model_list = None,
								 initial_cell_states = None,
								 return_dataframes = False,
								 return_sample_times = True):
	if global_species is None:
		global_species = []
	if interface_list is None:
		interface_list = []
	if model_list is None:
		model_list = []
	if initial_cell_states is None:
		initial_cell_states = []

	interface_list, simulator, initial_cell_states, global_species_inds = \
		py_set_up_InteractingLineage(global_species = global_species,
									 interface_list = interface_list,
									 model_list = model_list,
									 initial_cell_states = initial_cell_states,
									 simulator = simulator,
									 global_species_inds = global_species_inds,
									 global_volume_simulator = global_volume_simulator,
									 global_volume_model = global_volume_model,
									 t0 = timepoints[0])

	if isinstance(sample_times, int): #Return N=sample_times evenly spaced samples starting at the end of the simulation
		if sample_times == 1:
			sample_times = np.array(timepoints[-1])
		else:
			sample_times = timepoints[::-int(len(timepoints)/(sample_times-1))]
			sample_times = np.flip(sample_times) #reverse the order
	else:
		sample_times = np.array(sample_times, dtype = np.double) #convert sample_times into doubles

	final_cell_state_samples = simulator.py_PropagateInteractingCells(timepoints,
											interface_list, initial_cell_states,
											sample_times, global_sync_period,
											global_species_inds,
											global_volume,
											average_dist_threshold)

	if global_volume_model is not None:
		global_results = simulator.get_global_crn_results()
		#print("global_results returned")
	else:
		global_species_array = simulator.get_global_species_array()
		if return_dataframes:
			global_results = pandas.DataFrame(data = global_species_array, columns = global_species)
		else:
			global_results = global_species_array



	if return_sample_times:
		return final_cell_state_samples, sample_times, global_results, simulator
	else:
		return final_cell_state_samples, global_results, simulator


#Auxilary Python Function
def py_SimulateInteractingCellLineage(timepoints, global_sync_period,
	simulator = None, global_volume_simulator = "stochastic", global_volume_model = None,
	global_volume = 0, global_species = [], global_species_inds = None, average_dist_threshold = 1.0,
	interface_list = [], model_list = [], initial_cell_states = [], return_dataframes = False):

	interface_list, simulator, initial_cell_states, global_species_inds = py_set_up_InteractingLineage(global_species = global_species, interface_list = interface_list, model_list = model_list, initial_cell_states = initial_cell_states,
		simulator = simulator, global_species_inds = global_species_inds, global_volume_simulator = global_volume_simulator, global_volume_model = global_volume_model, t0 = timepoints[0])

	lineage_list = simulator.py_SimulateInteractingCellLineage(timepoints, interface_list, initial_cell_states, global_sync_period, global_species_inds, global_volume, average_dist_threshold)
	#print("auxilary func lineage_list returned")
	if global_volume_model is not None:
		global_results = simulator.get_global_crn_results()
		#print("global_results_returned")
	else:
		global_species_array = simulator.get_global_species_array()
		if return_dataframes:
			global_results = pandas.DataFrame(data = global_species_array, columns = global_species)
		else:
			global_results = global_species_array

	return lineage_list, global_results, simulator