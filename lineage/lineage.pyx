# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
cimport numpy as np
from bioscrape.simulator cimport VolumeCellState, VolumeSSAResult, DelayVolumeSSAResult, VolumeSplitter, CSimInterface, ModelCSimInterface, DelayVolumeSSASimulator, VolumeSSASimulator, DelayQueue, DelayVolumeCellState
from bioscrape.simulator import VolumeCellState, VolumeSSAResult, DelayVolumeSSAResult, VolumeSplitter, CSimInterface, ModelCSimInterface, DelayVolumeSSASimulator, VolumeSSASimulator, DelayQueue, DelayVolumeCellState

from bioscrape.types cimport Model, Volume, Schnitz, Lineage, Propensity, Term, Rule
#from bioscrape.types import Model, Volume, Schnitz, Lineage, Propensity, Term, Rule

from bioscrape.types import sympy_species_and_parameters, parse_expression
from bioscrape.random cimport normal_rv

cimport bioscrape.random as cyrandom
import bioscrape.random as cyrandom

from bioscrape.vector cimport vector 
#from vector cimport vector 

import warnings


#Events are general objects that happen with some internal propensity but are not chemical reactions
cdef class Event():

	cdef initialize(self, dict event_params, dict species_indices, dict parameter_indices):
		raise NotImplementedError("VolumeReactions Must be subclassed")

	def get_species_and_parameters(self, dict event_fields):
		raise NotImplementedError("VolumeReactions Must be subclassed")

	#cdef Propensity get_propensity(self):
	#	return <Propensity>self.propensity

	#Meant to be subclassed if an event is supposed to do something.
	cdef double evaluate_event(self, double *state, double *params, double volume, double time):
		return 0

#Volume Events are stochastic events which alter a single cell's volume.
cdef class VolumeEvent(Event):
	cdef double evaluate_event(self, double *state, double *params, double volume, double time):
		return self.get_volume(state, params,volume, time)

	cdef double get_volume(self, double *state, double *params, double volume, double time):
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

	cdef double get_volume(self, double *state, double *params, double volume, double time):
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

	cdef double get_volume(self, double *state, double *params, double volume, double time):
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

	cdef double get_volume(self, double *state, double *params, double volume, double time):
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
	cdef double evaluate_event(self, double *state, double *params, double volume, double time):
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
	cdef double evaluate_event(self, double *state, double *params, double volume, double time):
		return 1

#Dummy class to help with inheritance, compilation, and code simplification. Does nothing
cdef class LineageRule():
	def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
		raise NotImplementedError("LineageRule must be subclassed")

	def get_species_and_parameters(self, dict fields):
		raise NotImplementedError("LineageRule must be subclassed")

#Volume Rules occur every dt (determined by the simulation timepoints) and update the volume
cdef class VolumeRule(LineageRule):
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt):
		raise NotImplementedError("get_volume must be implemented in VolumeRule Subclasses")

cdef class LinearVolumeRule(VolumeRule):
	cdef unsigned has_noise
	cdef unsigned growth_rate_ind
	cdef unsigned noise_ind
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt):
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
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt):
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
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt):
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
	cdef double get_volume(self, double *state, double *params, double volume, double time, double dt):
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

	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume):
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
	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume):
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
	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume):
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
	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume):
		if self.has_noise > 0:
			if volume - initial_volume >= params[self.threshold_ind]+normal_rv(0, params[self.threshold_noise_ind]):
				return 1
			else:
				return 0
		else:
			if volume - initial_volume >= params[self.threshold_ind] - 1E-9:
				#print("time:", time, "volume", volume, "initial_time", initial_time, "initial_volume", initial_volume)
				#print("Cell Divided: DeltaV=", volume - initial_volume, "threshold=", params[self.threshold_ind])
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
	cdef unsigned check_divide(self, double *state, double *params, double time, double volume, double initial_time, double intial_volume):
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
	cdef unsigned check_dead(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume):
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
	cdef unsigned check_dead(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume):
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

	cdef unsigned check_dead(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume):
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
	cdef unsigned check_dead(self, double *state, double *params, double time, double volume, double initial_time, double initial_volume):
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

	def __init__(self, filename = None, species = [], reactions = [], parameters = [], rules = [], events = [], global_species = [], global_volume = None, sbml_filename = None, initial_condition_dict = None, input_printout = False, initialize_model = True):


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

		
		if len(global_species) > 0 and global_volume == None:
			warnings.warn("global species added to LineageModel without the global_volume keyword being passed in. Global volume is defaulting to 1 and no reaction rates will be rescaled.")
		elif len(global_species) == 0 and global_volume != None:
			warnings.warn("Setting global_volume without passing in the global_species keyword to LineageModel will do nothing unless you manually add global rections with the LineageModel.create_global_reaction function.")
		
		#Seperate reactions with global species inputs
		local_reactions = []
		global_reactions = []
		self.global_species = global_species
		if len(global_species) > 0:
			global_rxn_count = 0
			for rxn in reactions:
				reactants = rxn[0]
				if len(rxn) > 4:
					delay_reactants = rxn[5]
				else:
					delay_reactants = []
				if len([r for r in global_species if r in reactants]) > 0:
					global_reactions.append(rxn)
				else:
					local_reactions.append(rxn)
		else:
			local_reactions = reactions

		#Call super constructor
		super().__init__(filename = filename, species = species, reactions = local_reactions, parameters = parameters, rules = original_rules, initial_condition_dict = initial_condition_dict, sbml_filename = sbml_filename,  input_printout = input_printout, initialize_model = False)

		if global_volume == None:
			self.global_volume = 1.0
		else:
			self.global_volume = global_volume
		self._add_param("global_volume")
		self.set_parameter("global_volume", self.global_volume)

		#add global reactions (which are recaled by 1/volume_global^[# global species in reactants])
		global_rxn_count = 0
		for rxn in global_reactions:
			self.create_global_reaction(rxn, volume_param = "global_volume", volume_value = self.global_volume, identifier = global_rxn_count, global_species = self.global_species)
			global_rxn_count += 1

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

	def add_event(self, Event event_object, dict event_param_dict, Propensity prop_object, dict propensity_param_dict, str event_type = None, VolumeSplitter volume_splitter = None):
		self.initialized = False

		species_names_e, param_names_e = event_object.get_species_and_parameters(event_param_dict)
		species_names_p, param_names_p = prop_object.get_species_and_parameters(propensity_param_dict)

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
		if event_type in ["", "death", "DeathEvent", "death event", "Death Event", "default", "Default"]:
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
		if event_type in ["", "division", "Division", "DivisionEvent", "division event", "Division Event", "deafult"]:
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

		if event_type in ["linear", "Linear", "Linear Volume", "linear volume", "LinearVolume" "LinearVolumeEvent", "linear volume event", "Linear Volume Event"]:
			self._param_dict_check(event_params, "growth_rate", "DummyVar_LinearVolumeEvent")
			event_object = LinearVolumeEvent()
		elif event_type in ["multiplicative", "multiplicative", "Multiplicative Volume", "multiplicative volume", "MultiplicativeVolume", "MultiplicativeVolumeEvent", "Multiplicative Volume Event", "multiplicative volume event"]:
			self._param_dict_check(event_params, "growth_rate", "DummyVar_MultiplicativeVolumeEvent")
			event_object = MultiplicativeVolumeEvent()
		elif event_type in ["general", "General", "General Volume", "general volume", "GeneralVolume", "GeneralVolumeEvent", "General Volume Event", "general volume event"]:
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
		if "division" in rule_type or "Division" in rule_type:
			if volume_splitter == None:
				raise ValueError("DivisionRules must be added with a volume splitter object in add_lineage_rule.")
			rule_object.initialize(rule_param_dict, self.species2index,  self.params2index)
			self.division_rules_list.append((rule_object, volume_splitter))
		else:
			rule_object.initialize(rule_param_dict, self.species2index,  self.params2index)
			if "death" in rule_type or "Death" in rule_type:
				self.death_rules.append(rule_object)
			elif "volume" in rule_type or "Volume" in rule_type:
				self.volume_rules.append(rule_object)
			else:
				raise ValueError("add_lineage_rule only takes rules of type 'DeathRule', 'DivisionRule', and 'VolumeRule'. For Other rule types, consider trying Model.add_rule.")

	def create_death_rule(self, str rule_type, dict rule_param_dict):
		if rule_type in ["species", "Species", "SpeciesDeathRule"]:
			self._param_dict_check(rule_param_dict, "specie", "DummyVar_SpeciesDeathRule")
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_SpeciesDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_SpeciesDeathRule")
			rule_object = SpeciesDeathRule()
		elif rule_type in ["param", "parameter", "Param", "Parameter", "ParamDeathRule"]:
			self._param_dict_check(rule_param_dict, "param", "DummyVar_ParamDeathRule")
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_ParamDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_ParamDeathRule")
			rule_object = ParamDeathRule()
		elif rule_type in ["general", "General", "GeneralDeathRule"]:
			rule_object = GeneralDeathRule()
		else:
			raise ValueError("Unknown DeathRule type: "+str(rule_type))

		self.add_lineage_rule(rule_object, rule_param_dict, rule_type = "death")

	def create_division_rule(self, str rule_type, dict rule_param_dict, VolumeSplitter volume_splitter):
		if rule_type in ["time", "Time", "TimeDivisionRule"]:
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_TimeDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_TimeDeathRule")
			rule_object = TimeDivisionRule()
		elif rule_type in ["volume", "Volume", "VolumeDivisionRule"]:
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_VolumeDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_VolumeDeathRule")
			rule_object = VolumeDivisionRule()
		elif rule_type in ["delta", "Delta", "deltaV", "DeltaV", "DeltaVDivisionRule"]:
			self._param_dict_check(rule_param_dict, "threshold", "DummyVar_DeltaVDeathRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_DeltaVDeathRule")
			rule_object = DeltaVDivisionRule()
		elif rule_type in ["general", "General", "GeneralDivisionRule"]:
			rule_object = GeneralDivisionRule()
		else:
			raise ValueError("Unknown DivisionRule type: "+str(rule_type))
		self.add_lineage_rule(rule_object, rule_param_dict, rule_type = 'division', volume_splitter = volume_splitter)

	def create_volume_rule(self, str rule_type, dict rule_param_dict):
		if rule_type in ["linear", "Linear", "LinearVolumeRule"]:
			self._param_dict_check(rule_param_dict, "growth_rate", "DummyVar_LinearVolumeRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_LinearVolumeRule")
			rule_object = LinearVolumeRule()
		elif rule_type in ["multiplicative", "MultiplicativeVolume", "MultiplicativeVolumeRule"]:
			self._param_dict_check(rule_param_dict, "growth_rate", "DummyVar_MultiplicativeVolumeRule")
			if "noise" in rule_param_dict:
				self._param_dict_check(rule_param_dict, "noise", "DummyVar_MultiplicativeVolumeRule")
			rule_object = MultiplicativeVolumeRule()
		elif rule_type in ["assignment", "Assignment", "AssignmentVolumeRule"]:
			rule_object = AssignmentVolumeRule()
		elif rule_type in ["ode", "ODE", "ODEVolumeRule"]:
			rule_object = ODEVolumeRule()
		else:
			raise ValueError("Unknown VolumeRule type: "+str(rule_type))
		self.add_lineage_rule(rule_object, rule_param_dict, rule_type = 'volume')

	
	def create_global_reaction(self, rxn, volume_param = "global_volume", volume_value = 1, identifier = "", global_species = None):

		if len(rxn) == 4:
			reactants, products, propensity_type, propensity_param_dict = rxn
			delay_type, delay_reactants, delay_products, delay_param_dict = None, None,  None, None
		elif len(rxn) == 8:
			reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict = rxn
		else:
			raise ValueError("Reaction Tuple of the wrong length! Must be of length 4 (no delay) or 8 (with delays). See BioSCRAPE Model API for details.")
		
		if global_species == None:
			global_species = self.global_species

		if len(global_species) == 0:
			warnings.warn("No global species defined for this model or passed into create_global_reaction. Defaulting to non-global reaction.")
			self.create_reaction(reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict)
		elif "k" not in propensity_param_dict:
			warnings.warn("create_global_reaction only works with propensities that have a rate parameter 'k' in their param_dictionary. propensity_type="+propensity_type+" either doesn't have the proper parameter or is incompatible with automatic global reactions. This reaction will be added but not rescaled.") 
			self.create_reaction(reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict)
		else:
			old_val = propensity_param_dict["k"]
			try:
				float(old_val)
				float_val = True
			except ValueError:
				float_val = False

			rate_var = "global_reaction_"+propensity_type+"_k_rescaled_"+str(identifier)
			n_global = len([r for r in reactants if r in global_species])
			
			self._add_param(rate_var)
			if float_val:
				rule_equation = "_"+rate_var + "=" + str(old_val)+"/"+"(_"+volume_param+"^"+str(n_global)+")"
				self.set_parameter(rate_var, 1.*old_val/(volume_value**n_global))
			else:
				rule_equation = "_"+rate_var + "= _"+old_val+"/"+"(_"+volume_param+"^"+str(n_global)+")"
				self._add_param(old_val)
				self.set_parameter(rate_var, 1./(volume_value**n_global))
			propensity_param_dict["k"] = rate_var
			

			self._add_param(volume_param)
			self.set_parameter(volume_param, volume_value)
			self.create_reaction(reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict)
			print("rule_equation", rule_equation)
			self.create_rule("assignment", {"equation":rule_equation})





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
	cdef void compute_lineage_propensities(self, double *state, double *propensity_destination, double volume, double time):
		for ind in range(self.num_reactions):
			propensity_destination[ind] = (<Propensity> (self.c_propensities[0][ind])).get_stochastic_volume_propensity(state, self.c_param_values, volume, time)
		for ind in range(self.num_lineage_propensities):
			propensity_destination[self.num_reactions+ind] = (<Propensity> (self.c_lineage_propensities[0][ind])).get_stochastic_volume_propensity(state, self.c_param_values, volume, time)

	#Applies all rules that change the volume of a cell
	cdef double apply_volume_rules(self, double *state, double volume, double time, double dt):
		for ind in range(self.num_volume_rules):
			volume = (<VolumeRule>self.c_volume_rules[0][ind]).get_volume(state, self.c_param_values, volume, time, dt)
		return volume

	#Applies death rules in the order they were added to the model. Returns the index of the first death rule that returns True. -1 otherwise.
	cdef int apply_death_rules(self, double *state, double volume, double time, double start_volume, double start_time):
		cdef unsigned isdead = 0
		for ind in range(self.num_death_rules):
			isdead = (<DeathRule>self.c_death_rules[0][ind]).check_dead(state, self.c_param_values, time, volume, start_time, start_volume)
			if isdead > 0:
				return <int>ind
		return -1

	#Applies divison rules in the order they were added to the model. Returns the index of the first division rule that returns True. -1 otherwise
	cdef int apply_division_rules(self, double *state, double volume, double time, double start_volume, double start_time):
		cdef unsigned divided = 0
		for ind in range(self.num_division_rules):
			divided = (<DivisionRule>self.c_division_rules[0][ind]).check_divide(state, self.c_param_values, time, volume, start_time, start_volume)
			if divided > 0:
				#print("cell divided with time=", time, "volume=", volume, "start_time=", start_time, "start_volume=", start_volume)
				return <int>ind
		return -1

	#Applies a single volume event, determined by the index passed in
	cdef double apply_volume_event(self, unsigned event_index, double *state, double current_time, double current_volume):
		current_volume = (<VolumeEvent>self.c_volume_events[0][event_index]).get_volume(state, self.c_param_values, current_volume, current_time)
		return current_volume

	#Divides a single cell using a VolumeSplitter determined by vsplit_index
	cdef np.ndarray partition(self, unsigned vsplit_ind, LineageVolumeCellState parent):
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
		return super().py_set_state(state)

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

	cdef VolumeCellState get_final_cell_state(self):
		cdef unsigned final_index = (<np.ndarray[np.double_t,ndim=1]> self.timepoints).shape[0]-1
		cdef LineageVolumeCellState cs  = LineageVolumeCellState(t0 = self.timepoints[0], v0 = self.volume[0], state = self.simulation_result[final_index,:])
		cs.set_time(self.timepoints[final_index])
		cs.set_volume(self.volume[final_index])
		cs.set_divided(self.divided)
		cs.set_dead(self.dead)
		return cs


cdef class LineageVolumeSplitter(VolumeSplitter):
	cdef unsigned how_to_split_v
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

		#Figure out how volume will be split
		if "volume" not in options or options["volume"] == "binomial":
			self.how_to_split_v = 0
		elif options["volume"] == "duplicate":
			self.how_to_split_v = 1
		elif options["volume"] == "perfect":
			self.how_to_split_v = 2
		elif options["volume"] in custom_partition_functions:
			self.how_to_split_v = 3
			self.ind2customsplitter["volume"] = options["volume"]
		else:
			raise ValueError("Custom partition function key, "+str(options["volume"])+", for 'volume' not in custom_partition_functions")

		#Figure out how other species are split
		for s in M.get_species2index():
			index = M.get_species_index(s)
			if s not in options or options[s] == "binomial":
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


		# partition the states, copying already takes care of duplication replications.
		cdef np.ndarray[np.double_t,ndim=1] dstate = parent.get_state().copy()
		cdef np.ndarray[np.double_t,ndim=1] estate = parent.get_state().copy()
		cdef unsigned length = dstate.shape[0]

		cdef unsigned loop_index = 0
		cdef unsigned species_index = 0
		cdef unsigned amount = 0
		cdef unsigned amount2 = 0
		cdef double d_value = 0.0

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
		cdef LineageVolumeCellState d = LineageVolumeCellState(v0 = v0d, t0 = t0, state = dstate)
		cdef LineageVolumeCellState e = LineageVolumeCellState(v0 = v0e, t0 = t0, state = estate)
		ans[0] = d
		ans[1] = e

		return ans 

cdef class LineageSSASimulator():
	#SSA for a single cell. Simulates until it devides or dies using division / death rules and/or reactions.
	cdef SingleCellSSAResult SimulateSingleCell(self, LineageCSimInterface sim, LineageVolumeCellState v, np.ndarray timepoints):
		cdef np.ndarray[np.double_t,ndim=1] c_timepoints = timepoints.copy()

		cdef np.ndarray[np.double_t,ndim=1] c_current_state
		cdef np.ndarray[np.double_t,ndim=2] c_stoich = sim.get_update_array() + sim.get_delay_update_array()
		
		cdef unsigned num_species = c_stoich.shape[0]

		cdef unsigned num_reactions = c_stoich.shape[1]
		cdef unsigned num_volume_events = sim.get_num_volume_events()
		cdef unsigned num_death_events = sim.get_num_death_events()
		cdef unsigned num_division_events = sim.get_num_division_events()
		cdef unsigned num_volume_rules = sim.get_num_volume_rules()
		cdef unsigned num_death_rules = sim.get_num_death_rules()
		cdef unsigned num_division_rules = sim.get_num_division_rules()
		cdef unsigned num_propensities = num_reactions + num_volume_events + num_death_events + num_division_events
		cdef unsigned num_timepoints = len(timepoints)

		cdef double initial_time = v.get_initial_time()
		cdef double current_time = v.get_time()
		cdef double final_time = c_timepoints[num_timepoints-1]
		cdef double proposed_time = 0.0
		cdef double Lambda = 0.0
		cdef np.ndarray[np.double_t,ndim=2] c_results = np.zeros((num_timepoints,num_species))
		cdef np.ndarray[np.double_t,ndim=1] c_propensity = np.zeros(num_propensities)
		cdef np.ndarray[np.double_t,ndim=1] c_volume_trace = np.zeros(num_timepoints,)

		cdef unsigned current_index = 0
		cdef unsigned reaction_choice = 4294967295 # https://en.wikipedia.org/wiki/4,294,967,295
		cdef unsigned species_index = 4294967295

		cdef double delta_t = c_timepoints[1]-c_timepoints[0]
		cdef double next_queue_time = c_timepoints[current_index+1]
		cdef unsigned move_to_queued_time = 0

		initial_volume = v.get_initial_volume()
		cdef double current_volume = v.get_volume()
		cdef int cell_divided = -1
		cdef int cell_dead = -1

		cdef SingleCellSSAResult SCR
		cdef VolumeCellState d1, d2
		
		#Set Initial State
		if v.get_state_set() == 1:
			c_current_state = v.py_get_state().copy()
		else:
			warnings.warn("No initial state set (via LineageVolumeCellState v) in SingleCellSSAResuslt. Defaulting to the Model's initial state.")
			c_current_state = sim.get_initial_state().copy()
			v.py_set_state(c_current_state)


		#Warn user if delays are in the model (which will be converted to non-delay reactions)
		if (sim.py_get_delay_update_array() != np.zeros(sim.py_get_delay_update_array().shape)).any():
			warnings.warn("delay reactions found in the model. SingleCellSSASimulator will simulate these reactions without delay. To use delays, please use the SingleCellDelaySSASimulator.")

		# Do the SSA part now
		while current_index < num_timepoints:
			# Compute rules in place
			sim.apply_repeated_rules(<double*> c_current_state.data, current_time)

			#returns the index of the first DeathRule that returned True and -1 otherwise
			cell_dead = sim.apply_death_rules(<double*> c_current_state.data, current_volume, current_time, initial_volume, initial_time)
			#returns the index of the first DivisionRule that returned True and -1 otherwise
			cell_divided = sim.apply_division_rules(<double*> c_current_state.data, current_volume, current_time, initial_volume, initial_time)
			
			#Break the loop if we are at a queued time
			if cell_dead >= 0 and cell_divided >= 0:
				warnings.warn("Cell Death and Division Occured Simultaneously - Death Takes Precedent")
				cell_divided = -1
				break
			elif cell_dead >= 0:
				break
			elif cell_divided >= 0:
				break
			#Compute Reaction and Event propensities in-place
			sim.compute_lineage_propensities(<double*> (c_current_state.data), <double*> (c_propensity.data), current_volume, current_time)

			Lambda = cyrandom.array_sum(<double*> (c_propensity.data), num_propensities)
			#print("propensities computed")
			# Either we are going to move to the next queued time, or we move to the next reaction time.
			
			if Lambda == 0:
				proposed_time = final_time+1
			else:
				proposed_time = current_time + cyrandom.exponential_rv(Lambda)
			if next_queue_time < proposed_time:
				current_time = next_queue_time
				next_queue_time += delta_t
				move_to_queued_time = 1
			else:
				current_time = proposed_time
				move_to_queued_time = 0
			v.set_time(current_time)

			# Update the results array with the state for the time period that we just jumped through.
			while current_index < num_timepoints and c_timepoints[current_index] <= current_time:
				for species_index in range(num_species):
					c_results[current_index,species_index] = c_current_state[species_index]
				c_volume_trace[current_index] = current_volume
				current_index += 1

			# Now update the state accordingly.
			# IF the queue won, then update the volume and continue on or stop if the cell divided.
			if move_to_queued_time == 1:
				# Update the volume every dtyp
				current_volume = sim.apply_volume_rules(<double*>(c_current_state.data), current_volume, current_time, delta_t)
				v.set_volume(current_volume)

			# if an actual reaction happened, do the reaction and maybe update the queue as well.
			else:
				# select a reaction
				#print("num_propensities", num_propensities, "c_propensity", c_propensity, "lambda", Lambda )
				reaction_choice = cyrandom.sample_discrete(num_propensities, <double*> c_propensity.data, Lambda )

				#Propensities are Ordered:
				# Reactions, Divison Events, Volume Events, Death Events
				#print("rxn choince =", reaction_choice)
				#Propensity is a reaction
				if reaction_choice < num_reactions:
					# Do the reaction's initial stoichiometry.
					for species_index in range(num_species):
						c_current_state[species_index] += c_stoich[species_index,reaction_choice]
				#Propensity is a VolumeEvent
				elif reaction_choice >= num_reactions and reaction_choice < num_reactions + num_volume_events:
					current_volume = sim.apply_volume_event(reaction_choice - num_reactions, <double*>(c_current_state.data), current_time, current_volume)
					v.set_volume(current_volume)
				#Propensity is a DivisionEvent.
				elif reaction_choice >= num_reactions+num_volume_events and reaction_choice < num_reactions + num_volume_events+num_division_events:
					#Cell Divided = DivisionEvent Index + num_division_rules
					cell_divided = reaction_choice - num_reactions - num_volume_events + num_division_rules
					break
				#Propensity is a Death Event
				elif reaction_choice >= num_reactions + num_volume_events+num_division_events:
					#Cell Divided = DeathEvent Index + num_death_rules
					cell_dead = reaction_choice - num_reactions + num_volume_events+num_division_events+num_death_rules
					break
				else:
					raise ValueError("More reaction propensities than expected!")


		if cell_divided>=0 or cell_dead>=0:
			#Push current state to the nearest index
			if current_time < c_timepoints[current_index]:
				for species_index in range(num_species):
					c_results[current_index,species_index] = c_current_state[species_index]
				c_volume_trace[current_index] = current_volume
				current_index += 1
			c_timepoints = c_timepoints[:(current_index)]
			c_volume_trace = c_volume_trace[:(current_index)]
			c_results = c_results[:current_index,:]

		#vsr (SingleCellSSAResult) contains the simulation results until cell death / division or simualtion termination.
		#cell_divided and cell_dead are returend via vsr so the events/rules/VolumeSplitters can be called by the lineage simualtion loop.

		SCR = SingleCellSSAResult(c_timepoints,c_results,c_volume_trace, cell_divided >= 0)
		SCR.set_divided(cell_divided)
		SCR.set_dead(cell_dead)
		SCR.set_volume_object(v.get_volume_object())
		
		return SCR

	def py_SimulateSingleCell(self, np.ndarray timepoints, LineageModel Model = None, LineageCSimInterface interface = None, LineageVolumeCellState v = None):
		if Model == None and interface == None:
			raise ValueError('py_SimulateSingleCell requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
		elif interface == None:
			interface = LineageCSimInterface(Model)
			interface.py_set_initial_time(timepoints[0])
		if v == None:
			v = LineageVolumeCellState(v0 = 1, t0 = 0, state = Model.get_species_array())

		return self.SimulateSingleCell(interface, v, timepoints)

def py_SimulateSingleCell(np.ndarray timepoints, LineageModel Model = None, LineageCSimInterface interface = None, LineageVolumeCellState v = None):
	simulator = LineageSSASimulator()
	result = simulator.py_SimulateSingleCell(timepoints, Model = Model, interface = interface, v = v)
	return result

cdef Lineage SimulateCellLineage(LineageCSimInterface sim, list initial_cell_states, np.ndarray timepoints, LineageSSASimulator simulator):
	# Prepare a lineage structure to save the data output.
	cdef Lineage l = Lineage()
	cdef np.ndarray[np.double_t, ndim=1] c_timepoints = timepoints
	cdef np.ndarray[np.double_t, ndim=1] c_truncated_timepoints
	cdef double final_time = c_timepoints[c_timepoints.shape[0]-1]

	cdef unsigned i
	cdef unsigned list_index = 0
	cdef unsigned debug_index
	cdef list old_schnitzes = []
	cdef list old_cell_states = []

	cdef SingleCellSSAResult r
	cdef Schnitz s
	cdef LineageVolumeCellState cs
	cdef LineageVolumeCellState dcs
	cdef Schnitz daughter_schnitz1, daughter_schnitz2
	cdef LineageVolumeCellState d1, d2, d1final, d2final
	cdef np.ndarray daughter_cells

	# Simulate the first cell until death division or max time
	for i in range(len(initial_cell_states)):
		r = simulator.SimulateSingleCell(sim, initial_cell_states[i], c_timepoints)
		s = r.get_schnitz()
		s.set_parent(None)
		l.add_schnitz(s)
		old_schnitzes.append(s)
		old_cell_states.append(r.get_final_cell_state())

	#print('Entering Simulatin Queue')
	while list_index < len(old_cell_states):
		cs = old_cell_states[list_index]
		s = old_schnitzes[list_index]
		list_index += 1

		#If the cell has already simulated all its time, do nothing
		if cs.get_time() >= final_time- 1E-9:
			pass
		#If the cell is dead, do nothing
		elif cs.get_dead() >= 0:
			pass
		#If the cell has divided, apply the appropriate division rule
		elif cs.get_divided() >= 0:
			daughter_cells = sim.partition(cs.get_divided(), cs)
			
			d1 = <LineageVolumeCellState>(daughter_cells[0])
			d2 = <LineageVolumeCellState>(daughter_cells[1])

			#Create a new timepoint array and simulate the first daughter and queue if it doesn't reach final time.
			c_truncated_timepoints = c_timepoints[c_timepoints > cs.get_time()]
			r = simulator.SimulateSingleCell(sim, d1, c_truncated_timepoints)
			daughter_schnitz1 = r.get_schnitz()
			d1final = r.get_final_cell_state()

			# Add on the new daughter if final time wasn't reached.
			if d1final.get_time() < final_time + 1E-9:
				old_schnitzes.append(daughter_schnitz1)
				old_cell_states.append(d1final)
			else:
				warnings.warn("Daughter cell simulation went over the total time. Simulation has been discarded. Check for model errors.")

			c_truncated_timepoints = c_timepoints[c_timepoints > cs.get_time()]
			r = simulator.SimulateSingleCell(sim, d2, c_truncated_timepoints)
			daughter_schnitz2 = r.get_schnitz()
			d2final = r.get_final_cell_state()

			if d2final.get_time() < final_time + 1E-9:
				old_schnitzes.append(daughter_schnitz2)
				old_cell_states.append(d2final)
			else:
				warnings.warn("Daughter cell simulation went over the total time. Simulation has been discarded. Check for model errors.")

			# Set up daughters and parent appropriately.
			daughter_schnitz1.set_parent(s)
			daughter_schnitz2.set_parent(s)
			s.set_daughters(daughter_schnitz1,daughter_schnitz2)

			# Add daughters to the lineage
			l.add_schnitz(daughter_schnitz1)
			l.add_schnitz(daughter_schnitz2)

	return l

def py_SimulateCellLineage(timepoints, initial_cell_states = [], LineageCSimInterface interface = None, LineageModel Model = None, LineageSSASimulator simulator = None):
	if Model == None and interface == None:
		raise ValueError('py_SimulateCellLineage requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
	elif interface == None:
		interface = LineageCSimInterface(Model)
		interface.py_set_initial_time(timepoints[0])
	if len(initial_cell_states) == 0:
		initial_cell_states = [LineageVolumeCellState(v0 = 1, t0 = 0)]
	if simulator == None:
		simulator = LineageSSASimulator()

	return SimulateCellLineage(interface, initial_cell_states, timepoints, simulator)

#Simulates an ensemble of cells over some amount of time.
#Returns the final cell-states of all cells (and their offspring).
#  dead cells are included based upon the include_dead_cells parameter (1 = included, otherwise = excluded)
cdef list PropagateCells(LineageCSimInterface sim, list initial_cell_states, np.ndarray timepoints, unsigned include_dead_cells, LineageSSASimulator simulator):
	cdef np.ndarray c_timepoints = timepoints
	cdef np.ndarray c_truncated_timepoints
	cdef np.ndarray index_list
	cdef list final_cell_states = []
	cdef double final_time = c_timepoints[timepoints.shape[0]-1]
	cdef unsigned list_index = 0
	cdef list old_cell_states = list(initial_cell_states)
	cdef LineageVolumeCellState cs

	for i in range(len(initial_cell_states)):
		cs = simulator.SimulateSingleCell(sim, initial_cell_states[i], c_timepoints).get_final_cell_state()
		old_cell_states.append(cs)

	while list_index < len(old_cell_states):
		cs = old_cell_states[list_index]
		list_index += 1

		#If the cell has already simulated all its time, do nothing
		if cs.get_time() >= final_time- 1E-9:
			final_cell_states.append(cs)
		#If the cell is dead, do nothing
		elif cs.get_dead() >= 0:
			if include_dead_cells == 1:
				final_cell_states.append(cs)
		#If the cell has divided, apply the appropriate division rule
		elif cs.get_divided() >= 0:
			daughter_cells = sim.partition(cs.get_divided(), cs)
			
			d1 = <LineageVolumeCellState>(daughter_cells[0])
			d2 = <LineageVolumeCellState>(daughter_cells[1])

			#Create a new timepoint array and simulate the first daughter and queue if it doesn't reach final time.
			c_truncated_timepoints = c_timepoints[c_timepoints > cs.get_time()]
			d1final = simulator.SimulateSingleCell(sim, d1, c_truncated_timepoints).get_final_cell_state()

			# Add on the new daughter if final time wasn't reached.
			if d1final.get_time() < final_time + 1E-9:
				old_cell_states.append(d1final)
			else:
				warnings.warn("Daughter cell simulation went over the total time. Simulation has been discarded. Check for model errors.")

			c_truncated_timepoints = c_timepoints[c_timepoints > cs.get_initial_time()]
			d2final = simulator.SimulateSingleCell(sim, d2, c_truncated_timepoints).get_final_cell_state()

			if d2final.get_time() < final_time + 1E-9:
				old_cell_states.append(d2final)
			else:
				warnings.warn("Daughter cell simulation went over the total time. Simulation has been discarded. Check for model errors.")

	return final_cell_states


def  py_PropagateCells(timepoints, initial_cell_states = [], LineageModel Model = None, LineageCSimInterface interface = None, LineageSSASimulator simulator = None, include_dead_cells = 0):

	if Model == None and interface == None:
		raise ValueError('py_PropagateCells requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
	elif interface == None:
		interface = LineageCSimInterface(Model)
		interface.py_set_initial_time(timepoints[0])
	if len(initial_cell_states) == 0:
		initial_cell_states = [LineageVolumeCellState(v0 = 1, t0 = 0, state = Model.get_species_array())]
	if simulator == None:
		simulator = LineageSSASimulator()

	return PropagateCells(interface, initial_cell_states, timepoints, include_dead_cells, simulator)

#Propogates a single cell trajectory, ignoring half the daughters every division.
cdef list SingleCellLineage(LineageCSimInterface sim, LineageVolumeCellState initial_cell, np.ndarray timepoints, LineageSSASimulator simulator):
	cdef np.ndarray c_timepoints = timepoints
	cdef np.ndarray c_truncated_timepoints
	cdef np.ndarray index_list
	cdef list old_cell_states = []
	cdef list single_cell_trajectory = []
	cdef double final_time = c_timepoints[timepoints.shape[0]-1]
	cdef unsigned list_index = 0
	cdef unsigned ind
	cdef SingleCellSSAResult r
	cdef LineageVolumeCellState cs

	r = simulator.SimulateSingleCell(sim,  initial_cell, c_timepoints)
	old_cell_states.append(r.get_final_cell_state())
	single_cell_trajectory.append(r)

	while list_index < len(old_cell_states):
		cs = old_cell_states[list_index]
		list_index += 1

		#If the cell has already simulated all its time, do nothing
		if cs.get_time() >= final_time- 1E-9:
			pass
		#If the cell is dead, do nothing
		elif cs.get_dead() >= 0:
			pass
		#If the cell has divided, apply the appropriate division rule
		elif cs.get_divided() >= 0:
			daughter_cells = sim.partition(cs.get_divided(), cs)

			#Create a new timepoint array and simulate a random daughter and queue if it doesn't reach final time.
			ind = <unsigned>cyrandom.uniform_rv()>.5
			d = <LineageVolumeCellState>(daughter_cells[ind])
			c_truncated_timepoints = c_timepoints[c_timepoints > cs.get_time()]
			r = simulator.SimulateSingleCell(sim, d, c_truncated_timepoints)
			single_cell_trajectory.append(r)
			dfinal = r.get_final_cell_state()

			# Add on the new daughter if final time wasn't reached.
			if dfinal.get_time() < final_time + 1E-9:
				old_cell_states.append(dfinal)
				single_cell_trajectory.append(r)
			else:
				warnings.warn("Daughter cell simulation went over the total time. Simulation has been discarded. Check for model errors.")
	return single_cell_trajectory

def py_SingleCellLineage(np.ndarray timepoints, LineageVolumeCellState initial_cell = None,LineageModel Model = None, LineageCSimInterface interface = None, LineageSSASimulator simulator = None):
	if Model == None and interface == None:
		raise ValueError('py_SingleCellPropogation requires either a LineageModel Model or a LineageCSimInterface interface to be passed in as keyword parameters.')
	elif interface == None:
		interface = LineageCSimInterface(Model)
		interface.py_set_initial_time(timepoints[0])
	if initial_cell == None:
		initial_cell = LineageVolumeCellState(v0 = 1, t0 = 0, state = Model.get_species_array())
	if simulator == None:
		simulator = LineageSSASimulator()
	return SingleCellLineage(interface, initial_cell, timepoints, simulator)


#Inputs:
#  initial_cells: a list of LineageVolumeCellStates
#  sims: a list of LineageCsimInterfaces.
#  interface_inds is a mapping of from the initial_cell index --> LineageCSimInterface index [in the sims list]
#  global_sync_period: global species are synchronized between the cell ensemble every global_sync_period
#  global_species_indices: an (# global species x # interfaces) array of species indices for each model (interface)
# Returns a Lineage

"""How it works:
* cells are linked to an interface (to allow different cell types) via the interface_inds array:
** cell_state[i] --> sim_interfaces[interface_inds[i]]
** global_species_inds[i, j] corresponds to the ith global species index for the jth interface (model) which allows for different models to have different indices

* every T=global_sync_period:
** each living cell
** each cell is simulated by SimulateCellLineage with its appropriate interface
*** new living cell states are extracted from the returned Lineage
** global species are synchronized between living cells via the global_species_inds array:
***  each global species g_i is updated by all the cells: g_i <-- g_i + sum_cells [Change in g_i over cells life]
** returned lineages are appended to the lineage

	"""

cdef Lineage SimulateInteractingCellLineage(list sim_interfaces, list initial_cell_states, list interface_inds, np.ndarray timepoints, LineageSSASimulator simulator, double global_sync_period, np.ndarray global_species_inds):
	
	# Prepare a lineage structure to save the data output.
	cdef Lineage l = Lineage() #stores the final output lineages
	cdef np.ndarray[np.double_t, ndim=1] c_timepoints = timepoints #stores all timepoints
	cdef np.ndarray[np.double_t, ndim=1] c_period_timepoints #stores the time points for the given sim period
	cdef np.ndarray[np.double_t, ndim=1] c_truncated_timepoints #stores a truncated set of timepoints
	cdef double dt = c_timepoints[1] - c_timepoints[0]
	cdef double final_time = c_timepoints[c_timepoints.shape[0]-1] #when the entire simulation ends (adding dt onto the end for rounding reasons)
	cdef double period_time = c_timepoints[0]+global_sync_period
	cdef double current_time = c_timepoints[0]
	cdef LineageCSimInterface sim #Stores the LineageCSimInterface passed into SimulateCellLineage between interaction syncs
	cdef unsigned num_global_species = global_species_inds.shape[0]
	cdef np.ndarray[np.double_t, ndim=1] global_species = np.zeros(num_global_species)
	cdef np.ndarray[np.double_t, ndim=1] delta_global_species = np.zeros(num_global_species)
	cdef unsigned i = 0
	cdef unsigned j = 0
	cdef int spec_ind = 0
	cdef unsigned list_index = 0

	#Used for lineage simulation
	cdef list old_interface_inds = [] #stores the index of the interface for each simulation
	cdef list new_interface_inds = [] #stores the index of the interface for the next round of simulation
	cdef list old_schnitzes = [] #stores schnitzes that have been added to the lineage
	cdef list new_schnitzes = [] #stores schnitzes for merging in the next round of simulation
	cdef list old_cell_states = [] #stores cell states which have been simulated
	cdef list new_cell_states = [] #stores cell states to begin the next round of simulation
	cdef SingleCellSSAResult r, new_r, merge_r
	cdef Schnitz s, new_s, daughter_schnitz1, daughter_schnitz2
	cdef LineageVolumeCellState cs, new_cs, d1, d2, d1final, d2final
	cdef np.ndarray daughter_cells
	cdef np.ndarray state

	#Check that len sims == len initial_cells. As cells divide, they inherit their mother's CSimInterface
	if len(interface_inds) != len(initial_cell_states):
		raise ValueError("sims [a list of LineageCSimInterfaces] must be the same length as initial cells [a list of LineageVolumeCellStates].")
	if max(interface_inds) >= len(sim_interfaces):
		raise ValueError("interface_inds index too large. Must contain valid to the sims (LineagesCSimInterface list) argument.")

	#Check that global period is greater than dt
	if global_sync_period < dt:
		raise ValueError("global sync period must be larger than the timestep, dt, in the timepoints passed in for simulation.")
	#print("lineage sim from", current_time, "to", final_time)

	# Simulate the initial cells until death division or synch-time
	#print('Simulating Initial Cells from', current_time, "to", period_time)
	c_period_timepoints = c_timepoints[c_timepoints < period_time+1E-9]
	#print("c_period_timepoints:", c_period_timepoints[0], c_period_timepoints[len(c_period_timepoints)-1])
	#print("next timepoint:", c_timepoints[len(c_period_timepoints)])
	for list_index in range(len(initial_cell_states)):
		#print("simulating cell", list_index)
		i = interface_inds[list_index]
		sim = sim_interfaces[i]
		r = simulator.SimulateSingleCell(sim, initial_cell_states[list_index], c_period_timepoints)
		s = r.get_schnitz()
		s.set_parent(None)
		cs = r.get_final_cell_state()

		if cs.get_time() < c_timepoints[len(c_period_timepoints)] and cs.get_time() >= c_timepoints[len(c_period_timepoints)-1]:
			#print("cell requires next period")
			new_schnitzes.append(s)
			new_cell_states.append(cs)
			new_interface_inds.append(i)
		else:
			#print("cell requires a continuation simulation: cs.get_time()=", cs.get_time(), "period_time=", period_time)
			old_schnitzes.append(s)
			old_cell_states.append(cs)
			old_interface_inds.append(i)

	#Global Simulation Loop till final_time
	#print("Entering Queue Loop")
	while period_time <= final_time+dt+1E-8:
		#print("simulating from", current_time, "to", period_time)
		c_period_timepoints = c_timepoints[c_timepoints < period_time]
		current_time = period_time
		#Update the period time
		#If the total sim time isn't a multiple of the period time, update carefully
		if period_time + global_sync_period > final_time+dt and period_time < final_time:
			period_time = final_time+dt
		#Otherwise update normally
		else:
			period_time = period_time + global_sync_period
		
		
		#print("****Synchronizing Global Species!*****")
		#print("initial global_speices = ", global_species)
		for i in range(num_global_species):
			#Calculate the global change in species
			delta_global_species[i] = 0
			for list_index in range(len(old_cell_states)):
				j = old_interface_inds[list_index]
				spec_ind = global_species_inds[i, j]
				cs = old_cell_states[list_index]
				state = cs.get_state()
				#print('Cell', list_index, "global_species", i, "=", state[spec_ind])
				if spec_ind >= 0:
					delta_global_species[i] = delta_global_species[i] + state[spec_ind] - global_species[i]
			
			#print("delta_global_species", delta_global_species)
			if delta_global_species[i] >= -global_species[i]:
				global_species[i] = global_species[i] + delta_global_species[i]
			else:
				global_species[i] = 0

			#Synchronize all the cells
			for list_index in range(len(old_cell_states)):
				j = old_interface_inds[list_index]
				spec_ind = global_species_inds[i, j]
				cs = old_cell_states[list_index]
				state = cs.get_state()
				if spec_ind >= 0:
					cs.set_state_comp(global_species[i], spec_ind)
		#print("final global_speices = ", global_species)
		#enter simulation queue
		list_index = 0
		while list_index < len(old_cell_states):
			#print("Simulating list index:", list_index, "len(old_cell_states)", len(old_cell_states))
			cs = old_cell_states[list_index]
			s = old_schnitzes[list_index]
			i = old_interface_inds[list_index]
			sim = sim_interfaces[i]
			list_index += 1

			#If the cell is dead or it has simulated to the final time, add it to the lineage
			if cs.get_dead() >= 0 or cs.get_time() >= final_time - 1E-9:
				#print("simulation complete")
				# Add daughters to the lineage
				l.add_schnitz(s)
			#If the cell has divided, apply the appropriate division rule
			elif cs.get_divided() >= 0:
				#print('cell divided - adding daughters')

				c_truncated_timepoints = c_period_timepoints[c_period_timepoints > cs.get_time()]
				#sometimes a cell divides at the very end of a period - push it to the next period's queue.
				if len(c_truncated_timepoints) <= 1:
					new_schnitzes.append(s)
					new_interface_inds.append(i)
					new_cell_states.append(cs)

				#Otherwise, divide the cell
				else:
					#Add mother to the lineage
					l.add_schnitz(s)

					daughter_cells = sim.partition(cs.get_divided(), cs)
					
					d1 = <LineageVolumeCellState>(daughter_cells[0])
					d2 = <LineageVolumeCellState>(daughter_cells[1])

					#Create a new timepoint array and simulate the first daughter
					
					#print("c_truncated_timepoints", c_truncated_timepoints[0], c_truncated_timepoints[len(c_truncated_timepoints)-1])
					r = simulator.SimulateSingleCell(sim, d1, c_truncated_timepoints)
					daughter_schnitz1 = r.get_schnitz()
					d1final = r.get_final_cell_state()

					#simulate the second daughter
					c_truncated_timepoints = c_period_timepoints[c_period_timepoints > cs.get_time()]
					r = simulator.SimulateSingleCell(sim, d2, c_truncated_timepoints)
					daughter_schnitz2 = r.get_schnitz()
					d2final = r.get_final_cell_state()

					# Set up daughters and parent appropriately.
					daughter_schnitz1.set_parent(s)
					daughter_schnitz2.set_parent(s)
					s.set_daughters(daughter_schnitz1,daughter_schnitz2)

					# Add the new daughter1 cell to the next periods queue if the final period time hasn't been reached
					if d1final.get_time() < c_timepoints[len(c_period_timepoints)] and d1final.get_time() >= c_timepoints[len(c_period_timepoints)-1]:
						new_schnitzes.append(daughter_schnitz1)
						new_cell_states.append(d1final)
						new_interface_inds.append(i)
					# Otherwise continue in the same queue
					else:
						old_schnitzes.append(daughter_schnitz1)
						old_cell_states.append(d1final)
						old_interface_inds.append(i)

					# Add the new daughter2 cell to the next periods queue if the final period time hasn't been reached
					if d2final.get_time() < c_timepoints[len(c_period_timepoints)] and d2final.get_time() >= c_timepoints[len(c_period_timepoints)-1]:
						new_schnitzes.append(daughter_schnitz2)
						new_cell_states.append(d2final)
						new_interface_inds.append(i)
					# Otherwise continue in the same queue
					else:
						old_schnitzes.append(daughter_schnitz2)
						old_cell_states.append(d2final)
						old_interface_inds.append(i)

			#If the cell isn't dead or divided, simulate it more
			else:
				c_truncated_timepoints = c_period_timepoints[c_period_timepoints > cs.get_time()]
				#print("period_Time", period_time, "len(c_period_timepoints)", len(c_period_timepoints), "c_period_timepoints[0]=", c_period_timepoints[0], "c_period_timepoints[-1]=", c_period_timepoints[len(c_period_timepoints)-1])
				#print("cs.get_time()=", cs.get_time(), "len(c_truncated_timepoints)", len(c_truncated_timepoints), "c_truncated_timepoints[0]=", c_truncated_timepoints[0], "c_truncated_timepoints[-1]=", c_truncated_timepoints[len(c_truncated_timepoints)-1])

				#If there is only one timepoint left, push to the next period.
				if len(c_truncated_timepoints) <= 1:
					#print("only one timepoint left, push by deltaT")
					new_schnitzes.append(s)
					new_cell_states.append(cs)
					new_interface_inds.append(i)
				else:
					#print("continuation simulation")
					new_r = simulator.SimulateSingleCell(sim, cs, c_truncated_timepoints)
					#print("SSA Complete")
					new_cs = new_r.get_final_cell_state()
					new_cs.set_initial_vars(cs.get_initial_volume(), cs.get_initial_time())
					new_cs.set_time(cs.get_time())
					new_cs.set_volume(cs.get_volume())

					merge_r = SingleCellSSAResult(
						np.concatenate((s.get_time(), new_r.get_timepoints())), 
						np.concatenate((s.get_data(), new_r.get_result())), 
						np.concatenate((s.get_volume(), new_r.get_volume())), 
						new_r.get_divided() >= 0)

					merge_r.set_divided(new_r.get_divided())
					merge_r.set_dead(new_r.get_dead())
					new_s = merge_r.get_schnitz()
					new_s.set_parent(s.get_parent())
					
					#print("ssa results merged")

					#print("adding to lists")
					#Save final Schnitz
					if new_cs.get_time() >= final_time:
						#print("Saving Final Schnitz", "new_cs.get_time()", new_cs.get_time())
						l.add_schnitz(new_s)
					#move to the next time period
					elif new_cs.get_time() <= period_time-1E-8 and new_cs.get_time() >= c_truncated_timepoints[len(c_truncated_timepoints)-1]:
						#print("Add Schnitz to next Period", new_cs.get_time(), period_time)
						new_schnitzes.append(new_s)
						new_cell_states.append(new_cs)
						new_interface_inds.append(i)
					#stay in the same period
					else:
						#print('Adding Schnitz to to current queue', new_cs.get_time(), period_time)
						old_schnitzes.append(new_s)
						old_cell_states.append(new_cs)
						old_interface_inds.append(i)
						
					#print("continuation sim complete")

		#print("while loop complete: period_time=", period_time, "final_time=", final_time)
		#reset lists
		old_cell_states = new_cell_states
		new_cell_states = []
		old_schnitzes = new_schnitzes
		new_schnitzes = []
		old_interface_inds = new_interface_inds
		new_interface_inds = []


	return l


def py_SimulateInteractingCellLineage(timepoints, global_sync_period, global_species = [], sim_interfaces = [], models = [], initial_cell_states = [], interface_inds = [], simulator = None, global_species_inds = None):
	if simulator == None:
		simulator = LineageSSASimulator()
	if len(models) == 0 and len(sim_interfaces) == 0:
		raise ValueError("Missing Required Keyword Arguments:models = [LineageModel] or sim_interfaces = [LineageCSimInterface]")
	elif len(sim_interfaces) == 0:
		sim_interfaces = [LineageCSimInterface(m) for m in models]
		global_species_inds = np.zeros((len(global_species), len(models)))
		for i in range(len(global_species)):
			for j in range(len(models)):
				m = models[j]
				s = global_species[i]
				ind = m.get_species_index(s)
				if ind != None:
					global_species_inds[i, j] = ind
				else:
					global_species_inds[i, j] = -1
	elif len(models) != 0 and len(sim_interfaces) != 0:
		raise ValueError("Must call py_SimulateInteractingCellLineage with either the keyword argument models or the keyword argument sim_interfaces, not both.")
	elif len(models) == 0 and (global_species_inds == None or global_species_inds.shape[1] != len(sim_interfaces)):
		raise ValueError("When calling py_SimulateInteractingCellLineage with the keyword argument sim_interfaces, the argument global_species_inds is required where global_species_inds[i, j] corresponds the species index of the ith global species in the jth interface.")
	elif len(global_species) == 0 and global_species_inds == None:
		warnings.warn('Calling SimulateInteractintCellLineage without any global species defined. Use the global_species or global_species_inds keywords.')
		global_species_inds = np.array()


	if len(initial_cell_states) == len(models) and (not isinstance(initial_cell_states[0], VolumeCellState)):
		initial_cell_counts = initial_cell_states
		initial_cell_states = []
		interface_inds = []
		for i in range(len(models)):
			M = models[i]
			initial_state = M.get_species_array()
			for j in range(initial_cell_counts[i]):
				lvcs = LineageVolumeCellState(v0 = 1, t0 = 0, state = initial_state)
				initial_cell_states.append(lvcs)
				interface_inds.append(i)
	elif len(initial_cell_states) == 0 and len(models) > 0:
		warnings.warn("Calling py_SimulateInteractintCellLineage without any initial_cell_states. Defaulting to creating one initial cell for each model.")
		initial_cell_states = [LineageVolumeCellState(v0 = 1, t0 = 0, state = m.get_species_array()) for m in models]
		interface_inds = [i for i in range(len(models))]
	elif len(initial_cell_states) > len(models) and len(interface_inds) != len(initial_cell_states):
		raise ValueError("When passing in more initial_cell_states than models, the keyword argument interface_inds is also required where interface_inds[i] corresponds to the index of the interface/model beloning to initial_cell_state[i].")
	
	return SimulateInteractingCellLineage(sim_interfaces, initial_cell_states, interface_inds, timepoints,simulator, global_sync_period, global_species_inds)

