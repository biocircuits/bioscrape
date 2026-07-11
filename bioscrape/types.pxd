# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

# NOTE: class/method docstrings live in types.pyx, not here. If both
# files define a docstring for the same class/method, the .pyx one
# silently wins with no warning at compile time.

import numpy as np
cimport numpy as np
from vector cimport vector

##################################################                ####################################################
######################################              PROPENSITY TYPES                    ##############################
#################################################                     ################################################


"""
    Enumerated type for each type of propensity function.
    Currently this is used for nothing, but it might come in handy at some point.
"""
ctypedef enum PropensityType:
    unset = -1
    constitutive = 0
    unimolecular = 1
    bimolecular = 2
    general = 3
    hill_positive = 4
    hill_negative = 5
    proportional_hill_positive = 6
    proportional_hill_negative = 7
    mass_action = 8

cdef class Propensity:
    cdef PropensityType propensity_type

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)


cdef class ConstitutivePropensity(Propensity):
    # variables
    cdef unsigned rate_index
    # constructor

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class UnimolecularPropensity(Propensity):

    # variables
    cdef unsigned rate_index
    cdef unsigned species_index


    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)


cdef class BimolecularPropensity(Propensity):

    # variables
    cdef unsigned rate_index
    cdef unsigned s1_index
    cdef unsigned s2_index


    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class PositiveHillPropensity(Propensity):
    # variables
    cdef unsigned K_index
    cdef unsigned rate_index
    cdef unsigned n_index
    cdef unsigned s1_index

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class PositiveProportionalHillPropensity(Propensity):

    # variables
    cdef unsigned K_index
    cdef unsigned rate_index
    cdef unsigned n_index
    cdef unsigned s1_index
    cdef unsigned d_index

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class NegativeHillPropensity(Propensity):

    # variables
    cdef unsigned K_index
    cdef unsigned rate_index
    cdef unsigned n_index
    cdef unsigned s1_index

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class NegativeProportionalHillPropensity(Propensity):

    # variables
    cdef unsigned K_index
    cdef unsigned rate_index
    cdef unsigned n_index
    cdef unsigned s1_index
    cdef unsigned d_index

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class MassActionPropensity(Propensity):
    # variables
    cdef vector[int] sp_inds
    cdef vector[int] sp_counts
    cdef unsigned k_index
    cdef int num_species

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class Term:
    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class ConstantTerm(Term):
    cdef double value

    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class SpeciesTerm(Term):
    cdef unsigned index

    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class ParameterTerm(Term):
    cdef unsigned index

    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class VolumeTerm(Term):
    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

# Putting stuff together

cdef class BinaryTerm(Term):
    cdef vector[void*] terms
    cdef list terms_list

    cdef void add_term(self,Term trm)

    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class SumTerm(BinaryTerm):
    cdef void add_term(self,Term trm)

    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class ProductTerm(BinaryTerm):
    cdef void add_term(self,Term trm)
    cdef double evaluate(self, double *species, double *params, double time)

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class MaxTerm(BinaryTerm):
    cdef void add_term(self,Term trm)
    cdef double evaluate(self, double *species, double *params, double time)

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class MinTerm(BinaryTerm):
    cdef void add_term(self,Term trm)
    cdef double evaluate(self, double *species, double *params, double time)

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class PowerTerm(Term):
    cdef Term base
    cdef Term exponent

    cdef void set_base(self, Term base)
    cdef void set_exponent(self, Term exponent)
    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class ExpTerm(Term):
    cdef Term arg

    cdef void set_arg(self, Term arg)
    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class LogTerm(Term):
    cdef Term arg

    cdef void set_arg(self, Term arg)
    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)


cdef class StepTerm(Term):
    cdef Term arg

    cdef void set_arg(self, Term arg)
    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class AbsTerm(Term):
    cdef Term arg

    cdef void set_arg(self, Term arg)
    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class TimeTerm(Term):
    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class GeneralPropensity(Propensity):
    cdef Term term

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)


##################################################                ####################################################
######################################              DELAY  TYPES                        ##############################
#################################################                     ################################################
"""
enumerated type for types of Delay, currently not used for anything but might come in handy later.
"""
ctypedef enum DelayType:
    unset_delay = -1
    none = 0
    fixed = 1
    gaussian = 2
    gamma = 3


cdef class Delay:
    # delay type
    cdef DelayType delay_type

    # to be overriden
    cdef double get_delay(self, double* state, double* params)

cdef class NoDelay(Delay):

    cdef double get_delay(self, double* state, double* params)

cdef class FixedDelay(Delay):

    cdef unsigned delay_index

    cdef double get_delay(self, double* state, double* params)


cdef class GaussianDelay(Delay):

    cdef unsigned mean_index
    cdef unsigned std_index


    cdef double get_delay(self, double* state, double* params)

cdef class GammaDelay(Delay):
    cdef unsigned k_index
    cdef unsigned theta_index

    cdef double get_delay(self, double* state, double* params)

##################################################                ####################################################
######################################              RULE   TYPES                       ###############################
#################################################                     ################################################

cdef class Rule:
    cdef double frequency_flag #-1 if the rule always repeats, -2 if the rule repeats every dt, otherwise the rule runs at time t == frequency flag.

    cdef void rule_operation(self, double *state, double *params, double time, double dt)
    cdef void rule_volume_operation(self, double *state, double *params, double volume, double time, double dt)
    cdef void execute_rule(self, double *state, double *params, double time, double dt, unsigned rule_step)
    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time, double dt, unsigned rule_step)


cdef class AdditiveAssignmentRule(Rule):
    cdef vector[int] species_source_indices
    cdef unsigned dest_index

    cdef void rule_operation(self, double *state, double *params, double time, double dt)
    cdef void rule_volume_operation(self, double *state, double *params, double volume, double time, double dt)

cdef class GeneralAssignmentRule(Rule):
    cdef Term rhs
    cdef unsigned dest_index
    cdef int param_flag # 1 if the assigned thing is a parameter, 0 if it's a species


    cdef void rule_operation(self, double *state, double *params, double time, double dt)
    cdef void rule_volume_operation(self, double *state, double *params, double volume, double time, double dt)

cdef class GeneralODERule(Rule):
    cdef Term rhs
    cdef unsigned dest_index
    cdef int param_flag # 1 if the assigned thing is a parameter, 0 if it's a species

    cdef void rule_operation(self, double *state, double *params, double time, double dt)
    cdef void rule_volume_operation(self, double *state, double *params, double volume, double time, double dt)



##################################################                ####################################################
######################################              VOLUME  TYPES                       ##############################
#################################################                     ################################################

cdef class Volume:

    cdef double current_volume
    cdef double get_volume_step(self, double *state, double *params, double time, double volume, double dt)
    cdef void initialize(self, double *state, double *params, double time, double volume)
    cdef unsigned cell_divided(self, double *state, double *params, double time, double volume, double dt)
    cdef Volume copy(self)

    cdef inline void set_volume(self, double v):
        """
        Set the volume to some number > 0
        :param v: (double) the volume to set.
        :return: None
        """
        self.current_volume = v

    cdef inline double get_volume(self):
        """
        Get the current volume
        :return: (double) the current volume.
        """
        return self.current_volume


cdef class StochasticTimeThresholdVolume(Volume):
    cdef double division_time
    cdef double average_division_volume
    cdef double division_noise
    cdef double cell_cycle_time
    cdef double growth_rate

    cdef double get_volume_step(self, double *state, double *params, double time, double volume, double dt)
    cdef void initialize(self, double *state, double *params, double time, double volume)
    cdef unsigned cell_divided(self, double *state, double *params, double time, double volume, double dt)
    cdef Volume copy(self)


cdef class StateDependentVolume(Volume):
    cdef double division_volume
    cdef double average_division_volume
    cdef double division_noise
    cdef Term growth_rate

    cdef double get_volume_step(self, double *state, double *params, double time, double volume, double dt)
    cdef void initialize(self, double *state, double *params, double time, double volume)
    cdef unsigned cell_divided(self, double *state, double *params, double time, double volume, double dt)
    cdef Volume copy(self)


##################################################                ####################################################
######################################              MODEL   TYPES                       ##############################
#################################################                     ################################################


cdef class Model:
    ############################################################################
    # DEVELOPER WARNING 
    # 
    # In order to be copiable and usable with multiprocessing, Model must be 
    # picklable. To do that, Model implements a __getstate__ method and a 
    # __setstate__ method, which respectively compress all of the Model's state 
    # variables into a picklable tuple and use those tuples to make a new Model 
    # identical to the old one. 
    # 
    # IF YOU ADD, REMOVE, OR CHANGE ANY VARIABLES HERE, YOU MUST REFLECT THOSE 
    # CHANGES IN THE __getstate__ AND __setstate__ METHODS.
    #
    # This is especially important for newly-added variables. If you add 
    # variables but don't update the pickling methods, then you will introduce 
    # SILENT bugs whenever a user makes a copy of a Model or tries to use a 
    # Model in multiple threads/processes with multiprocessing. 
    ############################################################################
    

    cdef unsigned _next_species_index
    cdef unsigned _next_params_index
    cdef unsigned _dummy_param_counter
    cdef unsigned has_delay

    cdef vector[void*] c_propensities
    cdef list propensities

    cdef list delays
    cdef vector[void*] c_delays

    cdef list repeat_rules
    cdef vector[void*] c_repeat_rules

    cdef dict species2index
    cdef dict params2index
    cdef np.ndarray species_values
    cdef np.ndarray params_values

    cdef np.ndarray update_array
    cdef np.ndarray delay_update_array
    cdef list reaction_list
    cdef list reaction_updates
    cdef list delay_reaction_updates
    cdef int initialized

    cdef (vector[void*])* get_c_propensities(self)
    cdef (vector[void*])* get_c_delays(self)
    cdef (vector[void*])* get_c_repeat_rules(self)

    cdef np.ndarray get_update_array(self)
    cdef np.ndarray get_delay_update_array(self)
    cdef np.ndarray get_species_values(self)
    cdef np.ndarray get_params_values(self)

    cdef void _initialize(self)

    cdef dict txt_dict #Used to store text for writing bioscrape XML
    cdef list reaction_definitions #used to store reaction tuples for SBML
    cdef list rule_definitions #used to store rule tuples for SBML



##################################################                ####################################################
######################################              DATA    TYPES                       ##############################
#################################################                     ################################################

cdef class Schnitz:
    cdef np.ndarray data
    cdef np.ndarray time
    cdef np.ndarray volume
    cdef Schnitz parent
    cdef Schnitz daughter1
    cdef Schnitz daughter2

    cdef inline np.ndarray get_data(self):
        return self.data

    cdef inline np.ndarray get_time(self):
        return self.time

    cdef inline np.ndarray get_volume(self):
        return self.volume

    cdef inline Schnitz get_parent(self):
        return self.parent

    cdef inline Schnitz get_daughter_1(self):
        return self.daughter1

    cdef inline Schnitz get_daughter_2(self):
        return self.daughter2

    cdef inline void set_parent(self, Schnitz p):
        self.parent = p

    cdef inline void set_daughters(self,Schnitz d1, Schnitz d2):
        self.daughter1 = d1
        self.daughter2 = d2


cdef class Lineage:

    cdef list schnitzes
    cdef vector[void*] c_schnitzes


    cdef inline void add_schnitz(self, Schnitz s):
        """
        Add a schnitz to the lineage.
        :param s: (Schnitz) Add s to the lineage.
        :return: None
        """
        self.schnitzes.append(s)
        self.c_schnitzes.push_back(<void*> s)

    cdef inline unsigned size(self):
        """
        Get the total number of schnitzes in the lineage.
        :return: (int) size of lineage
        """
        return self.c_schnitzes.size()

    cdef inline Schnitz get_schnitz(self, unsigned index):
        """
        Get a specific schnitz from the lineage
        :param index: (unsigned) the Schnitz to retrieve 0 <= index < size()
        :return: (Schnitz) the requested Schnitz
        """
        return (<Schnitz> (self.c_schnitzes[index]))

cdef class ExperimentalLineage(Lineage):
    cdef dict species_dict

