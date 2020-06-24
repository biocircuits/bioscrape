# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

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
    """
    Class for defining a propensity. Must contain a propensity type as well as two functions.
    get_propensity returns a propensity given a state vector, parameter vector
    get_volume_propensity is the same but also requires a volume
    """
    cdef PropensityType propensity_type

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)


cdef class ConstitutivePropensity(Propensity):
    """
    A propensity class for constitutive propensities. (k)

    Attributes:
        rate_index (unsigned): the parameter index containing the rate k.
    """
    # variables
    cdef unsigned rate_index
    # constructor

    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class UnimolecularPropensity(Propensity):
    """
    A propensity class for a unimolecular propensity (k*x)

    Attributes:
        rate_index (unsigned): the parameter index containing the rate k
        species_index (unsigned): the species index containing the species x
    """


    # variables
    cdef unsigned rate_index
    cdef unsigned species_index


    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)


cdef class BimolecularPropensity(Propensity):
    """
    A propensity class for a bimolecular propensity (k*x1*x2)

    Attributes:
        rate_index (unsigned): the parameter index containing the rate k
        s1_index (unsigned): the species index containing the species x1
        s2_index (unsigned): the species index containing the species x2
    """

    # variables
    cdef unsigned rate_index
    cdef unsigned s1_index
    cdef unsigned s2_index


    cdef double get_propensity(self, double* state, double* params, double time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time)
    cdef double get_stochastic_propensity(self, double* state, double* params, double time)
    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time)

cdef class PositiveHillPropensity(Propensity):
    """
    A propensity class for an activating Hill function (   rate * (x1/K)**n / (1 + (x1/K)**n )   )

    Attributes:
        K_index (unsigned): parameter index of Hill constant K
        rate_index (unsigned): parameter index for rate value
        n_index (unsigned): parameter index for cooperativity value
        s1_index (unsigned): species index for x1
    """
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
    """
    A propensity class for an activating Hill function proportional to another second species
     (   rate * d * (x1/K)**n / (1 + (x1/K)**n )   )

    Attributes:
        K_index (unsigned): parameter index of Hill constant K
        rate_index (unsigned): parameter index for rate value
        n_index (unsigned): parameter index for cooperativity value
        s1_index (unsigned): species index for x1
        d_index (unsigned): species index for d
    """


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
    """
    A propensity class for a repressing Hill function (   rate * 1 / (1 + (x1/K)**n )   )

    Attributes:
        K_index (unsigned): parameter index of Hill constant K
        rate_index (unsigned): parameter index for rate value
        n_index (unsigned): parameter index for cooperativity value
        s1_index (unsigned): species index for x1
    """

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
    """
    A propensity class for a repressing Hill function proportional to another second species
     (   rate * d * 1 / (1 + (x1/K)**n )   )

    Attributes:
        K_index (unsigned): parameter index of Hill constant K
        rate_index (unsigned): parameter index for rate value
        n_index (unsigned): parameter index for cooperativity value
        s1_index (unsigned): species index for x1
        d_index (unsigned): species index for d
    """

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

cdef class SumTerm(Term):
    cdef vector[void*] terms
    cdef list terms_list

    cdef void add_term(self,Term trm)

    cdef double evaluate(self, double *species, double *params, double time)
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class ProductTerm(Term):
    cdef vector[void*] terms
    cdef list terms_list

    cdef void add_term(self,Term trm)
    cdef double evaluate(self, double *species, double *params, double time)

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class MaxTerm(Term):
    cdef vector[void*] terms
    cdef list terms_list

    cdef void add_term(self,Term trm)
    cdef double evaluate(self, double *species, double *params, double time)

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time)

cdef class MinTerm(Term):
    cdef vector[void*] terms
    cdef list terms_list

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
    """
    A delay class for different delay types. Must be overridden by subclasses.

    Contains one function get_delay that returns a delay as a function of the state and parameters.

    Attributes:
        delay_type (DelayType): the type of delay
        delay_params (dictionary): parameters for the delay distribution

    """
    # delay type
    cdef DelayType delay_type

    # to be overriden
    cdef double get_delay(self, double* state, double* params)

cdef class NoDelay(Delay):
    """
    A delay class for representing no delay. Will always return a delay of 0.0
    """

    cdef double get_delay(self, double* state, double* params)

cdef class FixedDelay(Delay):
    """
    A delay class for representing a fixed delay.

    Attributes:
        delay_index (unsigned): parameter index containing the fixed delay value.
    """

    cdef unsigned delay_index

    cdef double get_delay(self, double* state, double* params)


cdef class GaussianDelay(Delay):
    """
    A delay class for representing a Gaussian distributed delay with specified mean and std.

    Attributes:
        mean_index (unsigned): parameter index containing the mean delay
        std_index (unsigned): parameter index containing the standard deviation of the delay.
    """

    cdef unsigned mean_index
    cdef unsigned std_index


    cdef double get_delay(self, double* state, double* params)

cdef class GammaDelay(Delay):
    """
    A delay class for representing a gamma distributed delay with specified k and theta, where
    it is *like* k independent exponential variables with mean theta

    Mean is k*theta
    Variance is k*theta^2

    Attributes:
        k_index (unsigned): parameter index containing the k value
        theta_index (unsigned): parameter index containing the theta value
    """
    cdef unsigned k_index
    cdef unsigned theta_index

    cdef double get_delay(self, double* state, double* params)

##################################################                ####################################################
######################################              RULE   TYPES                       ###############################
#################################################                     ################################################

cdef class Rule:
    """
    A class for doing rules that must be done either at the beginning of a simulation or repeatedly at each step of
    the simulation.
    """
    cdef void execute_rule(self, double *state, double *params, double time)
    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time)


cdef class AdditiveAssignmentRule(Rule):
    """
    A class for assigning a species to a sum of a bunch of other species.
    """
    cdef vector[int] species_source_indices
    cdef unsigned dest_index

    cdef void execute_rule(self, double *state, double *params, double time)
    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time)

cdef class GeneralAssignmentRule(Rule):
    """
    A class for assigning a species to a product of a bunch of other species.
    """
    cdef Term rhs
    cdef unsigned dest_index
    cdef int param_flag # 1 if the assigned thing is a parameter, 0 if it's a species

    cdef void execute_rule(self, double *state, double *params, double time)
    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time)


cdef class GrowthRule(Rule):
    """
    A class for assigning rules to govern volume growth in cell division simulations
    """
    cdef void execute_rule(self, double *state, double *params, double time)
    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time)


cdef class DivisionRule(Rule):
    """
    A class for assigning rules to govern cell division
    """
    cdef void execute_rule(self, double *state, double *params, double time)
    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time)


cdef class DeathRule(Rule):
    """
    A class for assigning rules to govern cell death
    """
    cdef void execute_rule(self, double *state, double *params, double time)
    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time)



##################################################                ####################################################
######################################              VOLUME  TYPES                       ##############################
#################################################                     ################################################

cdef class Volume:
    """
    A volume class for different volume update types. Must be overridden by subclasses implementing volume models.


    Attributes:
        current_volume (double): the current volume of the cell. (should be positive)

    Subclasses must implement get_volume_step, initialize, cell_divided
    """

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
    """
    A volume class for a cell that grows at a deterministic exponential rate and divides after a random normally
    distributed time.


    The relevant parameters for this model are the average division volume of the cell, the cell growth rate, and the
    specified noise in the division time.

    Attributes:
        division_time (double): the time at which the cell will divide.
        average_division_volume (double): the average volume at which to divide.
        division_noise (double):  << 1, the noise in the division time.
        cell_cycle_time (double): the average cell cycle time for these cells.
        growth_rate (double): the volume growth rate g (dV/dt = g*V)
    """
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
    """
    A volume class for a cell where growth rate depends on state and the division volume is chosen randomly
    ahead of time with some randomness.

    Attributes:
        division_volume (double): the volume at which the cell will divide.
        average_division_volume (double): the average volume at which to divide.
        division_noise (double):  << 1, the noise in the division time (c.o.v.)
        growth_rate (Term): the growth rate evaluated based on the state
    """
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
    """
    A class for keeping track of a chemical reaction model with delays. Does not support volumes, those are supplied
    externally.

    Attributes:
        _next_species_index (unsigned): used internally for indexing
        _next_params_index (unsigned): used internally for indexing
        c_propensities (void *): vector of pointers to the model propensity objects. Must cast back to Propensity type
        c_delays (void *): vector of pointers to model delay objects. Must cast back to Delay type
        propensities (list): List of propensity objects. same contents as c_propensities but slower access
        delays (list): List of delay objects. Same contents as c_delays but slower access
        species2index (dict:str -> int): maps species names to index in species state vector.
        params2index (dict:str -> int): maps parameter names to index in parameter vector.
        species_values (np.ndarray): array of the initial species values
        params_values (np.ndarray): array of the parameter values
        update_array (np.ndarray): 2-D array containing stoichiometric matrix num_species x num_reactions
        delay_update_array (np.ndarray): 2-D array containing delay stoich matrix (i.e. updates that occur after delay)


    """
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
    """
    A class for representing the data acquired from a single cell trajectory.

    Attributes:
        data (np.ndarray): A 2-D array containing one column for each measured output, 1 row for each timepoint
        time (np.ndarray): A 1-D array containing the time points
        volume (np.ndarray): A 1-D array with the volume at each time point
        parent (Schnitz): The Schnitz for the parent of this cell (if it exists, None if not)
        daughter1, daughter2 (Schnitz): The Schnitz's for the daughters of this cell (if they exist, None if not)
    """
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
    """
    A class for keeping track of cell lineages consisting of many Schnitz's

    Attributes:
        schnitzes (list): A list of all the Schnitz's in the lineage.
        c_schnitzes (vector[void*]): A vector containing void pointers to all the Schnitz's in the lineage
    """

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

