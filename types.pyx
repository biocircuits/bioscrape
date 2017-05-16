# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
from bs4 import BeautifulSoup
cimport numpy as np
cimport random as cyrandom
from vector cimport vector
import re

from libc.math cimport log, sqrt, cos, round, exp

##################################################                ####################################################
######################################              PROPENSITY TYPES                    ##############################
#################################################                     ################################################


cdef class Propensity:
    def __init__(self):
        """
        Set the propensity type enum variable.
        """
        self.propensity_type = PropensityType.unset
    def py_get_propensity(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params):
        """
        Calculate propensity in pure python given a state and parameter vector.
        :param state: (np.ndarray) state vector of doubles
        :param params: (np.ndarray) parameter vector of doubles.
        :return: (double) computed propensity, should be non-negative
        """
        return self.get_propensity(<double*> state.data, <double*> params.data)

    # must be overriden by the daughter class
    cdef double get_propensity(self, double* state, double* params):
        """
        Compute the propensity given state and parameters (MUST be overridden, this returns -1.0)
        :param state: (double *) pointer to state vector
        :param params: (double *) pointer to parameter vector
        :return: (double) computed propensity, should be non-negative
        """
        return -1.0

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        """
        Compute the propensity given state and parameters and volume. (MUST be overridden)
        :param state: (double *) pointer to state vector
        :param params:(double *) pointer to parameter vector
        :param volume: (double) the cell volume
        :return: (double) computed propensity, should be non-negative
        """
        return -1.0

    def py_get_volume_propensity(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                                 double volume):
        return self.get_volume_propensity(<double*> state.data, <double*> params.data, volume)



    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        """
        Initializes the parameters and species to look at the right indices in the state
        :param dictionary: (dict:str--> str) the fields for the propensity 'k','s1' etc map to the actual parameter
                                             and species names
        :param species_indices: (dict:str-->int) map species names to entry in species vector
        :param parameter_indices: (dict:str-->int) map param names to entry in param vector
        :return: nothing
        """
        pass

    def get_species_and_parameters(self, dict fields):
        """
        get which fields are species and which are parameters
        :param dict(str-->str) dictionary containing the XML attributes for that propensity to process.
        :return: (list(string), list(string)) First entry is the names of species, second entry is the names of parameters
        """
        return (None,None)



cdef class ConstitutivePropensity(Propensity):
    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.constitutive

    cdef double get_propensity(self, double* state, double* params):
        return params[self.rate_index]

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        return params[self.rate_index] * volume

    def initialize(self, dict dictionary, species_indices, parameter_indices):
        for key,value in dictionary.items():
            if key == 'k':
                self.rate_index = parameter_indices[value]
            elif key == 'species':
                pass
            else:
                print('Warning! Useless field for constitutive reaction', key)
    def get_species_and_parameters(self, dict fields):
        return ([],[fields['k']])




cdef class UnimolecularPropensity(Propensity):
    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.unimolecular

    cdef double get_propensity(self, double* state, double* params):
        return params[self.rate_index] * state[self.species_index]

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        return params[self.rate_index] * state[self.species_index]


    def initialize(self, dict dictionary, species_indices, parameter_indices):
        for key,value in dictionary.items():
            if key == 'species':
                self.species_index = species_indices[value]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                print('Useless field for unimolecular reaction', key)

    def get_species_and_parameters(self, dict fields):
        return ([ fields['species'] ],[ fields['k'] ])



cdef class BimolecularPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.bimolecular

    cdef double get_propensity(self, double* state, double* params):
        return params[self.rate_index] * state[self.s1_index] * state[self.s2_index]

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        return params[self.rate_index] * state[self.s1_index] * state[self.s2_index] / volume


    def initialize(self, dict dictionary, species_indices, parameter_indices):
        for key,value in dictionary.items():
            if key == 'species':
                species_names = [x.strip() for x in value.split('*')]
                species_names = [x for x in species_names if x != '']
                assert(len(species_names) == 2)
                self.s1_index = species_indices[species_names[0]]
                self.s2_index = species_indices[species_names[1]]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                print('Useless field for bimolecular reaction', key)

    def get_species_and_parameters(self, dict fields):
        return ([ x.strip() for x in fields['species'].split('*') ],[ fields['k'] ])


cdef class PositiveHillPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.hill_positive

    cdef double get_propensity(self, double* state, double* params):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        return rate * (X / K) ** n / (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        return rate * (X / K) ** n / (1 + (X/K)**n)

    def initialize(self, dict dictionary, species_indices, parameter_indices):
        for key,value in dictionary.items():
            if key == 's1':
                self.s1_index = species_indices[value]
            elif key == 'K':
                self.K_index = parameter_indices[value]
            elif key == 'n':
                self.n_index = parameter_indices[value]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                print('Warning! Useless field for Hill propensity', key)

    def get_species_and_parameters(self, dict fields):
        return ([ fields['s1'] ],[ fields['K'],fields['n'],fields['k'] ])


cdef class PositiveProportionalHillPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.proportional_hill_positive

    cdef double get_propensity(self, double* state, double* params):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        cdef double d = state[self.d_index]
        return rate * d *  (X / K) ** n / (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double d = state[self.d_index]
        cdef double rate = params[self.rate_index]
        return d * rate * (X / K) ** n / (1 + (X/K)**n)


    def initialize(self, dict dictionary, species_indices, parameter_indices):
        for key,value in dictionary.items():
            if key == 's1':
                self.s1_index = species_indices[value]
            elif key == 'd':
                self.d_index = species_indices[value]
            elif key == 'K':
                self.K_index = parameter_indices[value]
            elif key == 'n':
                self.n_index = parameter_indices[value]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                print('Warning! Useless field for proportional Hill propensity', key)


    def get_species_and_parameters(self, dict fields):
        return ([ fields['s1'], fields['d'] ],[ fields['K'],fields['n'],fields['k'] ])



cdef class NegativeHillPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.hill_negative

    cdef double get_propensity(self, double* state, double* params):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        return rate * 1 / (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        return rate * 1 / (1 + (X/K)**n)

    def initialize(self, dict dictionary, species_indices, parameter_indices):
        for key,value in dictionary.items():
            if key == 's1':
                self.s1_index = species_indices[value]
            elif key == 'K':
                self.K_index = parameter_indices[value]
            elif key == 'n':
                self.n_index = parameter_indices[value]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                print('Warning! Useless field for Hill propensity', key)

    def get_species_and_parameters(self, dict fields):
        return ([ fields['s1'] ],[ fields['K'],fields['n'],fields['k'] ])



cdef class NegativeProportionalHillPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.proportional_hill_negative

    cdef double get_propensity(self, double* state, double* params):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        cdef double d = state[self.d_index]
        return rate * d *  1/ (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double d = state[self.d_index]
        cdef double rate = params[self.rate_index]
        return d * rate * 1 / (1 + (X/K)**n)


    def initialize(self, dict dictionary, species_indices, parameter_indices):
        for key,value in dictionary.items():
            if key == 's1':
                self.s1_index = species_indices[value]
            elif key == 'd':
                self.d_index = species_indices[value]
            elif key == 'K':
                self.K_index = parameter_indices[value]
            elif key == 'n':
                self.n_index = parameter_indices[value]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                print('Warning! Useless field for proportional Hill propensity', key)

    def get_species_and_parameters(self, dict fields):
        return ([ fields['s1'], fields['d'] ],[ fields['K'],fields['n'],fields['k'] ])

    def set_species(self, species, species_indices):
        for key in species:
            if key == 's1':
                self.s1_index = species_indices[species['s1']]
            elif key == 'd':
                self.d_index = species_indices[species['d']]
            else:
                print('Warning! Useless species for Hill propensity', key)
    def set_parameters(self,parameters, parameter_indices):
        for key in parameters:
            if key == 'K':
                self.K_index = parameter_indices[parameters[key]]
            elif key == 'n':
                self.n_index = parameter_indices[parameters[key]]
            elif key == 'k':
                self.rate_index = parameter_indices[parameters[key]]
            else:
                print('Warning! Useless parameter for Hill propensity', key)



cdef class MassActionPropensity(Propensity):
    def __init__(self):
        self.propensity_type = PropensityType.mass_action

    cdef double get_propensity(self, double* state, double* params):
        cdef double ans = params[self.k_index]
        cdef int i
        for i in range(self.num_species):
            ans *= state[self.sp_inds[i]]

        return ans

    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        cdef double ans = params[self.k_index]
        cdef int i
        for i in range(self.num_species):
            ans *= state[self.sp_inds[i]]
        if self.num_species == 1:
            return ans
        elif self.num_species == 2:
            return ans / volume
        elif self.num_species == 0:
            return ans * volume
        else:
            return ans / (volume ** (self.num_species - 1) )


    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        for key, value in dictionary.items():
            if key == 'species':
                species_names = [s.strip() for s in value.split('*')]
                for species_name in species_names:
                    if species_name == '':
                        continue
                    self.sp_inds.push_back(species_indices[species_name])
                self.num_species = self.sp_inds.size()
            elif key == 'k':
                self.k_index = parameter_indices[value]
            else:
                print('Warning: useless field for mass action propensity', key)

    def get_species_and_parameters(self, dict fields):
        species_list = [x.strip()   for x in fields['species'].split('*') ]
        species_list = [x for x in species_list if x != '']

        return (species_list, [ fields['k'] ])


##################################################                ####################################################
######################################              PARSING                             ##############################
#################################################                     ################################################


cdef class Term:
    cdef double evaluate(self, double *species, double *params):
        raise SyntaxError('Cannot make Term base object')

    cdef double volume_evaluate(self, double *species, double *params, double vol):
        raise SyntaxError('Cannot make Term base object')


    def py_evaluate(self, np.ndarray species, np.ndarray params):
        return self.evaluate(<double*> species.data, <double*> params.data)

    def py_volume_evaluate(self, np.ndarray species, np.ndarray params,
                           double vol):
        return self.volume_evaluate(<double*> species.data, <double*> params.data,
                                    vol)

# Base building blocks

cdef class ConstantTerm(Term):

    def __init__(self, double val):
        self.value = val

    cdef double evaluate(self, double *species, double *params):
        return self.value
    cdef double volume_evaluate(self, double *species, double *params, double vol):
        return self.value

cdef class SpeciesTerm(Term):


    def __init__(self, unsigned ind):
        self.index = ind

    cdef double evaluate(self, double *species, double *params):
        return species[self.index]
    cdef double volume_evaluate(self, double *species, double *params, double vol):
        return species[self.index]

cdef class ParameterTerm(Term):

    def __init__(self, unsigned ind):
        self.index = ind

    cdef double evaluate(self, double *species, double *params):
        return params[self.index]
    cdef double volume_evaluate(self, double *species, double *params, double vol):
        return params[self.index]

cdef class VolumeTerm(Term):
    cdef double evaluate(self, double *species, double *params):
        return 1.0
    cdef double volume_evaluate(self, double *species, double *params, double vol):
        return vol

# Putting stuff together

cdef class SumTerm(Term):


    def __init__(self):
        self.terms_list = []

    cdef void add_term(self,Term trm):
        self.terms.push_back(<void*> trm)
        self.terms_list.append(trm)

    cdef double evaluate(self, double *species, double *params):
        cdef double ans = 0.0
        cdef unsigned i
        for i in range(self.terms.size()):
            ans += (<Term>(self.terms[i])).evaluate(species, params)
        return ans

    cdef double volume_evaluate(self, double *species, double *params, double vol):
        cdef double ans = 0.0
        cdef unsigned i
        for i in range(self.terms.size()):
            ans += (<Term>(self.terms[i])).volume_evaluate(species,params,vol)
        return ans

cdef class ProductTerm(Term):


    def __init__(self):
        self.terms_list = []

    cdef void add_term(self,Term trm):
        self.terms.push_back(<void*> trm)
        self.terms_list.append(trm)

    cdef double evaluate(self, double *species, double *params):
        cdef double ans = 1.0
        cdef unsigned i
        for i in range(self.terms.size()):
            ans *= (<Term>(self.terms[i])).evaluate(species, params)
        return ans

    cdef double volume_evaluate(self, double *species, double *params, double vol):
        cdef double ans = 1.0
        cdef unsigned i
        for i in range(self.terms.size()):
            ans *= (<Term>(self.terms[i])).volume_evaluate(species,params,vol)
        return ans

cdef class PowerTerm(Term):


    cdef void set_base(self, Term base):
        self.base = base
    cdef void set_exponent(self, Term exponent):
        self.exponent = exponent

    cdef double evaluate(self, double *species, double *params):
        return self.base.evaluate(species,params) ** \
               self.exponent.evaluate(species,params)

    cdef double volume_evaluate(self, double *species, double *params, double vol):
        return self.base.volume_evaluate(species,params,vol) ** \
               self.exponent.volume_evaluate(species,params,vol)


# Parsing part is here.

def locations_outside_parentheses(c, s):
    paren_count = 0

    locs = []
    for i in range(len(s)):
        if s[i] == '(':
            paren_count += 1
        elif s[i] == ')':
            paren_count -= 1
        elif s[i] == c and paren_count == 0:
            locs.append(i)

        if paren_count < 0:
            raise SyntaxError('Strange parentheses: ', s)

    return locs

def is_single_term(s):
    if s[0] == '(' and s[len(s)-1] == ')':
        min_paren_count = 100
        paren_count = 0
        for i in range(len(s)-1):
            if s[i] == '(':
                paren_count += 1
            elif s[i] == ')':
                paren_count -= 1

            if paren_count < min_paren_count:
                min_paren_count = paren_count

        if min_paren_count > 0:
            return True
        else:
            return False
    return False



def is_a_number(s):
    try:
        a = float(s)
        return True
    except ValueError:
        return False

def parse_expression(instring, species2index, params2index):
    instring = instring.strip()

    print('Calling it with: %s' % (instring))

    # try parentheses first
    if is_single_term(instring):
        return parse_expression(instring[1:len(instring)-1],
                                species2index,
                                params2index)

    # try plus next
    l = locations_outside_parentheses('+', instring)
    cdef SumTerm sumterm
    cdef ProductTerm productterm
    cdef PowerTerm powerterm

    if len(l) > 0:
        sumterm = SumTerm()
        l.insert(0,-1)
        l.append(len(instring))
        for i in range(len(l)-1):
            sumterm.add_term(parse_expression(instring[ l[i]+1 : l[i+1] ],
                                              species2index, params2index))
        return sumterm

    # try product term next
    l = locations_outside_parentheses('*', instring)
    if len(l) > 0:
        productterm = ProductTerm()
        l.insert(0,-1)
        l.append(len(instring))
        for i in range(len(l)-1):
            productterm.add_term(parse_expression(instring[ l[i]+1 : l[i+1] ],
                                              species2index, params2index))
        return productterm

    # try exponential term next
    l = locations_outside_parentheses('^', instring)
    if len(l) > 0:
        assert (len(l) == 1)
        powerterm = PowerTerm()
        powerterm.set_base( parse_expression(instring[0:l[0]],
                                            species2index,
                                            params2index) )
        powerterm.set_exponent( parse_expression(instring[l[0]+1:],
                                            species2index,
                                            params2index) )
        return powerterm

    # otherwise must be a building block
    if is_a_number(instring):
        return ConstantTerm(float(instring))

    if instring == 'volume':
        return VolumeTerm()

    if instring[0] == '|':
        return ParameterTerm(params2index[ instring[1:] ])

    if instring in species2index:
        return SpeciesTerm(species2index[instring])

    raise SyntaxError('Something is horribly wrong with your expression!',
                       instring)


cdef class GeneralPropensity(Propensity):

    cdef double get_propensity(self, double* state, double* params):
        return self.term.evaluate(state,params)
    cdef double get_volume_propensity(self, double *state, double *params, double volume):
        return self.term.volume_evaluate(state,params,volume)

    def __init__(self):
        self.propensity_type = PropensityType.general

    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        instring = dictionary['rate']

        self.term = parse_expression(instring, species_indices, parameter_indices)


    def get_species_and_parameters(self, dict fields):
        instring = fields['rate'].strip()
        # split on all the operators and then strip down all whitespace
        names = re.split('[*^+()]',  instring)
        names = [x.strip() for x in names]
        names = [x for x in names if x != '']

        species_names = []
        param_names = []

        for name in names:
            if is_a_number(name) or name == 'volume':
                continue
            elif name[0] == '|':
                param_names.append(name[1:])
            else:
                species_names.append(name)

        return species_names, param_names






##################################################                ####################################################
######################################              DELAY TYPES                        ##############################
#################################################                     ################################################

cdef class Delay:
    def __init__(self):
        """
        Set the delay_type attribute to the appropriate enum value.
        """

        self.delay_type = DelayType.unset_delay

    def py_get_delay(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params):
        """
        Return the delay given the state and parameter vector
        :param state: (np.ndarray) the state vector
        :param params: (np.ndarray) the parameters vector
        :return: (double) the computed delay

        This function should NOT be overridden by subclases. It is just a Python wrapped of the cython delay function.
        """
        return self.get_delay(<double*> state.data, <double*> params.data)



    cdef double get_delay(self, double* state, double* params):
        """
        Compute a delay given the state and parameters vectors.
        :param state: (double *) the array containing the state vector
        :param params: (double *) the array containing the parameters vector
        :return: (double) the computed delay.

        This function must be overridden by subclasses.
        """

        return -1.0

    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        """
        Initializes the parameters and species to look at the right indices in the state
        :param dictionary: (dict:str--> str) the fields for the propensity 'k','s1' etc map to the actual parameter
                                             and species names
        :param species_indices: (dict:str-->int) map species names to entry in species vector
        :param parameter_indices: (dict:str-->int) map param names to entry in param vector
        :return: nothing
        """
        pass

    def get_species_and_parameters(self, dict fields):
        """
        get which fields are species and which are parameters
        :return: (list(string), list(string)) First entry is the fields that are species, second entry is the fields
                                              that are parameters
        """
        return [],[]


cdef class NoDelay(Delay):
    def __init__(self):
        self.delay_type = DelayType.none

    cdef double get_delay(self, double* state, double* params):
        return 0.0

cdef class FixedDelay(Delay):

    def __init__(self):
        self.delay_type = DelayType.fixed

    cdef double get_delay(self, double* state, double* params):
        return params[self.delay_index]

    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        for key,value in dictionary.items():
            if key == 'delay':
                self.delay_index = parameter_indices[value]
            else:
                print('Warning! Useless field for fixed delay', key)

    def get_species_and_parameters(self, dict fields):
        return [], [fields['delay']]

cdef class GaussianDelay(Delay):

    def __init__(self):
        self.delay_type = DelayType.gaussian

    cdef double get_delay(self, double* state, double* params):
        return cyrandom.normal_rv(params[self.mean_index],params[self.std_index])


    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        for key,value in dictionary.items():
            if key == 'mean':
                self.mean_index = parameter_indices[value]
            elif key == 'std':
                self.std_index = parameter_indices[value]
            else:
                print('Warning! Useless field for gaussian delay', key)

    def get_species_and_parameters(self, dict fields):
        return [],[fields['mean'], fields['std']]



cdef class GammaDelay(Delay):

    def __init__(self):
        self.delay_type = DelayType.gamma

    cdef double get_delay(self, double* state, double* params):
        return cyrandom.gamma_rv(params[self.k_index],params[self.theta_index])

    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        for key,value in dictionary.items():
            if key == 'k':
                self.k_index = parameter_indices[value]
            elif key == 'theta':
                self.theta_index = parameter_indices[value]
            else:
                print('Warning! Useless field for gamma delay', key)

    def get_species_and_parameters(self, dict fields):
        return [],[fields['k'], fields['theta']]

##################################################                ####################################################
######################################              RULE   TYPES                       ###############################
#################################################                     ################################################

cdef class Rule:
    """
    A class for doing rules that must be done either at the beginning of a simulation or repeatedly at each step of
    the simulation.
    """
    cdef void execute_rule(self, double *state, double *params):
        raise NotImplementedError('Creating base Rule class. This should be subclassed.')

    cdef void execute_volume_rule(self, double *state, double *params, double volume):
        self.execute_rule(state, params)

    def py_execute_rule(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params):
        self.execute_rule(<double*> state.data, <double*> params.data)

    def py_execute_volume_rule(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                               double volume ):
        self.execute_volume_rule(<double*> state.data, <double*> params.data, volume)

    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        """
        Initializes the parameters and species to look at the right indices in the state
        :param dictionary: (dict:str--> str) the fields for the propensity 'k','s1' etc map to the actual parameter
                                             and species names
        :param species_indices: (dict:str-->int) map species names to entry in species vector
        :param parameter_indices: (dict:str-->int) map param names to entry in param vector
        :return: nothing
        """
        pass

    def get_species_and_parameters(self, dict fields):
        """
        get which fields are species and which are parameters
        :param dict(str-->str) dictionary containing the XML attributes for that propensity to process.
        :return: (list(string), list(string)) First entry is the names of species, second entry is the names of parameters
        """
        return (None,None)


cdef class AdditiveAssignmentRule(Rule):
    """
    A class for assigning a species to a sum of a bunch of other species.
    """

    cdef void execute_rule(self, double *state, double *params):
        cdef unsigned i = 0
        cdef double answer = 0.0
        for i in range(self.species_source_indices.size()):
            answer += state[ self.species_source_indices[i] ]

        state[self.dest_index] = answer

    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices):
        equation = dictionary['equation']
        split_eqn = [s.strip() for s in equation.split('=') ]
        assert(len(split_eqn) == 2)
        dest_name = split_eqn[0]
        src_names = [s.strip() for s in split_eqn[1].split('+')]

        self.dest_index = species_indices[dest_name]

        for string in src_names:
            self.species_source_indices.push_back(  species_indices[string]  )

    def get_species_and_parameters(self, dict fields):
        # Add the species names
        equation = fields['equation']
        split_eqn = [s.strip() for s in equation.split('=') ]
        assert(len(split_eqn) == 2)
        dest_name = split_eqn[0]
        species_names = [s.strip() for s in split_eqn[1].split('+')]
        species_names.append(dest_name)
        return species_names, []



cdef class GeneralAssignmentRule(Rule):
    """
    A class for doing rules that must be done either at the beginning of a simulation or repeatedly at each step of
    the simulation.
    """
    cdef void execute_rule(self, double *state, double *params):
        state[self.dest_index] = self.rhs.evaluate(state,params)

    cdef void execute_volume_rule(self, double *state, double *params, double volume):
        state[self.dest_index] = self.rhs.volume_evaluate(state,params,volume)

    def initialize(self, dict fields, species2index, params2index):
        self.rhs = parse_expression(fields['equation'].split('=')[1], species2index, params2index)
        self.dest_index = species2index[ fields['equation'].split('=')[0].strip() ]

    def get_species_and_parameters(self, dict fields):
        instring = fields['equation'].strip()

        dest_name = instring.split('=')[0].strip()

        instring = instring.split('=')[1]

        # split on all the operators and then strip down all whitespace
        names = re.split('[*^+()]',  instring)
        names = [x.strip() for x in names]
        names = [x for x in names if x != '']

        species_names = [dest_name]
        param_names = []

        for name in names:
            if is_a_number(name) or name == 'volume':
                continue
            elif name[0] == '|':
                param_names.append(name[1:])
            else:
                species_names.append(name)

        return species_names, param_names


##################################################                ####################################################
######################################              VOLUME TYPES                        ##############################
#################################################                     ################################################

cdef class Volume:
    cdef double get_volume_step(self, double *state, double *params, double time, double volume, double dt):
        """
        Return the volume change in a time step of dt ending at time t given the state, parameters, and volume at t-d

        Must be overridden by subclass

        :param state: (double *) pointer to state vector
        :param params: (double *) pointer to parameter vector
        :param time: (double) ending time after the volume step has occurred
        :param volume: (double) the volume before the volume step occurs
        :param dt: (double) the time step over which you want the volume change
        :return: (double) the change in cell volume from time - dt to time
        """

        return 0.0


    def py_get_volume_step(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                           double time, double volume, double dt):
        return self.get_volume_step(<double*> state.data, <double*> params.data, time, volume, dt)


    cdef void initialize(self, double *state, double *params, double time, double volume):
        """
        Initialize the volume object given a new initial time and volume and the current state and parameters.

        This is required in order to handle non-memoryless properties, like the cell division time, which can be
        pre-sampled once in the initialize() function and then simply queried later.

        Must be overridden by subclass.

        :param state: (double *) pointer to the state vector
        :param params: (double *) pointer to the parameter vector
        :param time: (double) current initial time
        :param volume: (double) initial volume
        :return: None
        """

        pass

    def py_initialize(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                      double time, double volume):
        self.initialize(<double*> state.data, <double*> params.data, time, volume)

    cdef unsigned cell_divided(self, double *state, double *params, double time, double volume, double dt):
        """
        Return true or false if the cell divided during the time interval between time-dt and time. Note, in order
        to compute this the cell should have already been updated up to time t first.

        Must be overridden by subclass.

        :param state: (double *) pointer to the state vector
        :param params: (double *) pointer to parameter vector
        :param time: (double) the ending time of the time step
        :param volume: (double) the volume AFTER the time step occurred
        :param dt: (double) the width of the time step
        :return: 1 if the cell divided in [time-dt, time] or 0 if it did not divide.
        """

        return 0

    def py_cell_divided(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params, double time,
                        double volume, double dt):
        return self.cell_divided(<double*> state.data, <double*> params.data, time, volume, dt)


    def py_set_volume(self, double v):
        self.set_volume(v)
    def py_get_volume(self):
        return self.get_volume()




cdef class StochasticTimeThresholdVolume(Volume):
    def __init__(self, double cell_cycle_time, double average_division_volume, double division_noise):
        """
        Initialize the class with the cell cycle time, average division volume, and noise parameter.

        :param cell_cycle_time: (double) cell cycle time on average
        :param average_division_volume: (double) average volume at division
        :param division_noise: (double) noise in the cell cycle time as a relative c.o.v.
        """
        self.cell_cycle_time = cell_cycle_time
        self.average_division_volume = average_division_volume
        self.division_noise = division_noise
        self.division_time = -1.0

        # Compute growth rate yourself.
        self.growth_rate = 0.69314718056 / cell_cycle_time # log(2) / cycle time

    cdef double get_volume_step(self, double *state, double *params, double time, double volume, double dt):
        """
        Compute a deterministic volume step that is independent of state and parameters.

        :param state: (double *) the state vector. not used here
        :param params: (double *) the parameter vector. not used here
        :param time: (double) the final time.
        :param volume: (double) the volume at time - dt
        :param dt: (double) the time step
        :return:
        """

        return ( exp(self.growth_rate*dt) - 1.0) * volume

    cdef void initialize(self, double *state, double *params, double time, double volume):
        """
        Initialize the volume by setting initial time and volume and sampling the division time ahead of time with
        the the time left to division being the deterministic time left given the growth rate, cell cycle time,
        average division volume, and current volume. Then the actual time left to division is normal(1.0, noise) * T,
         where noise is the division noise parameter. This sets the future division time.

        :param state: (double *) the state vector. not used here
        :param params: (double *) the parameter vector. not used here
        :param time: (double) current time
        :param volume: (double) current volume
        :return:
        """

        self.set_volume(volume)
        cdef double time_left = log(self.average_division_volume / volume) / self.growth_rate
        time_left = cyrandom.normal_rv(1.0, self.division_noise) * time_left
        self.division_time = time + time_left
        #print("Volume:", volume, "Division Time:", self.division_time )

    cdef unsigned  cell_divided(self, double *state, double *params, double time, double volume, double dt):
        """
        Check if the cell has divided in the interval time-dt to time. Does not depend on any of the parameters for
         this volume type.

        :param state: (double *) the state vector. not used here
        :param params: (double *) the parameter vector. not used here
        :param time: (double) current time
        :param volume: (double) current volume
        :param dt: (double) time step
        :return: 1 if cell divided, 0 otherwise
        """


        if self.division_time > time - dt and self.division_time <= time:
            return 1
        return 0



##################################################                ####################################################
######################################              MODEL   TYPES                       ##############################
#################################################                     ################################################

cdef class Model:
    def __init__(self, filename):
        """
        Read in a model from a file using XML format for the model.

        :param filename: (str) the file to read the model
        """

        self._next_species_index = 0
        self._next_params_index = 0
        self.parse_model(filename)

    def _add_species(self, species):
        """
        Helper function for putting together the species vector (converting species names to indices in vector)

        If the species has already been added, then do nothing. otherwise give it a new index, and increase
        the next_species_index by 1

        :param species: (str) the species name
        :return: None
        """

        if species not in self.species2index:
            self.species2index[species] = self._next_species_index
            self._next_species_index += 1


    def _add_param(self, param):
        """
        Helper function for putting together the parameter vector (converting parameter names to indices in vector)

        If the parameter name has already been seen, then do nothing. Otherwise, give it a new index, and increase the
        next_params_index by 1.

        :param param: (str) the parameter name
        :return: None
        """

        if param not in self.params2index:
            self.params2index[param] = self._next_params_index
            self._next_params_index += 1


    def parse_model(self, filename):
        """
        Parse the model from the file filling in all the local variables (propensities, delays, update arrays). Also
        maps the species and parameters to indices in a species and parameters vector.

        :param filename: (str) the model file.
        :return: None
        """


        # open XML file from the filename and use BeautifulSoup to parse it
        xml_file = open(filename,'r')
        xml = BeautifulSoup(xml_file,features="xml")
        xml_file.close()
        # Go through the reactions and parse them 1 by 1 keeping track of species and reactions

        # Brief Outline
        #
        # Any time a species or parameter name is seen, add it to the index mapping names to indices if it has not
        # already been added.
        #
        # 1. For each reaction XML tag, parse the text to get the reactants and products. create a dictionary for each
        #    reaction that maps the species involved in each reaction to its update in the reaction i.e. for TX, you
        #    would have reaction_update_dict['mRNA'] = +1.0
        # 2. For each reaction, also do the same thing for the delayed updates.
        # 3. Parse the propensity and delay for each reaction and create the appropriate object for each and initialize
        #    by calling set_species and set_parameters for each.
        # 4. At the very end, with the params2index and species2index fully populated, use the saved updated dicts to re-
        #    construct the update array and delay update array.
        # 5. Read in the intial species and parameters values. If a species is not set, print a warning and set to 0.
        #    If a parameter is not set, throw an error.



        self._next_species_index = 0
        self._next_params_index = 0
        self.species2index = {}
        self.params2index = {}
        self.propensities = []
        self.delays = []
        self.repeat_rules = []

        reaction_updates = []
        delay_reaction_updates = []
        reaction_index = 0

        Reactions = xml.find_all('reaction')
        for reaction in Reactions:
            # create a new set of updates
            reaction_update_dict = {}

            # Parse the stoichiometry
            text = reaction['text']
            reactants = [s for s in [r.strip() for r in text.split('--')[0].split('+')] if s]
            products = [s for s in [r.strip() for r in text.split('--')[1].split('+')] if s]

            for r in reactants:
                # if the species hasn't been seen add it to the index
                self._add_species(r)
                # update the update array
                if r not in reaction_update_dict:
                    reaction_update_dict[r] = 0
                reaction_update_dict[r]  -= 1

            for p in products:
                # if the species hasn't been seen add it to the index
                self._add_species(p)
                # update the update array
                if p not in reaction_update_dict:
                    reaction_update_dict[p] = 0
                reaction_update_dict[p]  += 1

            reaction_updates.append(reaction_update_dict)


            # parse the delayed part of the reaction the same way as we did before.
            delay_reaction_update_dict = {}

            if reaction.has_attr('after'):
                text = reaction['after']
                reactants = [s for s in [r.strip() for r in text.split('--')[0].split('+')] if s]
                products = [s for s in [r.strip() for r in text.split('--')[1].split('+')] if s]

                for r in reactants:
                    # if the species hasn't been seen add it to the index
                    self._add_species(r)
                    # update the update array
                    if r not in delay_reaction_update_dict:
                        delay_reaction_update_dict[r] = 0
                    delay_reaction_update_dict[r]  -= 1

                for p in products:
                    # if the species hasn't been seen add it to the index
                    self._add_species(p)
                    # update the update array
                    if p not in delay_reaction_update_dict:
                        delay_reaction_update_dict[p] = 0
                    delay_reaction_update_dict[p]  += 1

            delay_reaction_updates.append(delay_reaction_update_dict)


            # Then look at the propensity and set up a propensity object
            propensity = reaction.find_all('propensity')
            if len(propensity) != 1:
                raise SyntaxError('Incorrect propensity tags in XML model\n' + propensity)
            propensity = propensity[0]
            # go through propensity types

            init_dictionary = propensity.attrs
            
            if propensity['type'] == 'hillpositive':
                prop_object = PositiveHillPropensity()

            elif propensity['type'] == 'proportionalhillpositive':
                prop_object = PositiveProportionalHillPropensity()

            elif propensity['type'] == 'hillnegative':
                prop_object = NegativeHillPropensity()

            elif propensity['type'] == 'proportionalhillnegative':
                prop_object = NegativeProportionalHillPropensity()

            elif propensity['type'] == 'massaction':
                species_names = [s.strip() for s in propensity['species'].split('*') ]
                species_names = [x for x in species_names if x != '']

                # if mass action propensity has less than 3 things, then use consitutitve, uni, bimolecular for speed.
                if len(species_names) == 0:
                    prop_object = ConstitutivePropensity()
                elif len(species_names) == 1:
                    prop_object = UnimolecularPropensity()
                elif len(species_names) == 2:
                    prop_object = BimolecularPropensity()
                else:
                    prop_object = MassActionPropensity()

            elif propensity['type'] == 'general':
                prop_object = GeneralPropensity()

            else:
                raise SyntaxError('Propensity Type makes no sense: ' + propensity['type'])

            species_names, param_names = prop_object.get_species_and_parameters(init_dictionary)

            for species_name in species_names:
                self._add_species(species_name)
            for param_name in param_names:
                self._add_param(param_name)
            
            init_dictionary.pop('type')
            prop_object.initialize(init_dictionary,self.species2index,self.params2index)

            self.propensities.append(prop_object)
            self.c_propensities.push_back(<void*> prop_object)


            # Then look at the delay and set up a delay object
            delay = reaction.find_all('delay')
            if len(delay) != 1:
                raise SyntaxError('Incorrect delay spec')
            delay = delay[0]
            init_dictionary = delay.attrs

            if delay['type'] == 'none':
                delay_object = NoDelay()

            elif delay['type'] == 'fixed':
                delay_object = FixedDelay()

            elif delay['type'] == 'gaussian':
                delay_object = GaussianDelay()

            elif delay['type'] == 'gamma':
                delay_object = GammaDelay()

            else:
                raise SyntaxError('Unknown delay type: ' + delay['type'])

            species_names, param_names = delay_object.get_species_and_parameters(init_dictionary)

            for species_name in species_names:
                self._add_species(species_name)
            for param_name in param_names:
                self._add_param(param_name)

            init_dictionary.pop('type',None)
            delay_object.initialize(init_dictionary,self.species2index,self.params2index)

            self.delays.append(delay_object)
            self.c_delays.push_back(<void*> delay_object)


        # Parse through the rules

        Rules = xml.find_all('rule')
        cdef Rule rule_object
        for rule in Rules:
            init_dictionary = rule.attrs
            # Parse the rule by rule type
            if rule['type'] == 'additive':
                rule_object = AdditiveAssignmentRule()
            elif rule['type'] == 'assignment':
                rule_object = GeneralAssignmentRule()
            else:
                raise SyntaxError('Invalid type of Rule: ' + rule['type'])

            # Add species and params to model
            species_names, params_names = rule_object.get_species_and_parameters(init_dictionary)
            for s in species_names: self._add_species(s)
            for p in params_names: self._add_param(p)

            # initialize the rule
            init_dictionary.pop('type')
            rule_object.initialize(init_dictionary,self.species2index,self.params2index)
            # Add the rule to the right place
            if rule['frequency'] == 'repeated':
                self.repeat_rules.append(rule_object)
                self.c_repeat_rules.push_back(<void*> rule_object)
            else:
                raise SyntaxError('Invalid Rule Frequency: ' + rule['frequency'])

        # With all reactions read in, generate the update array

        num_species = len(self.species2index.keys())
        num_reactions = len(Reactions)
        self.update_array = np.zeros((num_species, num_reactions))
        self.delay_update_array = np.zeros((num_species,num_reactions))
        for reaction_index in range(num_reactions):
            reaction_update_dict = reaction_updates[reaction_index]
            delay_reaction_update_dict = delay_reaction_updates[reaction_index]
            for sp in reaction_update_dict:
                self.update_array[self.species2index[sp],reaction_index] = reaction_update_dict[sp]
            for sp in delay_reaction_update_dict:
                self.delay_update_array[self.species2index[sp],reaction_index] = delay_reaction_update_dict[sp]


        # Generate species values and parameter values
        self.params_values = np.empty(len(self.params2index.keys()), )
        self.params_values.fill(np.nan)
        Parameters = xml.find_all('parameter')
        for param in Parameters:
            param_value = float(param['value'])
            param_name = param['name']
            if param_name not in self.params2index:
                print ('Warning! Useless parameter '+ param_name)
            else:
                param_index = self.params2index[param_name]
                self.params_values[param_index] = param_value

        for p in self.params_values:
            if np.isnan(p):
                raise SyntaxError('Did not specify all parameter values!')

        self.species_values = np.empty(len(self.species2index.keys()), )
        self.species_values.fill(np.nan)
        Species = xml.find_all('species')
        for species in Species:
            species_value = float(species['value'])
            species_name = species['name']
            if species_name not in self.species2index:
                print ('Warning! Useless species value ' + species_name)
            else:
                species_index = self.species2index[species_name]
                self.species_values[species_index] = species_value

        for index in range(len(self.species_values)):
            if np.isnan(self.species_values[index]):
                print('Warning! Did not specify species value, setting it to zero!!! index =', index)

        self.species_values[np.isnan(self.species_values)] = 0.0


        print(self.species2index)
        print(self.params2index)
        print(self.update_array)
        print(self.delay_update_array)

    def get_params(self):
        """
        Get the set of parameter names.
        :return: (dict_keys str) the parameter names
        """

        return self.params2index.keys()

    def get_species(self):
        """
        Get the set of species names.
        :return: (dict_keys str) the species names
        """

        return self.species2index.keys()


    def set_params(self, param_dict):
        """
        Set parameter values

        :param param_dict: (dict:str -> double) Dictionary containing the parameters to set mapped to desired values.
        :return: None
        """

        param_names = set(self.params2index.keys())
        for p in param_dict:
            if p in param_names:
                self.params_values[self.params2index[p]] = param_dict[p]
            else:
                print('Useless parameter', p)


    def set_species(self, species_dict):
        """
        Set initial species values

        :param species_dict: (dict:str -> double) Dictionary containing the species to set mapped to desired values.
        :return: None
        """
        species_names = set(self.species2index.keys())
        for s in species_dict:
            if s in species_names:
                self.species_values[self.species2index[s]] = species_dict[s]
            else:
                print('Useless species level', s)

    cdef (vector[void*])* get_c_repeat_rules(self):
        """
        Get the set of rules to implement as a set of void pointers. Must be cast back to a Rule object to be used.
        This is much faster than accessing the Rules as a list though.
        :return: (vector[void*])* pointer to the vector of Rule objects
        """
        return & self.c_repeat_rules

    def get_propensities(self):
        """
        Get the propensities list.

        :return: (list) List of the propensities for each reaction.
        """
        return self.propensities

    def get_delays(self):
        """
        Get the delays list

        :return: (list) List of the delay objects for each reaction.
        """
        return self.delays

    cdef np.ndarray get_species_values(self):
        """
        Get the species values as an array
        :return: (np.ndarray) the species values
        """

        return self.species_values

    cdef np.ndarray get_params_values(self):
        """
        Get the parameter values as an array
        :return: (np.ndarray) the parameter values
        """
        return self.params_values

    cdef (vector[void*])* get_c_propensities(self):
        """
        Get the propensity objects for each reaction as a vector of void pointers. Must be cast back to Propensity
        type to use. This is much faster than accessing the list of propensities.
        :return: (vector[void*]*) Pointer to a vector of void pointers, where i-th void pointer points to propensity i
        """
        return & self.c_propensities

    cdef (vector[void*])* get_c_delays(self):
        """
        Get the delay objects for each reaction as a vector of void pointers. Must be cast back to Delay.
        :return: (vector[void*] *) Pointer to vector of void *, where the i-th void pointer points to delay for rxn i
        """
        return & self.c_delays

    cdef np.ndarray get_update_array(self):
        """
        Get the stoichiometric matrix for changes that occur immdeiately.
        :return: (np.ndarray) A 2-D array with 1 row per species, 1 column for each reaction.
        """
        return self.update_array

    def py_get_update_array(self):
        return self.update_array

    cdef np.ndarray get_delay_update_array(self):
        """
        Get the stoichiometric matrix for changes that occur after a delay.
        :return: (np.ndarray) A 2-D array with 1 row per species, 1 column for each reaction.
        """
        return self.delay_update_array

    def py_get_delay_update_array(self):
        return self.delay_update_array


    def get_param_index(self, param_name):
        if param_name in self.params2index:
            return self.params2index[param_name]
        return -1


    def get_species_index(self, species_name):
        if species_name in self.species2index:
            return self.species2index[species_name]
        return -1

    def get_param_value(self, param_name):
        if param_name in self.params2index:
            return self.params_values[self.params2index[param_name]]
        else:
            raise LookupError('No parameter with name '+ param_name)

    def get_species_value(self, species_name):
        if species_name in self.species2index:
            return self.species_values[self.species2index[species_name]]
        else:
            raise LookupError('No species with name '+ species_name)


##################################################                ####################################################
######################################              DATA    TYPES                       ##############################
#################################################                     ################################################

cdef class Schnitz:
    def __init__(self, time, data, volume):
        """
        Create a Schnitz with the provided time, data, and volume arrays. Parents and daughters are left as None and
        must be set later if required.

        :param time: (np.ndarray) 1-D array with time points
        :param data: (np.ndarray) 2-D array with one row for each time point, one column for each measured output
        :param volume: (np.ndarray) 1-D array with volume at each time point
        """
        self.parent = None
        self.daughter1 = None
        self.daughter2 = None
        self.time = time
        self.volume = volume
        self.data = data



    def py_get_data(self):
        return self.data

    def py_get_time(self):
        return self.time

    def py_get_volume(self):
        return self.volume

    def py_get_parent(self):
        return self.parent

    def py_get_daughters(self):
        return (self.daughter1, self.daughter2)



cdef class Lineage:
    def __init__(self):
        """
        Creates a lineage object.
        """
        self.schnitzes = []

    def py_size(self):
        """
        Get the total number of schnitzes in the lineage.
        :return: (int) size of lineage
        """
        return self.c_schnitzes.size()

    def py_get_schnitz(self, unsigned index):
        """
        Get a specific schnitz from the lineage
        :param index: (unsigned) the Schnitz to retrieve 0 <= index < size()
        :return: (Schnitz) the requested Schnitz
        """
        return (<Schnitz> (self.c_schnitzes[index]))





