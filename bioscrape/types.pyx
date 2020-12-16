# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
from bs4 import BeautifulSoup
cimport numpy as np
cimport random as cyrandom
from vector cimport vector
import re
import sympy
from sympy.abc import _clash1
import warnings
import logging
import libsbml
from bioscrape.sbmlutil import add_species, add_parameter, add_reaction, add_rule, create_sbml_model, import_sbml

from libc.math cimport log, sqrt, cos, round, exp, fabs

##################################################                ####################################################
######################################              PROPENSITY TYPES                    ##############################
#################################################                     ################################################


cdef class Propensity:
    def __init__(self):
        """
        Set the propensity type enum variable.
        """
        self.propensity_type = PropensityType.unset

    def py_get_propensity(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                          double time = 0.0):
        """
        Calculate propensity in pure python given a state and parameter vector.
        :param state: (np.ndarray) state vector of doubles
        :param params: (np.ndarray) parameter vector of doubles.
        :return: (double) computed propensity, should be non-negative
        """
        return self.get_propensity(<double*> state.data, <double*> params.data, time)

    # must be overriden by the daughter class
    cdef double get_propensity(self, double* state, double* params, double time):
        """
        Compute the propensity given state and parameters (MUST be overridden, this returns -1.0)
        :param state: (double *) pointer to state vector
        :param params: (double *) pointer to parameter vector
        :param time: (double) the current time
        :return: (double) computed propensity, should be non-negative
        """
        return -1.0

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        """
        Compute the propensity given state and parameters and volume. (MUST be overridden)
        :param state: (double *) pointer to state vector
        :param params:(double *) pointer to parameter vector
        :param volume: (double) the cell volume
        :param time: (double) the current time
        :return: (double) computed propensity, should be non-negative
        """
        return -1.0


    cdef double get_stochastic_propensity(self, double* state, double* params, double time):
        """
        By default, stochastic propensities are the same as deterministic propensities but can be overwritten for specific propensity types.
        """
        return self.get_propensity(state, params, time)

    cdef double get_stochastic_volume_propensity(self, double* state, double* params, double volume, double time):
        """
        By default, stochastic propensities are the same as deterministic propensities but can be overwritten for specific propensity types.
        """
        return self.get_volume_propensity(state, params, volume, time)

    def py_get_volume_propensity(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                                 double volume, double time = 0.0):
        return self.get_volume_propensity(<double*> state.data, <double*> params.data, volume, time)



    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
        """
        Initializes the parameters and species to look at the right indices in the state
        :param dictionary: (dict:str--> str) the fields for the propensity 'k','s1' etc map to the actual parameter
                                             and species names
        :param species_indices: (dict:str-->int) map species names to entry in species vector
        :param parameter_indices: (dict:str-->int) map param names to entry in param vector
        :return: nothing
        """
        pass

    def get_species_and_parameters(self, dict fields, **keywords):
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

    cdef double get_propensity(self, double* state, double* params, double time):
        return params[self.rate_index]

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        return params[self.rate_index] * volume

    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
        for key,value in param_dictionary.items():
            if key == 'k':
                self.rate_index = parameter_indices[value]
            elif key == 'species':
                pass
            elif key == "propensity_type":
                pass
            else:
                logging.info('Warning! Useless field for ConstitutivePropensity'+str(key))

    def get_species_and_parameters(self, dict fields, **keywords):
        return ([],[fields['k']])


cdef class UnimolecularPropensity(Propensity):
    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.unimolecular

    cdef double get_propensity(self, double* state, double* params, double time):
        return params[self.rate_index] * state[self.species_index]

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        return params[self.rate_index] * state[self.species_index]


    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):

        for key,value in param_dictionary.items():
            if key == 'species':
                self.species_index = species_indices[value]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            elif key == "propensity_type":
                pass
            else:
                logging.info('Warning! Useless field for UnimolecularPropensity '+str(key))

    def get_species_and_parameters(self, dict fields, **keywords):
        return ([ fields['species'] ],[ fields['k'] ])



cdef class BimolecularPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.bimolecular

    cdef double get_propensity(self, double* state, double* params, double time):
        return params[self.rate_index] * state[self.s1_index] * state[self.s2_index]

    cdef double get_stochastic_propensity(self, double* state, double* params, double time):
        if self.s1_index != self.s2_index:
            return params[self.rate_index] * state[self.s1_index] * state[self.s2_index]
        else:
            return params[self.rate_index]*state[self.s1_index]*max(state[self.s1_index]-1, 0)


    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        return params[self.rate_index] * state[self.s1_index] * state[self.s2_index] / volume

    cdef double get_stochastic_volume_propensity(self, double* state, double* params, double volume, double time):
        if self.s1_index != self.s2_index:
            return params[self.rate_index] * state[self.s1_index] * state[self.s2_index] / volume
        else:
            return params[self.rate_index]*state[self.s1_index]*max(state[self.s1_index]-1, 0) / volume


    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):

        for key,value in param_dictionary.items():
            if key == 'species':
                species_names = [x.strip() for x in value.split('*')]
                species_names = [x for x in species_names if x != '']
                assert(len(species_names) == 2)
                self.s1_index = species_indices[species_names[0]]
                self.s2_index = species_indices[species_names[1]]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                logging.info('Warning! Useless field for BimolecularPropensity'+str(key))

    def get_species_and_parameters(self, dict fields, **keywords):
        return ([ x.strip() for x in fields['species'].split('*') ],[ fields['k'] ])


cdef class PositiveHillPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.hill_positive

    cdef double get_propensity(self, double* state, double* params, double time):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        return rate * (X / K) ** n / (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        return rate * (X / K) ** n / (1 + (X/K)**n)

    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):

        for key,value in param_dictionary.items():
            if key == 's1':
                self.s1_index = species_indices[value]
            elif key == 'K':
                self.K_index = parameter_indices[value]
            elif key == 'n':
                self.n_index = parameter_indices[value]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                logging.info('Warning! Useless field for PositiveHillPropensity '+str(key))

    def get_species_and_parameters(self, dict fields, **keywords):
        return ([ fields['s1'] ],[ fields['K'],fields['n'],fields['k'] ])


cdef class PositiveProportionalHillPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.proportional_hill_positive

    cdef double get_propensity(self, double* state, double* params, double time):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        cdef double d = state[self.d_index]
        return rate * d *  (X / K) ** n / (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double d = state[self.d_index]
        cdef double rate = params[self.rate_index]
        return d * rate * (X / K) ** n / (1 + (X/K)**n)


    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):

        for key,value in param_dictionary.items():
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
                logging.info('Warning! Useless field for PositiveProportionalHillPropensity '+str(key))


    def get_species_and_parameters(self, dict fields, **keywords):
        return ([ fields['s1'], fields['d'] ],[ fields['K'],fields['n'],fields['k'] ])



cdef class NegativeHillPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.hill_negative

    cdef double get_propensity(self, double* state, double* params, double time):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        return rate * 1 / (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        return rate * 1 / (1 + (X/K)**n)

    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):

        for key,value in param_dictionary.items():
            if key == 's1':
                self.s1_index = species_indices[value]
            elif key == 'K':
                self.K_index = parameter_indices[value]
            elif key == 'n':
                self.n_index = parameter_indices[value]
            elif key == 'k':
                self.rate_index = parameter_indices[value]
            else:
                logging.info('Warning! Useless field for NegativeHillPropensity '+str(key))

    def get_species_and_parameters(self, dict fields, **keywords):
        return ([ fields['s1'] ],[ fields['K'],fields['n'],fields['k'] ])



cdef class NegativeProportionalHillPropensity(Propensity):

    # constructor
    def __init__(self):
        self.propensity_type = PropensityType.proportional_hill_negative

    cdef double get_propensity(self, double* state, double* params, double time):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        cdef double d = state[self.d_index]
        return rate * d *  1/ (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double d = state[self.d_index]
        cdef double rate = params[self.rate_index]
        return d * rate * 1 / (1 + (X/K)**n)


    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):

        for key,value in param_dictionary.items():
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
                logging.info('Warning! Useless field for NegativeProportionalHillPropensity '+str(key))

    def get_species_and_parameters(self, dict fields, **keywords):
        return ([ fields['s1'], fields['d'] ],[ fields['K'],fields['n'],fields['k'] ])

    def set_species(self, species, species_indices):
        for key in species:
            if key == 's1':
                self.s1_index = species_indices[species['s1']]
            elif key == 'd':
                self.d_index = species_indices[species['d']]
            else:
                logging.info('Warning! Useless field for NegativeProportionalHillPropensity '+str(key))
    def set_parameters(self, parameters, parameter_indices):
        for key in parameters:
            if key == 'K':
                self.K_index = parameter_indices[parameters[key]]
            elif key == 'n':
                self.n_index = parameter_indices[parameters[key]]
            elif key == 'k':
                self.rate_index = parameter_indices[parameters[key]]
            else:
                logging.info('Warning! Useless field for NegativeProportionalHillPropensity '+str(key))



cdef class MassActionPropensity(Propensity):
    def __init__(self):
        self.propensity_type = PropensityType.mass_action

    cdef double get_propensity(self, double* state, double* params, double time):
        cdef double ans = params[self.k_index]
        cdef int i
        for i in range(len(self.sp_inds)):
            ans *= state[self.sp_inds[i]]

        return ans

    cdef double get_stochastic_propensity(self, double* state, double* params, double time):
        cdef double ans = params[self.k_index]
        cdef int i
        for i in range(len(self.sp_inds)):
            for j in range(self.sp_counts[i]):
                ans *= max(state[self.sp_inds[i]]-j, 0)

        return ans

    cdef double get_volume_propensity(self, double *state, double *params,
                                      double volume, double time):
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

    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time):

        cdef double ans = self.get_stochastic_propensity(state, params, time)
        if self.num_species == 0:
            return ans*volume
        elif self.num_species == 1:
            return ans
        elif self.num_species == 2:
            return ans / volume
        else:
            return ans / (volume ** (self.num_species - 1))


    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):

        sp_ind_dict = {}
        sp_ind_counter = 0
        for key,value in param_dictionary.items():
            if key == 'species':
                if '+' in value or '-' in value:
                    raise SyntaxError('Plus or minus character in mass action propensity string.')
                species_names = [s.strip() for s in value.split('*')]
                for species_name in species_names:
                    if species_name == '':
                        continue
                    if species_name not in sp_ind_dict:
                        self.sp_inds.push_back(species_indices[species_name])
                        self.sp_counts.push_back(1)
                        sp_ind_dict[species_name] = sp_ind_counter
                        sp_ind_counter += 1
                    else:
                        sp_ind =sp_ind_dict[species_name]
                        self.sp_counts[sp_ind] += 1

                self.num_species = int(sum(self.sp_counts))
            elif key == 'k':
                self.k_index = parameter_indices[value]
            else:
                logging.info('Warning! Useless field for MassActionPropensity '+str(key))



    def get_species_and_parameters(self, dict fields, **keywords):
        species_list = [x.strip()   for x in fields['species'].split('*') ]
        species_list = [x for x in species_list if x != '']

        return (species_list, [ fields['k'] ])


##################################################                ####################################################
######################################              PARSING                             ##############################
#################################################                     ################################################


cdef class Term:
    cdef double evaluate(self, double *species, double *params, double time):
        raise SyntaxError('Cannot make Term base object')

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        raise SyntaxError('Cannot make Term base object')


    def py_evaluate(self, np.ndarray species, np.ndarray params, double time=0.0):
        return self.evaluate(<double*> species.data, <double*> params.data, time)

    def py_volume_evaluate(self, np.ndarray species, np.ndarray params,
                           double vol, double time=0.0):
        return self.volume_evaluate(<double*> species.data, <double*> params.data,
                                    vol, time)

# Base building blocks

cdef class ConstantTerm(Term):

    def __init__(self, double val):
        self.value = val

    cdef double evaluate(self, double *species, double *params, double time):
        return self.value
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return self.value

cdef class SpeciesTerm(Term):


    def __init__(self, unsigned ind):
        self.index = ind

    cdef double evaluate(self, double *species, double *params, double time):
        return species[self.index]
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return species[self.index]

cdef class ParameterTerm(Term):

    def __init__(self, unsigned ind):
        self.index = ind

    cdef double evaluate(self, double *species, double *params, double time):
        return params[self.index]
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return params[self.index]

cdef class VolumeTerm(Term):
    cdef double evaluate(self, double *species, double *params, double time):
        return 1.0
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return vol

# Putting stuff together

cdef class SumTerm(Term):


    def __init__(self):
        self.terms_list = []

    cdef void add_term(self,Term trm):
        self.terms.push_back(<void*> trm)
        self.terms_list.append(trm)

    cdef double evaluate(self, double *species, double *params, double time):
        cdef double ans = 0.0
        cdef unsigned i
        for i in range(self.terms.size()):
            ans += (<Term>(self.terms[i])).evaluate(species, params, time)
        return ans

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        cdef double ans = 0.0
        cdef unsigned i
        for i in range(self.terms.size()):
            ans += (<Term>(self.terms[i])).volume_evaluate(species,params,vol, time)
        return ans

cdef class ProductTerm(Term):
    def __init__(self):
        self.terms_list = []

    cdef void add_term(self,Term trm):
        self.terms.push_back(<void*> trm)
        self.terms_list.append(trm)

    cdef double evaluate(self, double *species, double *params, double time):
        cdef double ans = 1.0
        cdef unsigned i
        for i in range(self.terms.size()):
            ans *= (<Term>(self.terms[i])).evaluate(species, params,time)
        return ans

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        cdef double ans = 1.0
        cdef unsigned i
        for i in range(self.terms.size()):
            ans *= (<Term>(self.terms[i])).volume_evaluate(species,params,vol,time)
        return ans

cdef class MaxTerm(Term):
    def __init__(self):
        self.terms_list = []

    cdef void add_term(self,Term trm):
        self.terms.push_back(<void*> trm)
        self.terms_list.append(trm)

    cdef double evaluate(self, double *species, double *params, double time):
        cdef double ans = (<Term>(self.terms[0])).evaluate(species, params,time)
        cdef unsigned i
        cdef double temp = 0
        for i in range(1,self.terms.size()):
            temp =  (<Term>(self.terms[i])).evaluate(species, params,time)
            if temp > ans:
                ans = temp

        return ans

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        cdef double ans = (<Term>(self.terms[0])).volume_evaluate(species,params,vol,time)
        cdef unsigned i
        cdef double temp = 0
        for i in range(1,self.terms.size()):
            temp = (<Term>(self.terms[i])).volume_evaluate(species,params,vol,time)
            if temp > ans:
                ans = temp
        return ans

cdef class MinTerm(Term):
    def __init__(self):
        self.terms_list = []

    cdef void add_term(self,Term trm):
        self.terms.push_back(<void*> trm)
        self.terms_list.append(trm)

    cdef double evaluate(self, double *species, double *params, double time):
        cdef double ans = (<Term>(self.terms[0])).evaluate(species, params,time)
        cdef unsigned i
        cdef double temp = 0
        for i in range(1,self.terms.size()):
            temp =  (<Term>(self.terms[i])).evaluate(species, params,time)
            if temp < ans:
                ans = temp

        return ans

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        cdef double ans = (<Term>(self.terms[0])).volume_evaluate(species,params,vol,time)
        cdef unsigned i
        cdef double temp = 0
        for i in range(1,self.terms.size()):
            temp = (<Term>(self.terms[i])).volume_evaluate(species,params,vol,time)
            if temp < ans:
                ans = temp
        return ans

cdef class PowerTerm(Term):


    cdef void set_base(self, Term base):
        self.base = base
    cdef void set_exponent(self, Term exponent):
        self.exponent = exponent

    cdef double evaluate(self, double *species, double *params, double time):
        return self.base.evaluate(species,params,time) ** \
               self.exponent.evaluate(species,params,time)

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return self.base.volume_evaluate(species,params,vol,time) ** \
               self.exponent.volume_evaluate(species,params,vol,time)


cdef class ExpTerm(Term):
    cdef void set_arg(self, Term arg):
        self.arg = arg

    cdef double evaluate(self, double *species, double *params, double time):
        return exp(self.arg.evaluate(species,params,time))

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return exp(self.arg.volume_evaluate(species,params,vol,time))

cdef class LogTerm(Term):
    cdef void set_arg(self, Term arg):
        self.arg = arg

    cdef double evaluate(self, double *species, double *params, double time):
        return log(self.arg.evaluate(species,params,time))

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return log(self.arg.volume_evaluate(species,params,vol,time))


cdef class StepTerm(Term):
    cdef void set_arg(self, Term arg):
        self.arg = arg

    cdef double evaluate(self, double *species, double *params, double time):
        if self.arg.evaluate(species,params,time) >= 0:
            return 1.0
        return 0

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        if self.arg.volume_evaluate(species,params,vol,time) >= 0:
            return 1.0
        return 0

cdef class AbsTerm(Term):
    cdef void set_arg(self, Term arg):
        self.arg = arg

    cdef double evaluate(self, double *species, double *params, double time):
        return fabs( self.arg.evaluate(species,params,time) )

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return fabs( self.arg.volume_evaluate(species,params,vol,time) )


cdef class TimeTerm(Term):
    cdef double evaluate(self, double *species, double *params, double time):
        return time

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return time


def sympy_species_and_parameters(instring, species2index, params2index):
    instring = instring.replace('^','**')
    instring = instring.replace('|','_')
    root = sympy.sympify(instring, _clash1)
    nodes = [root]
    index = 0
    while index < len(nodes):
        node = nodes[index]
        index += 1
        nodes.extend(node.args)

    

    #Old Way
    #names = [str(n) for n in nodes if type(n) == sympy.Symbol]
    #species_names = [s for s in names if (s[0] != '_' and s != 'volume' and s != 't')]
    #param_names = [s[1:] for s in names if s[0] == '_']

    #New Way
    #remove leading "_" if there is one.
    names = [str(n) for n in nodes if type(n) == sympy.Symbol if str(n)[0] != "_"]+[str(n)[1:] for n in nodes if type(n) == sympy.Symbol if str(n)[0] == "_"]
    species_names = [s for s in names if s in species2index]
    param_names = [s for s in names if (s not in species2index and s != 'volume' and s != 't')]

    return species_names, param_names

def sympy_recursion(tree, species2index, params2index):
    cdef SumTerm sumterm
    cdef ProductTerm productterm
    cdef PowerTerm powerterm
    cdef ExpTerm expterm
    cdef LogTerm logterm
    cdef StepTerm stepterm
    cdef AbsTerm absterm
    cdef MaxTerm maxterm
    cdef MinTerm minterm

    root = tree.func
    args = tree.args
    # check if symbol
    if type(tree) == sympy.Symbol:
        name = str(tree)

        #remove initial underscores in names
        if name[0] == '_':
            name = name[1:]

        #New method: check params based upon being in the dictionary
        if name in species2index:
            return SpeciesTerm(species2index[ name ])
        elif name in params2index:
            return ParameterTerm(params2index[ name ])
        elif name == 'volume':
            return VolumeTerm()
        elif name == 't':
            return TimeTerm()
        else:
            raise ValueError(f"Unknown term {name} not found in Species, Parameters, or built-in-terms.")
    # check if addition
    elif type(tree) == sympy.Add:
        sumterm = SumTerm()
        for a in args:
            sumterm.add_term(  sympy_recursion(a,species2index,params2index)   )
        return sumterm

    # check multiplication
    elif type(tree) == sympy.Mul:
        productterm = ProductTerm()
        for a in args:
            productterm.add_term(sympy_recursion(a,species2index,params2index))
        return productterm

    # check exponential

    elif type(tree) == sympy.Pow:
        powerterm = PowerTerm()
        powerterm.set_base( sympy_recursion(args[0],species2index,params2index) )
        powerterm.set_exponent( sympy_recursion(args[1], species2index,params2index) )
        return powerterm

    # check exp and log
    elif type(tree) == sympy.exp:
        expterm = ExpTerm()
        expterm.set_arg( sympy_recursion(args[0],species2index,params2index) )
        return expterm

    elif type(tree) == sympy.log:
        logterm = LogTerm()
        logterm.set_arg( sympy_recursion(args[0],species2index,params2index) )
        return logterm

    # check Heaviside
    elif type(tree) == sympy.Heaviside:
        stepterm = StepTerm()
        stepterm.set_arg( sympy_recursion(args[0],species2index,params2index) )
        return stepterm

    # check absolute value

    elif type(tree) == sympy.Abs:
        absterm = AbsTerm()
        absterm.set_arg( sympy_recursion(args[0],species2index,params2index) )
        return absterm

    # check for min and max

    elif type(tree) == sympy.Max:
        maxterm = MaxTerm()
        for a in args:
            maxterm.add_term(sympy_recursion(a,species2index,params2index))
        return maxterm

    elif type(tree) == sympy.Min:
        minterm = MinTerm()
        for a in args:
            minterm.add_term(sympy_recursion(a,species2index,params2index))
        return minterm

    # if nothing else, then it should be a number

    else:
        try:
            return ConstantTerm(float( tree.evalf() ))
        except:
            raise SyntaxError('This should be a number: ' + str(tree))



def parse_expression(instring, species2index, params2index):
    instring = instring.strip()
    instring = instring.replace('^','**')
    instring = instring.replace('|', '_')
    instring = instring.replace('heaviside', 'Heaviside')

    try:
        parse_tree = sympy.sympify(instring, _clash1)
    except:
        raise SyntaxError('Sympy unable to parse: ' + instring)

    return sympy_recursion(parse_tree,species2index,params2index)


cdef class GeneralPropensity(Propensity):

    cdef double get_propensity(self, double* state, double* params, double time):
        return self.term.evaluate(state,params,time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        return self.term.volume_evaluate(state,params,volume,time)

    def __init__(self):
        self.propensity_type = PropensityType.general

    def initialize(self, dict dictionary, dict species2index, dict params2index):
        instring = dictionary['rate']
        self.term = parse_expression(instring, species2index, params2index)

    def get_species_and_parameters(self, dict fields, dict species2index, dict params2index):
        instring = fields['rate'].strip()
        return sympy_species_and_parameters(instring, species2index, params2index)



##################################################                ####################################################
######################################              DELAY TYPES                        ##############################
#################################################                     ################################################

cdef class Delay:
    def __init__(self):
        """
        Set the delay_type attribute to the appropriate enum value.
        """
        self.delay_type = DelayType.unset_delay

    def py_get_delay(self, np.ndarray[np.double_t,ndim=1] state,
                     np.ndarray[np.double_t,ndim=1] params):
        """
        Return the delay given the state and parameter vector
        :param state: (np.ndarray) the state vector
        :param params: (np.ndarray) the parameters vector
        :return: (double) the computed delay

        This function should NOT be overridden by subclases. It is just a Python
        wrapped of the cython delay function.
        """
        return self.get_delay(<double*> state.data, <double*> params.data)

    def __eq__(self, Delay other):
        return self.delay_type == other.delay_type

    cdef double get_delay(self, double* state, double* params):
        """
        Compute a delay given the state and parameters vectors.
        :param state: (double *) the array containing the state vector
        :param params: (double *) the array containing the parameters vector
        :return: (double) the computed delay.

        This function must be overridden by subclasses.
        """

        return -1.0

    def initialize(self, dict param_dictionary, dict species_indices,
                   dict parameter_indices):

        """
        Initializes the parameters and species to look at the right indices in
            the state
        :param dictionary: (dict:str--> str) the fields for the propensity 'k',
                                            's1' etc map to the actual parameter
                                             and species names
        :param species_indices: (dict:str-->int) map species names to entry in
                                species vector
        :param parameter_indices: (dict:str-->int) map param names to entry in
                                  param vector
        :return: nothing
        """
        pass

    def get_species_and_parameters(self, dict fields, **keywords):
        """
        get which fields are species and which are parameters
        :return: (list(string), list(string)) First entry is the fields that are
                                            species, second entry is the fields
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

    def initialize(self, dict param_dictionary, dict species_indices,
                   dict parameter_indices):

        for key,value in param_dictionary.items():
            if key == 'delay':
                self.delay_index = parameter_indices[value]
            else:
                logging.info('Warning! Useless field for fixed delay', key)

    def get_species_and_parameters(self, dict fields, **keywords):
        return [], [fields['delay']]

cdef class GaussianDelay(Delay):

    def __init__(self):
        self.delay_type = DelayType.gaussian

    cdef double get_delay(self, double* state, double* params):
        return cyrandom.normal_rv(params[self.mean_index],params[self.std_index])


    def initialize(self, dict param_dictionary, dict species_indices,
                   dict parameter_indices):

        for key,value in param_dictionary.items():
            if key == 'mean':
                self.mean_index = parameter_indices[value]
            elif key == 'std':
                self.std_index = parameter_indices[value]
            else:
                logging.info('Warning! Useless field for gaussian delay', key)

    def get_species_and_parameters(self, dict fields, **keywords):
        return [],[fields['mean'], fields['std']]



cdef class GammaDelay(Delay):

    def __init__(self):
        self.delay_type = DelayType.gamma

    cdef double get_delay(self, double* state, double* params):
        return cyrandom.gamma_rv(params[self.k_index],params[self.theta_index])

    def initialize(self, dict param_dictionary, dict species_indices,
                   dict parameter_indices):

        for key,value in param_dictionary.items():
            if key == 'k':
                self.k_index = parameter_indices[value]
            elif key == 'theta':
                self.theta_index = parameter_indices[value]
            else:
                logging.info('Warning! Useless field for gamma delay', key)

    def get_species_and_parameters(self, dict fields, **keywords):
        return [],[fields['k'], fields['theta']]

##################################################                ####################################################
######################################              RULE   TYPES                       ###############################
#################################################                     ################################################

cdef class Rule:
    """
    A class for doing rules that must be done either at the beginning of a simulation or repeatedly at each step of
    the simulation.
    """
    cdef void execute_rule(self, double *state, double *params, double time):
        raise NotImplementedError('Creating base Rule class. This should be subclassed.')

    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time):
        self.execute_rule(state, params, time)

    def py_execute_rule(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                        double time = 0.0):
        self.execute_rule(<double*> state.data, <double*> params.data,time)

    def py_execute_volume_rule(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                               double volume, double time=0.0 ):
        self.execute_volume_rule(<double*> state.data, <double*> params.data, volume,time)

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

    def get_species_and_parameters(self, dict fields, **keywords):
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

    cdef void execute_rule(self, double *state, double *params, double time):
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

    def get_species_and_parameters(self, dict fields, **keywords):
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
    cdef void execute_rule(self, double *state, double *params, double time):
        if self.param_flag > 0:
            params[self.dest_index] = self.rhs.evaluate(state,params,time)
        else:
            state[self.dest_index] = self.rhs.evaluate(state,params,time)

    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time):
        if self.param_flag > 0:
            params[self.dest_index] = self.rhs.volume_evaluate(state,params,volume, time)
        else:
            state[self.dest_index] = self.rhs.volume_evaluate(state,params,volume, time)

    def initialize(self, dict fields, species2index, params2index):
        self.rhs = parse_expression(fields['equation'].split('=')[1], species2index, params2index)

        dest_name = fields['equation'].split('=')[0].strip()

        #if dest_name[0] == '_' or dest_name[0] == '|':
        if dest_name[0] == '_':
            dest_name = dest_name[1:]
        if dest_name in params2index:
            self.param_flag = 1
            self.dest_index = params2index[dest_name]
        else:
            self.param_flag = 0
            self.dest_index = species2index[dest_name]

    def get_species_and_parameters(self, dict fields, dict species2index, dict params2index):
        instring = fields['equation'].strip()
        dest_name = instring.split('=')[0].strip()
        instring = instring.split('=')[1]

        species_names, param_names = sympy_species_and_parameters(instring, species2index, params2index)

        if dest_name[0] == '_' or dest_name[0] == '|':
            dest_name = dest_name[1:]

        if dest_name in species2index:
            species_names.append(dest_name)
        else:
            param_names.append(dest_name)

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

    cdef Volume copy(self):
        """
        Returns a deep copy of the volume object
        """
        raise NotImplementedError('Need to implement copy for population simulations')

    def py_copy(self):
        """
        Copy function for deep copying
        :return: a deep copy of the volume object
        """
        return self.copy()

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


    cdef Volume copy(self):
        cdef StochasticTimeThresholdVolume v = StochasticTimeThresholdVolume(self.cell_cycle_time,
                                                                             self.average_division_volume,
                                                                             self.division_noise)
        v.division_time = self.division_time
        v.current_volume = self.current_volume
        return v

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

    def __init__(self):
        pass

    def setup(self, double average_division_volume, double division_noise, growth_rate, Model m):
        self.average_division_volume = average_division_volume
        self.division_noise = division_noise
        self.growth_rate = m.parse_general_expression(growth_rate)


    cdef double get_volume_step(self, double *state, double *params, double time, double volume, double dt):
        cdef double gr = self.growth_rate.evaluate(state,params, time)
        return ( exp(gr*dt) - 1.0) * volume

    cdef void initialize(self, double *state, double *params, double time, double volume):
        self.py_set_volume(volume)
        # Must choose division volume.
        self.division_volume = self.average_division_volume * cyrandom.normal_rv(1.0, self.division_noise)
        if self.division_noise > volume:
            raise RuntimeError('Division occurs before initial volume - change your parameters!')


    cdef unsigned cell_divided(self, double *state, double *params, double time, double volume, double dt):
        if volume > self.division_volume:
            return 1
        return 0

    cdef Volume copy(self):
        cdef StateDependentVolume sv = StateDependentVolume()
        sv.division_noise = self.division_noise
        sv.division_volume = self.division_volume
        sv.growth_rate = self.growth_rate
        sv.average_division_volume = self.average_division_volume
        sv.current_volume = self.current_volume
        return sv


###############################                #################################
###################              MODEL   TYPES                       ###########
##############################                     #############################

cdef class Model:
    def __init__(self, filename = None, species = [], reactions = [], parameters = [], rules = [], 
                initial_condition_dict = None, sbml_filename = None, input_printout = False, 
                initialize_model = True, **kwargs):
        """
        Read in a model from a file using XML format for the model.

        :param filename: (str) the file to read the model
        """
        self._next_species_index = 0
        self._next_params_index = 0
        self._dummy_param_counter = 0

        self.has_delay = False #Does the Model contain any delay reactions? 
                               #Updated in _add_reaction.

        self.species2index = {}
        self.params2index = {}
        self.propensities = []
        self.delays = []
        self.repeat_rules = []
        self.params_values = np.array([])
        self.species_values = np.array([])
        self.txt_dict = {'reactions':"", 'rules':""} # A dictionary to store XML 
                                                     #txt to write bioscrape xml
        self.reaction_definitions = [] # List of reaction tuples useful for writing SBML
        self.rule_definitions = [] #A list of rule tuples useful for writing SBML

        # These must be updated later
        self.update_array = None
        self.delay_update_array = None
        self.reaction_updates = []
        self.delay_reaction_updates = []
        # Set to True when the stochiometric matrices are created and model 
        # checked by the initialize() function
        self.initialized = False 
        self.reaction_list = [] # A list used to store tuples (propensity, 
                                # delay, update_array, delay_update_array) for 
                                # each reaction

        if filename != None and sbml_filename != None:
            raise ValueError("Cannot load both a bioSCRAPE xml file and an " 
                             "SBML file. Please choose just one.")
        elif filename != None:
            self.parse_model(filename, input_printout = input_printout)
        elif sbml_filename != None:
            import_sbml(sbml_filename, bioscrape_model = self, input_printout = input_printout, **kwargs)

        for species in species:
            self._add_species(species)

        for rxn in reactions:
            if len(rxn) == 4:
                reactants, products, propensity_type, propensity_param_dict = rxn
                delay_type, delay_reactants, delay_products, delay_param_dict =\
                        None, None,  None, None
            elif len(rxn) == 8:
                reactants, products, propensity_type, propensity_param_dict, \
                delay_type, delay_reactants, delay_products, delay_param_dict = rxn
            else:
                raise ValueError("Reaction Tuple of the wrong length! Must be "
                                 "of length 4 (no delay) or 8 (with delays). "
                                 "See BioSCRAPE Model API for details.")
            self.create_reaction(reactants, products, propensity_type, 
                                 propensity_param_dict, delay_type, 
                                 delay_reactants, delay_products, 
                                 delay_param_dict, 
                                 input_printout = input_printout)

        if isinstance(parameters, dict):
            parameters = parameters.items()

        for param, param_val in parameters:
                self._add_param(param)
                self.set_parameter(param, param_val)


        for rule in rules:
            if len(rule) == 2:
                rule_type, rule_attributes = rule
                self.create_rule(rule_type, rule_attributes, 
                                 input_printout = input_printout)
            elif len(rule) == 3:
                rule_type, rule_attributes, rule_frequency = rule
                self.create_rule(rule_type, rule_attributes, 
                                 rule_frequency = rule_frequency, 
                                 input_printout = input_printout)
            else:
                raise ValueError("Rules must be a tuple: (rule_type (string), "
                                 "rule_attributes (dict), rule_frequency "
                                 "(optional))")

        if initial_condition_dict != None:
            for specie in initial_condition_dict:
                self._add_species(specie)
            self.set_species(initial_condition_dict)

        if initialize_model:
            self._initialize()

    cdef void _initialize(self):
        #creates C vector objects
        self._create_vectors()

        #Create Stochiometric Matrices
        self._create_stochiometric_matrices()

        #Check for unspecified parameters
        self.check_parameters()

        #Check for species without intial conditions.
        #Set these initial conditions to 0 and issue a warning.
        self.check_species()

        self.initialized = True

    def py_initialize(self):
        self._initialize()

    def _create_vectors(self):
        #Create c-vectors of different objects
        self.propensities = []
        self.c_propensities.clear()
        self.delays = []
        self.c_delays.clear()
        for rxn in self.reaction_list:
            prop_object, delay_object, update_array, delay_update_array = rxn
            self.propensities.append(prop_object)
            self.c_propensities.push_back(<void*> prop_object)
            self.delays.append(delay_object)
            self.c_delays.push_back(<void*> delay_object)

        self.c_repeat_rules.clear()
        for rule_object in self.repeat_rules:
            self.c_repeat_rules.push_back(<void*> rule_object)

    def py_initialize(self):
        self._initialize()

    def __eq__(self, Model other):
        if other is None or not isinstance(other, Model):
            return False
        # Casting as a set means order doesn't matter.
        # Sets can only hold an element once, so this could give weird results
        # if the same reaction or rule definition appears multiple times.
        # 
        # If reaction/rule definitions are the same, that implies that many of 
        # the other attributes of the Model must be the same. 
        if sorted(self.reaction_definitions) != sorted(other.reaction_definitions):
            return False
        if sorted(self.rule_definitions) != sorted(other.rule_definitions):
            return False
        if self.species2index != other.species2index:
            return False
        if not np.array_equal(self.species_values, other.species_values):
            return False
        return True

    def __neq__(self, Model other):
        return not self.__eq__(other)

    def _add_species(self, species):
        """
        Helper function for putting together the species vector (converting species names to indices in vector)

        If the species has already been added, then do nothing. otherwise give it a new index, and increase
        the next_species_index by 1

        :param species: (str) the species name
        :return: None
        """
        self.initialized = False
        if species not in self.species2index and species is not None:
            self.species2index[species] = self._next_species_index
            self._next_species_index += 1
            self.species_values = np.concatenate((self.species_values, np.array([-1])))

    def _set_species_value(self, specie, value):
        if specie not in self.species2index:
            self._add_species(specie)
        self.species_values[self.species2index[specie]] = value

    #Helper function to add a reaction to the model
    #Inputs:
    #   reaction_update_dict (dictionary): species_index --> change in count. Species not in the products or reactants can be omitted
    #   propensity_object: an instance of a propensity_object
    #   propensity_param_dict: a dictionary containing the parameters of the propensity
    #   delay_reaction_update_dict: same as reaction_dict but for the delayed part of a reaction
    #   delay_object: an instance of one of a delay_object
    #   delay_param_dict: a dictionary containing the parameters of the delay distribution

    def _add_reaction(self, reaction_update_dict, propensity_object, propensity_param_dict,
        delay_reaction_update_dict = {}, delay_object = None, delay_param_dict = {}):
        self.initialized = False

        species_names, param_names = propensity_object.get_species_and_parameters(propensity_param_dict, species2index = self.species2index, params2index = self.params2index)

        for species_name in species_names:
            #self._add_species(species_name)
            #Now no species should be added here
            pass
        for param_name in param_names:
            self._add_param(param_name)

        self.reaction_updates.append(reaction_update_dict)
        propensity_object.initialize(propensity_param_dict, self.species2index, self.params2index)

        if delay_object == None:
           delay_object = NoDelay()
        elif not type(delay_object) == type(NoDelay()):
            self.has_delay = True

        species_names, param_names = delay_object.get_species_and_parameters(delay_param_dict, species2index = self.species2index, params2index = self.params2index)

        for species_name in species_names:
            #self._add_species(species_name)
            #Now anything not declared as a Species will be interpreted as a parameter
            pass
        for param_name in param_names:
            self._add_param(param_name)

        #Moved to Model._initialize
        #self.delays.append(delay_object)
        #self.c_delays.push_back(<void*> delay_object)
        self.delay_reaction_updates.append(delay_reaction_update_dict)
        delay_object.initialize(delay_param_dict, self.species2index, self.params2index)
        self.reaction_list.append((propensity_object, delay_object, reaction_update_dict, delay_reaction_update_dict))


    def create_propensity(self, propensity_type, propensity_param_dict, print_out = False):
        if print_out:
            warnings.warn("Creating Propensity: prop_type="+str(propensity_type)+" params="+str(propensity_param_dict))
        if 'type' in propensity_param_dict:
            propensity_param_dict.pop('type')
        #Create propensity object
        if propensity_type == 'hillpositive':
            #Check required propensity parameters and convert numeric parameters to dummy variables.
            self._param_dict_check(propensity_param_dict, "k", "DummyVar_PositiveHillPropensity")
            self._param_dict_check(propensity_param_dict, "K", "DummyVar_PositiveHillPropensity")
            self._param_dict_check(propensity_param_dict, "s1", "DummyVar_PositiveHillPropensity")
            self._param_dict_check(propensity_param_dict, "n", "DummyVar_PositiveHillPropensity")
            prop_object = PositiveHillPropensity()


        elif propensity_type == 'proportionalhillpositive':
            self._param_dict_check(propensity_param_dict, "k", "DummyVar_PositiveProportionalHillPropensity")
            self._param_dict_check(propensity_param_dict, "K", "DummyVar_PositiveProportionalHillPropensity")
            self._param_dict_check(propensity_param_dict, "s1", "DummyVar_PositiveProportionalHillPropensity")
            self._param_dict_check(propensity_param_dict, "d", "DummyVar_PositiveProportionalHillPropensity")
            self._param_dict_check(propensity_param_dict, "n", "DummyVar_PositiveProportionalHillPropensity")
            prop_object = PositiveProportionalHillPropensity()

        elif propensity_type == 'hillnegative':
            self._param_dict_check(propensity_param_dict, "k", "DummyVar_NegativeHillPropensity")
            self._param_dict_check(propensity_param_dict, "K", "DummyVar_NegativeHillPropensity")
            self._param_dict_check(propensity_param_dict, "s1", "DummyVar_NegativeHillPropensity")
            self._param_dict_check(propensity_param_dict, "n", "DummyVar_NegativeHillPropensity")
            prop_object = NegativeHillPropensity()

        elif propensity_type == 'proportionalhillnegative':
            self._param_dict_check(propensity_param_dict, "k", "DummyVar_NegativeProportionalHillPropensity")
            self._param_dict_check(propensity_param_dict, "K", "DummyVar_NegativeProportionalHillPropensity")
            self._param_dict_check(propensity_param_dict, "s1", "DummyVar_NegativeProportionalHillPropensity")
            self._param_dict_check(propensity_param_dict, "d", "DummyVar_NegativeProportionalHillPropensity")
            self._param_dict_check(propensity_param_dict, "n", "DummyVar_NegativeProportionalHillPropensity")
            prop_object = NegativeProportionalHillPropensity()

        elif propensity_type == 'massaction':
            species_string = propensity_param_dict['species']

            # if mass action propensity has less than 3 things, then use consitutitve, uni, bimolecular for speed.
            if species_string in ["0", "", '', None, 0]:
                prop_object = ConstitutivePropensity()
                self._param_dict_check(propensity_param_dict, "k", "DummyVar_ConstitutivePropensity")
            else:
                species_names = species_names = [x.strip() for x in species_string.split('*') if x.strip() not in ["", '']]

                if len(species_names) == 1:
                    prop_object = UnimolecularPropensity()
                    self._param_dict_check(propensity_param_dict, "k", "DummyVar_UnimolecularPropensity")

                elif len(species_names) == 2:
                    prop_object = BimolecularPropensity()
                    self._param_dict_check(propensity_param_dict, "k", "DummyVar_BimolecularPropensity")

                else:
                    prop_object = MassActionPropensity()
                    self._param_dict_check(propensity_param_dict, "k", "DummyVar_MassActionPropensity")

        elif propensity_type == 'general':
            prop_object = GeneralPropensity()
        else:
            raise SyntaxError('Propensity Type is not supported: ' + propensity_type)

        return prop_object
    #A function to programatically create a reaction (and add automatically add it to the model).
    #   Supports all native propensity types and delay types.
    #Required Inputs:
    #   reactants (list): a list of reactant specie names (strings)
    #   products (list): a list of product specie names (strings)
    #   propensity_type: a string indicating the type of propensity
    #       Supported types: "massaction", "hillpositive", "proportionalhillpositive", "hillnegative", "proportionalhillnegative", "general"
    #   propensity_param_dict: a dictionary of parameters for the given propensity type
    #Optional Inputs:
    #   delay_type: a string indicating the type of delay
    #   delay_reactants (list): a list of delay reaction reactant specie names (strings)
    #   delay_products: a list of delay reaction products specie names (strings)
    #   delay_param_dict: a dictionary of the parameters for the delay distribution
    def create_reaction(self, reactants, products, propensity_type, 
                        propensity_param_dict, delay_type = None, 
                        delay_reactants = None, delay_products = None, 
                        delay_param_dict = None, input_printout = False):

        if input_printout:
            warnings.warn("creating reaction with:"+
                "\n\tPropensity_type="+str(propensity_type)+" Inputs="+str(reactants)+" Outputs="+str(products)+
                "\n\tpropensity_param_dict="+str(propensity_param_dict)+
                "\n\tDelay_type="+str(delay_type)+" delay inputs ="+str(delay_reactants)+" delay outputs="+str(delay_products)+
                "\n\tdelay_param_dict="+str(delay_param_dict))
        self.initialized = False

        #Copy dictionaries so they aren't altered if they are being used by external code
        propensity_param_dict = dict(propensity_param_dict)

        #Remove "propensity_type" key which may be leftover from SBML annotations
        if 'propensity_type' in propensity_param_dict:
            propensity_param_dict.pop("propensity_type")

        if delay_param_dict != None:
            delay_param_dict = dict(delay_param_dict)


        #Reaction Reactants and Products stored in a dictionary
        reaction_update_dict = {}
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

        if 'species' not in propensity_param_dict and propensity_type == "massaction":
                reactant_string = ""
                for s in reactants:
                    if s is not None:
                        reactant_string += s+"*"
                propensity_param_dict['species'] = reactant_string[:len(reactant_string)-1]

        prop_object = self.create_propensity(propensity_type, propensity_param_dict, print_out = input_printout)

        #Create Delay Object
        #Delay Reaction Reactants and Products Stored in a Dictionary
        delay_reaction_update_dict = {}
        if delay_reactants != None:
            for r in delay_reactants:
                # if the species hasn't been seen add it to the index
                self._add_species(r)
                # update the update array
                if r not in delay_reaction_update_dict:
                    delay_reaction_update_dict[r] = 0
                delay_reaction_update_dict[r]  -= 1
        else:
            delay_reactants = []

        if delay_products != None:
            for p in delay_products:
                # if the species hasn't been seen add it to the index
                self._add_species(p)
                # update the update array
                if p not in delay_reaction_update_dict:
                    delay_reaction_update_dict[p] = 0
                delay_reaction_update_dict[p]  += 1
        else:
            delay_products = []


        if delay_type == 'none' or delay_type == None:
            delay_object = NoDelay()
            delay_param_dict = {}
        elif delay_type == 'fixed':
            self._param_dict_check(delay_param_dict, "delay", "DummyVar_FixedDelay")
            delay_object = FixedDelay()
        elif delay_type == 'gaussian':
            self._param_dict_check(delay_param_dict, "mean", "DummyVar_GaussianDelay")
            self._param_dict_check(delay_param_dict, "std", "DummyVar_GaussianDelay")
            delay_object = GaussianDelay()
        elif delay_type == 'gamma':
            self._param_dict_check(delay_param_dict, "k", "DummyVar_GammaDelay")
            self._param_dict_check(delay_param_dict, "theta", "DummyVar_GammaDelay")
            delay_object = GammaDelay()
        else:
            raise SyntaxError('Unknown delay type: ' + delay_type)
        delay_param_dict.pop('type',None)

        self._add_reaction(reaction_update_dict, prop_object, propensity_param_dict, delay_reaction_update_dict, delay_object, delay_param_dict)
        self.write_rxn_txt(reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict)
        self.reaction_definitions.append((reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict))
    
    def write_rxn_txt(self, reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict):
        #Write bioscrape XML and save it to the xml dictionary
        rxn_txt = '<reaction text= "'
        for r in reactants:
            if r is not None:
                rxn_txt += r +" + "
        if len(reactants)>0:
            rxn_txt = rxn_txt[:-2]
        rxn_txt += "-- "
        for p in products:
            if p is not None:
                rxn_txt += p+" + "
        if len(products)>0:
            rxn_txt = rxn_txt[:-2]
        rxn_txt +='"'
        if len(delay_reactants) > 0 or len(delay_products)> 0:
            rxn_txt += ' after= "'
            if len(delay_reactants) > 0:
                for r in delay_reactants:
                    if r is not None:
                        rxn_txt += r +" + "
                if len(delay_reactants) > 0:
                    rxn_txt = rxn_txt[:-2]
            rxn_txt += "-- "
            if len(delay_products)> 0:
                for p in delay_products:
                    if p is not None:
                        rxn_txt += p+" + "
                if len(delay_products)>0:
                    rxn_txt = rxn_txt[:-2]
                rxn_txt +='"'
        rxn_txt += '>\n\t<propensity type="'
        rxn_txt += propensity_type+'" '
        for k in propensity_param_dict:
            rxn_txt+=k+'="'+str(propensity_param_dict[k])+'" '
        rxn_txt += '/>\n\t<delay type="'
        if delay_type == None:
            rxn_txt += 'none" />'
        else:
            rxn_txt += delay_type+'" '
            for k in delay_param_dict:
                rxn_txt += 'k="'+str(delay_param_dict[k])+'" '
            rxn_txt+='/>'
        rxn_txt += '\n</reaction>\n'
        self.txt_dict['reactions']+=rxn_txt



    def _add_param(self, param_name):
        """
        Helper function for putting together the parameter vector (converting parameter names to indices in vector)

        If the parameter name has already been seen, then do nothing. Otherwise, give it a new index, and increase the
        next_params_index by 1.

        :param param: (str) the parameter name
        :return: None
        """
        self.initialized = False
        if param_name in self.species2index:
            raise ValueError(f"param_name {param_name} is the same as the name of a species!")
        if param_name not in self.params2index:
            self.params2index[param_name] = self._next_params_index
            self._next_params_index += 1
            self.params_values = np.concatenate((self.params_values, np.array([np.nan])))

    #Creates a rule and adds it to the model.
    #Inputs:
    #   rule_type (str): The type of rule. Supported: "additive" and "assignment"
    #   rule_attributes (dict): A dictionary of rule parameters / attributes.
    #       NOTE: the only attributes used by additive/assignment rules are 'equation'
    #   rule_frequency: must be 'repeated'
    #Rule Types Supported:
    def create_rule(self, rule_type, rule_attributes, rule_frequency = "repeated", input_printout = False):
        if input_printout:
            warnings.warn("Rule Created with \n\trule_type = "+str(rule_type)+"\n\trule_attributes="+str(rule_attributes)+"\n\trule_frequence="+str(rule_frequency))

        self.initialized = False

        # Parse the rule by rule type
        if rule_type == 'additive':
            rule_object = AdditiveAssignmentRule()
        elif rule_type == 'assignment':
            rule_object = GeneralAssignmentRule()
        else:
            raise SyntaxError('Invalid type of Rule: ' + rule_type)

        # Add species and params to model
        species_names, params_names = rule_object.get_species_and_parameters(rule_attributes, species2index = self.species2index, params2index=self.params2index)
        for s in species_names: pass #self._add_species(s) No species should be added here
        for p in params_names: self._add_param(p)

        # initialize the rule
        if 'type' in rule_attributes:
            rule_attributes.pop('type')
        rule_object.initialize(rule_attributes,self.species2index,self.params2index)
        # Add the rule to the right place
        if rule_frequency == 'repeated':
            self.repeat_rules.append(rule_object)
        else:
            raise SyntaxError('Invalid Rule Frequency: ' + str(rule_frequency))


        self.write_rule_txt(rule_type, rule_attributes, rule_frequency)
        self.rule_definitions.append((rule_type, rule_attributes, rule_frequency))

    def write_rule_txt(self, rule_type, rule_attributes, rule_frequency):
        rule_txt = '<rule type="'+rule_type+'" frequency="'+rule_frequency+'" '
        for k in rule_attributes:
            rule_txt += k+'="'+rule_attributes[k]+'" '

        rule_txt += " />\n"
        self.txt_dict["rules"]+=rule_txt

    #Sets the value of a parameter in the model
    def set_parameter(self, param_name, param_value):
        if param_name not in self.params2index:
            logging.info('Warning! parameter '+ param_name+" does not show up in any currently defined reactions or rules.")
            self._add_param(param_name)

        param_index = self.params2index[param_name]
        self.params_values[param_index] = param_value

    def create_parameter(self, param_name, param_value):
        self._add_param(param_name)
        self.set_parameter(param_name, param_value)

    #Checks that all parameters have values
    def check_parameters(self):
        error_string = "Unspecified Parameters: "
        unspecified_parameters = False
        for p in self.params2index:
            i = self.params2index[p]
            if np.isnan(self.params_values[i]):
                unspecified_parameters = True
                error_string += p+"="+str(self.params_values[i])+', '

        if unspecified_parameters:
            raise ValueError(error_string[:-2])

    #Checks that species' values are all set. Unset values default to 0 and warning is raised.
    def check_species(self):
        uninitialized_species = False
        warning_txt = "The following species are uninitialized and their value has been defaulted to 0: "
        for s in self.species2index.keys():
            i = self.species2index[s]
            if self.species_values[i] == -1:
                uninitialized_species = True
                warning_txt += s+", "
                self.species_values[i] = 0
        if uninitialized_species:
            logging.info(warning_txt)

    #Checks if the dictionary dic contains the keyword key.
    #if dic[key] = str: do nothing
    #if dic[key] = float (or a string that can be cast to a float without an error):
    #   create a dummy parameter and set its value to float then set dict[key] = dummy_param
    def _param_dict_check(self, dic, key, param_object_name):
        if key not in dic:
            raise ValueError("param dictionary does not contain required key: "+str(key)+" for param object "+param_object_name)
        else:
            try:
                val = float(dic[key])
                float_val = True
            except ValueError:
                float_val = False

            if float_val:
                dummy_var = param_object_name+"_"+str(key)+"_"+str(self._dummy_param_counter)

                if dummy_var in self.params2index:
                    raise ValueError("Trying to create a dummy parameter that already exists. Dummy Param Name: "+dummy_var+". Please don't name your parameters like this to avoid errors.")
                self._add_param(dummy_var)
                self.set_parameter(dummy_var, val)
                dic[key] = dummy_var
                self._dummy_param_counter = self._dummy_param_counter + 1

    #Helper Function to Create Stochiometric Matrices for Reactions and Delay Reactions
    def _create_stochiometric_matrices(self):
        # With all reactions read in, generate the update array
        num_species = len(self.species2index.keys())
        num_reactions = len(self.reaction_list)
        self.update_array = np.zeros((num_species, num_reactions))
        self.delay_update_array = np.zeros((num_species,num_reactions))
        for reaction_index in range(num_reactions):
            prop_object, delay_object, reaction_update_dict, delay_reaction_update_dict = self.reaction_list[reaction_index]
            #reaction_update_dict = self.reaction_updates[reaction_index]
            #delay_reaction_update_dict = self.delay_reaction_updates[reaction_index]
            for sp in reaction_update_dict:
                self.update_array[self.species2index[sp],reaction_index] = reaction_update_dict[sp]

            for sp in delay_reaction_update_dict:
                self.delay_update_array[self.species2index[sp],reaction_index] = delay_reaction_update_dict[sp]

        return self.update_array, self.delay_update_array

    def parse_model(self, filename, input_printout = False):
        """
        Parse the model from an XML file filling in all the local variables (propensities, delays, update arrays). Also
        maps the species and parameters to indices in a species and parameters vector.

        :param filename: (str or file) the model file. if a string, the file is opened. otherwise, it is assumed
                         that a file handle was passed in.
        :return: None
        """
        # open XML file from the filename and use BeautifulSoup to parse it
        warnings.warn("Depricated Warning: Bioscrape XML is being replaced by SBML and will no longer be supported in a future version of the software.")

        if type(filename) == str:
            xml_file = open(filename,'r')
        else:
            xml_file = filename
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


        # check for model tag at beginning.

        Model = xml.find_all('model')
        if len(Model) != 1:
            raise SyntaxError('Did not include global model tag in XML file')

        Species = xml.find_all('species')
        for species in Species:
            species_value = float(species['value'])
            species_name = species['name']
            self._set_species_value(species_name, species_value)

        Reactions = xml.find_all('reaction')
        for reaction in Reactions:
            # Parse the stoichiometry
            text = reaction['text']
            reactants = [s for s in [r.strip() for r in text.split('--')[0].split('+')] if s]
            products = [s for s in [r.strip() for r in text.split('--')[1].split('+')] if s]

            for s in reactants + products:
                if s not in self.species2index:
                    raise ValueError(f"Species {s} found in a reaction but not declared in Species. All Species must be declared for proper parsing.")

            # parse the delayed part of the reaction the same way as we did before.
            if reaction.has_attr('after'):
                text = reaction['after']
                delay_reactants = [s for s in [r.strip() for r in text.split('--')[0].split('+')] if s]
                delay_products = [s for s in [r.strip() for r in text.split('--')[1].split('+')] if s]
            else:
                delay_reactants = None
                delay_products = None

            # Then look at the propensity and set up a propensity object
            propensity = reaction.find_all('propensity')
            if len(propensity) != 1:
                raise SyntaxError('Incorrect propensity tags in XML model\n' + propensity)
            propensity = propensity[0]
            # go through propensity types
            propensity_param_dict = propensity.attrs

            # Then look at the delay and set up a delay object
            delay = reaction.find_all('delay')
            if len(delay) != 1:
                raise SyntaxError('Incorrect delay spec')
            delay = delay[0]
            delay_param_dict = delay.attrs
            delay_type = delay['type']

            self.create_reaction(reactants = reactants, products = products, propensity_type = propensity['type'], propensity_param_dict = propensity_param_dict,
                delay_reactants=delay_reactants, delay_products=delay_products, delay_param_dict = delay_param_dict, input_printout = input_printout)


        # Parse through the rules
        Rules = xml.find_all('rule')
        for rule in Rules:
            rule_attrs = rule.attrs
            rule_type = rule['type']
            rule_frequency = rule['frequency']
            self.create_rule(rule_type = rule_type, rule_attributes = rule_attrs, rule_frequency=rule_frequency, input_printout = input_printout)

        # Generate species values and parameter values
        unspecified_param_names = set(self.params2index.keys())
        Parameters = xml.find_all('parameter')
        for param in Parameters:
            param_value = float(param['value'])
            param_name = param['name']
            self.set_parameter(param_name = param_name, param_value = param_value)


    def get_params2index(self):
        return self.params2index

    def get_species2index(self):
        return self.species2index

    def get_species_list(self):
        l = [None] * self.get_number_of_species()
        for s in self.species2index:
            l[self.species2index[s]] = s
        return l

    def get_species_array(self):
        return np.array(self.species_values)

    def get_param_list(self):
        l = [None] * self.get_number_of_params()
        for p in self.params2index:
            l[self.params2index[p]] = p
        return l

    def get_params(self):
        """
        Get the set of parameter names.
        :return: (dict_keys str) the parameter names
        """

        return self.params2index.keys()

    def get_species(self):
        """
        Get the set of species names.
        :return: (dict_keys str) the parameter names
        """
        return self.species2index.keys()

    def get_number_of_params(self):
        return len(self.params2index.keys())

    def get_parameter_values(self):
        return self.params_values

    def get_parameter_dictionary(self):
        param_dict = {}
        keys = self.get_params()
        values = self.get_parameter_values()
        for (key, value) in zip(keys, values):
            param_dict[key] = value
        return param_dict

    def get_species(self):
        """
        Get the set of species names.
        :return: (dict_keys str) the species names
        """

        return self.species2index.keys()

    def get_species_dictionary(self):
        """
        Get a dictionary {"species":value}
        """
        A = self.get_species_array()
        species_dict = {}
        for s in self.species2index:
            species_dict[s] = A[self.species2index[s]]
        return species_dict
        # return {(s, A[self.species2index[s]]) for s in self.species2index}

    def get_number_of_species(self):
        return len(self.species2index.keys())


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
                warnings.warn('Trying to set parameter that is not in model: %s'  % p)


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
                warnings.warn('Trying to set species that is not in model: %s' % s)

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

    def get_reactions(self):
        return self.reaction_list

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


    def parse_general_expression(self, instring):
        return parse_expression(instring,self.species2index,self.params2index)


    def write_bioscrape_xml(self, file_name):
        warnings.warn("Depricated Warning: Bioscrape XML is being replaced by SBML and will no longer be supported in a future version of the software.")
        #Writes Bioscrape XML
        txt = "<model>\n"
        species = self.get_species_list()

        #Write the Species
        for s in species:
            v = self.get_species_value(s)
            txt+='<species name="'+s+'" value="'+str(v)+'" />\n'
        txt+='\n'
        parameters = self.get_param_list()
        for p in parameters:
            v = self.get_param_value(p)
            txt+='<parameter name="'+p+'" value="'+str(v)+'" />\n'
        txt+='\n'
        txt += self.txt_dict["reactions"]
        txt+='\n'
        txt += self.txt_dict["rules"]
        txt += "</model>"

        f = open(file_name, 'w')
        f.write(txt)
        f.close()

    #Generates an SBML Model
    def generate_sbml_model(self, stochastic_model = False, **keywords):
        if self.has_delay:
            raise NotImplementedError("Writing SBML for bioscrape models with delay has not been implemented.")

        # Create an empty SBMLDocument object to hold the bioscrape model
        document, model = create_sbml_model(**keywords)

        sorted_params = list(self.get_param_list())
        sorted_params.sort()
        for p in sorted_params:
            val = self.get_param_value(p)
            if p[0] == '_':
                # Remove the underscore at the beginning of the parameter name
                p = p.replace('_','',1)
            add_parameter(model = model, param_name=p, param_value = val)

        sorted_species = list(self.get_species())
        sorted_species.sort()
        for s in sorted_species:
            add_species(model = model, compartment=model.getCompartment(0),
                        species=s, initial_concentration=self.get_species_value(s))

        rxn_count = 0
        for rxn_tuple in self.reaction_definitions:
            rxn_id = "r" + str(rxn_count)

            (reactants, products, propensity_type, propensity_param_dict,
             delay_type, delay_reactants, delay_products, delay_param_dict) = rxn_tuple

            add_reaction(model, reactants, products, rxn_id, propensity_type,
                         propensity_param_dict, stochastic = stochastic_model)
            rxn_count += 1

        rule_count = 0
        for rule_tuple in self.rule_definitions:
            rule_id = "rule" + str(rule_count)
            # Syntax of rule_tuple = (rule_type, rule_dict, rule_frequency)
            (rule_type, rule_dict, rule_frequency) = rule_tuple
            # Extract the rule variable id from rule_dict:
            equation = rule_dict['equation']
            split_eqn = [s.strip() for s in equation.split('=') ]
            assert(len(split_eqn) == 2) # Checking rule_dict equation structure.
            # Extract the rule formula for the variable above from rule_dict:
            rule_formula = split_eqn[1]
            rule_variable = split_eqn[0]
            add_rule(model, rule_id, rule_type, rule_variable, rule_formula)
            rule_count += 1

        if document.getNumErrors():
            warnings.warn('SBML model generated has errors. Use document.getErrorLog() to print all errors.')
        return document, model

    #write an SBML Model
    def write_sbml_model(self, file_name, stochastic_model = False, **keywords):
        document, _ = self.generate_sbml_model(stochastic_model = stochastic_model, **keywords)
        sbml_string = libsbml.writeSBMLToString(document)
        with open(file_name, 'w') as f:
            f.write(sbml_string)
        return True

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

    def py_set_data(self, data):
        self.data = data

    def py_get_time(self):
        return self.time

    def py_get_volume(self):
        return self.volume

    def py_get_dataframe(self, Model = None):
        try:
            import pandas
            if Model == None:
                warnings.warn("No Model passed into py_get_dataframe. No species names will be attached to the data frame.")
                df = pandas.DataFrame(data = self.py_get_data())
            else:
                columns = Model.get_species_list()
                df = pandas.DataFrame(data = self.py_get_data(), columns = columns)
            df['time'] = self.time
            df['volume'] = self.volume
            return df

        except ModuleNotFoundError:
            warnings.warn("py_get_dataframe requires the pandas Module to return a Pandas Dataframe object. Numpy array being returned instead.")
            return self.py_get_result()

    def py_get_parent(self):
        return self.parent

    def py_get_daughters(self):
        return (self.daughter1, self.daughter2)


    def py_set_parent(self, Schnitz p):
        self.set_parent(p)

    def py_set_daughters(self,Schnitz d1, Schnitz d2):
        self.set_daughters(d1,d2)


    def get_sub_lineage(self, dict species_dict = None):
        cdef Lineage out
        if species_dict is None:
            out = Lineage()
        else:
            out = ExperimentalLineage(species_dict.copy())


        cdef list schnitzes_to_add = [self]
        cdef unsigned index = 0
        cdef Schnitz s = None


        while index < len(schnitzes_to_add):
            s = schnitzes_to_add[index]
            if s.get_daughter_1() is not None:
                schnitzes_to_add.append(s.get_daughter_1())
            if s.get_daughter_2() is not None:
                schnitzes_to_add.append(s.get_daughter_2())

            index += 1

        for index in range(len(schnitzes_to_add)):
            out.add_schnitz(schnitzes_to_add[index])

        return out


    def copy(self):
        cdef Schnitz temp = Schnitz(self.time.copy(),self.data.copy(),self.volume.copy())
        temp.daughter1 = self.daughter1
        temp.daughter2 = self.daughter2
        temp.parent = self.parent
        return temp

    def __setstate__(self,state):
        self.parent = state[0]
        self.daughter1 = state[1]
        self.daughter2 = state[2]
        self.time = state[3]
        self.volume = state[4]
        self.data = state[5]

    def __getstate__(self):
        return (self.parent,self.daughter1,self.daughter2, self.time, self.volume, self.data)

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

    def py_add_schnitz(self, Schnitz s):
        self.add_schnitz(s)

    def __setstate__(self, state):
        self.schnitzes = []
        self.c_schnitzes.clear()
        for s in state[0]:
            self.add_schnitz(s)

    def __getstate__(self):
        return (self.schnitzes,)

    def truncate_lineage(self,double start_time, double end_time):
        cdef Schnitz sch, new_sch
        cdef dict sch_dict = {}
        cdef Lineage new_lineage = Lineage()
        cdef np.ndarray newtime, newvolume, newdata, indices_to_keep

        sch_dict[None] = None
        for i in range(self.size()):
            sch = self.get_schnitz(i)

            newtime = sch.get_time().copy()
            newvolume = sch.get_volume().copy()
            newdata = sch.get_data().copy()

            # if the final time of this lineage is before the starting time
            # or the first time is before the end time, then skip it
            if newtime[newtime.shape[0]-1] < start_time or newtime[0] > end_time:
                sch_dict[sch] = None
                continue

            indices_to_keep = (sch.get_time() >= start_time) & (sch.get_time() <= end_time)
            newtime = newtime[indices_to_keep]
            newvolume = newvolume[indices_to_keep]
            newdata = newdata[indices_to_keep]

            sch_dict[sch] = Schnitz(newtime, newdata, newvolume)

        for i in range(self.size()):
            sch = self.get_schnitz(i)
            if sch_dict[sch] is not None:
                new_lineage.add_schnitz(sch_dict[sch])
                sch_dict[sch].py_set_parent( sch_dict[sch.get_parent()] )
                sch_dict[sch].py_set_daughters( sch_dict[sch.get_daughter_1()] , sch_dict[sch.get_daughter_2()] )

        return new_lineage

    #Returns the same tree as above
    def get_schnitzes_by_generation(self):
        #print("Creating Schnitz Tree")
        sch_tree = [[]]
        sch_tree_length = 1
        for i in range(self.py_size()):
            sch = self.get_schnitz(i)
            if sch.py_get_parent() is None:
                sch_tree[0].append(sch)
            else:
                for j in range(len(sch_tree)):
                    parent = sch.py_get_parent()
                    if parent in sch_tree[j]:
                        if len(sch_tree) <= j+1:
                            sch_tree.append([])
                            sch_tree_length += 1
                        sch_tree[j+1].append(sch)
        return sch_tree


cdef class ExperimentalLineage(Lineage):
    def __init__(self, dict species_indices={}):
        super().__init__()
        self.species_dict = species_indices

    def py_set_species_indices(self, dict species_indices):
        self.species_dict = species_indices.copy()

    def py_get_species_index(self, str species):
        if species in self.species_dict:
            return self.species_dict[species]
        warnings.warn('Species not found in experimental lineage: %s\n' % species)
        return -1

    def __setstate__(self, state):
        super().__setstate__(state[:len(state)-1])
        self.species_dict = state[len(state)-1]

    def __getstate__(self):
        return super().__getstate__() + (self.species_dict,)

    def truncate_lineage(self,double start_time, double end_time):
        cdef Schnitz sch, new_sch
        cdef dict sch_dict = {}
        cdef ExperimentalLineage new_lineage = ExperimentalLineage()
        cdef np.ndarray newtime, newvolume, newdata, indices_to_keep

        sch_dict[None] = None
        for i in range(self.size()):
            sch = self.get_schnitz(i)

            newtime = sch.get_time().copy()
            newvolume = sch.get_volume().copy()
            newdata = sch.get_data().copy()

            # if the final time of this lineage is before the starting time
            # or the first time is before the end time, then skip it
            if newtime[newtime.shape[0]-1] < start_time or newtime[0] > end_time:
                sch_dict[sch] = None
                continue

            indices_to_keep = (sch.get_time() >= start_time) & (sch.get_time() <= end_time)
            newtime = newtime[indices_to_keep]
            newvolume = newvolume[indices_to_keep]
            newdata = newdata[indices_to_keep]

            sch_dict[sch] = Schnitz(newtime, newdata, newvolume)

        for i in range(self.size()):
            sch = self.get_schnitz(i)
            if sch_dict[sch] is not None:
                new_lineage.add_schnitz(sch_dict[sch])
                sch_dict[sch].py_set_parent( sch_dict[sch.get_parent()] )
                sch_dict[sch].py_set_daughters( sch_dict[sch.get_daughter_1()] , sch_dict[sch.get_daughter_2()] )

        new_lineage.py_set_species_indices(self.species_dict.copy())

        return new_lineage










