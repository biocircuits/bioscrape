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
import cython
from bioscrape.sbmlutil import add_species, add_parameter, add_reaction, add_rule, create_sbml_model, import_sbml

from libc.math cimport log, sqrt, cos, round, exp, fabs

# Define static epsilon for handling rounding errors
from libcpp.limits cimport numeric_limits
cdef float float_epsilon = numeric_limits[float].epsilon()

##################################################                ####################################################
######################################              PROPENSITY TYPES                    ##############################
#################################################                     ################################################


cdef class Propensity:
    """
    Base class for reaction propensity (rate law) functions.

    Subclasses implement `get_propensity` and `get_volume_propensity`
    (both `cdef`, called from Python via `py_get_propensity` and
    `py_get_volume_propensity`). By default, volume propensities fall
    back to the non-volume propensity, and stochastic propensities
    fall back to the deterministic ones; a subclass overrides either
    where the behavior actually differs (e.g.
    `BimolecularPropensity` and `MassActionPropensity`, whose
    stochastic propensities account for discrete copy-number effects
    when a species reacts with itself).
    """
    def __init__(self):
        """See class docstring."""
        self.propensity_type = PropensityType.unset

    def py_get_propensity(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                          double time = 0.0):
        """
        Compute the propensity at a given state and parameter vector.

        Parameters
        ----------
        state : numpy.ndarray
            The state vector.
        params : numpy.ndarray
            The parameter vector.
        time : float, optional
            The current time (default 0.0).

        Returns
        -------
        float
            The computed propensity (non-negative).
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
        #By default, volume propensitiesa are the same as regular propensities, unless otherwise noted
        return self.get_propensity(state, params, time)


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
        """
        Compute the volume-dependent propensity.

        Parameters
        ----------
        state : numpy.ndarray
            The state vector.
        params : numpy.ndarray
            The parameter vector.
        volume : float
            The cell volume.
        time : float, optional
            The current time (default 0.0).

        Returns
        -------
        float
            The computed propensity (non-negative).
        """
        return self.get_volume_propensity(<double*> state.data, <double*> params.data, volume, time)



    def initialize(self, dict param_dictionary, dict species_indices, dict parameter_indices):
        """
        Bind this propensity's field names to state/parameter indices.

        Called once by `~bioscrape.types.Model` when the model is
        initialized; not normally called directly. The recognized keys
        in `param_dictionary` (e.g. 'k', 's1', 'K', 'n') depend on the
        propensity subclass -- see each subclass's docstring.

        Parameters
        ----------
        param_dictionary : dict
            Maps field names (e.g. 'k', 's1') to species or parameter
            names.
        species_indices : dict
            Maps species names to their index in the state vector.
        parameter_indices : dict
            Maps parameter names to their index in the parameter
            vector.
        """
        pass

    def get_species_and_parameters(self, dict fields, **keywords):
        """
        The species and parameter names referenced by this propensity.

        Parameters
        ----------
        fields : dict
            The propensity's field dictionary (as passed to
            `initialize`).
        **keywords
            Ignored here; `~Model` calls every propensity's
            `get_species_and_parameters` uniformly with
            `species2index`/`params2index` (dicts mapping species/
            parameter names to their state/parameter vector index),
            for the benefit of subclasses (e.g. `GeneralPropensity`)
            that need to resolve a symbolic expression instead of
            just reading names out of `fields`.

        Returns
        -------
        tuple of list of str
            `(species_names, parameter_names)`.
        """
        return (None,None)



cdef class ConstitutivePropensity(Propensity):
    r"""
    Constant (zeroth-order) propensity, rate = k.

    Selected automatically by `~Model.create_reaction` for
    ``propensity_type='massaction'`` reactions with zero reactant
    species (the fastest of the specialized mass-action propensities);
    not normally constructed directly. Recognized
    `propensity_param_dict` key: 'k' (the rate parameter name).

    Attributes
    ----------
    rate_index : int
        The parameter index of the rate k.
    """
    # constructor
    def __init__(self):
        """See class docstring."""
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

    def __str__(self):
        return " {p" + str(self.rate_index) + "} "


cdef class UnimolecularPropensity(Propensity):
    r"""
    First-order propensity, rate = k*x.

    Selected automatically by `~Model.create_reaction` for
    ``propensity_type='massaction'`` reactions with one reactant
    species; not normally constructed directly. Recognized
    `propensity_param_dict` keys: 'k' (the rate parameter name) and
    'species' (the reactant species name x).

    Attributes
    ----------
    rate_index : int
        The parameter index of the rate k.
    species_index : int
        The species index of x.
    """
    # constructor
    def __init__(self):
        """See class docstring."""
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

    def __str__(self):
        return " {p" + str(self.rate_index) + "} * {" + str(self.species_index) + "} "



cdef class BimolecularPropensity(Propensity):
    r"""
    Second-order propensity, rate = k*x1*x2.

    Selected automatically by `~Model.create_reaction` for
    ``propensity_type='massaction'`` reactions with two reactant
    species; not normally constructed directly. Recognized
    `propensity_param_dict` keys: 'k' (the rate parameter name) and
    'species' (the two reactant species names, e.g. 'A*B'). If the
    two species are the same (e.g. 'A*A'), the stochastic
    propensity accounts for the discrete copy-number effect of a
    species reacting with itself: rate = k*x1*(x1 - 1).

    Attributes
    ----------
    rate_index : int
        The parameter index of the rate k.
    s1_index : int
        The species index of x1.
    s2_index : int
        The species index of x2.
    """
    # constructor
    def __init__(self):
        """See class docstring."""
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

    def __str__(self):
        return " {p" + str(self.rate_index) + "} * {s" + str(self.s1_index) + "} * {" + str(self.s2_index) + "} "


cdef class PositiveHillPropensity(Propensity):
    r"""
    Activating Hill propensity, rate = k*(x1/K)^n / (1 + (x1/K)^n).

    Constructed via `~Model.create_reaction` with
    ``propensity_type='hillpositive'``. Recognized
    `propensity_param_dict` keys: 'k' (rate), 'K' (Hill constant),
    'n' (cooperativity/Hill coefficient), and 's1' (the species x1).

    Attributes
    ----------
    K_index : int
        The parameter index of the Hill constant K.
    rate_index : int
        The parameter index of the rate k.
    n_index : int
        The parameter index of the cooperativity n.
    s1_index : int
        The species index of x1.
    """
    # constructor
    def __init__(self):
        """See class docstring."""
        self.propensity_type = PropensityType.hill_positive

    cdef double get_propensity(self, double* state, double* params, double time):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        if X < 0 and abs(X) < float_epsilon:
            X = 0.        # Fix rounding errors
        return rate * (X / K) ** n / (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        if X < 0 and abs(X) < float_epsilon:
            X = 0.        # Fix rounding errors
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

    def __str__(self):
        return " {p" + str(self.rate_index) + "} * ({s" + str(self.s1_index) + "} / {p" + str(self.K_index) + \
                "}) ^ {p" + str(self.n_index) + "} / (1 + ({s" + str(self.s1_index) + "}/{p" + str(self.K_index) + \
                "})^{p" + str(self.n_index) + "})"


cdef class PositiveProportionalHillPropensity(Propensity):
    r"""
    Activating Hill propensity proportional to a second species.

    rate = k*d*(x1/K)^n / (1 + (x1/K)^n). Constructed via
    `~Model.create_reaction` with
    ``propensity_type='proportionalhillpositive'``. Recognized
    `propensity_param_dict` keys: 'k' (rate), 'K' (Hill constant),
    'n' (cooperativity/Hill coefficient), 's1' (the species x1),
    and 'd' (the proportional species d, e.g. an enzyme or
    polymerase whose abundance scales the rate).

    Attributes
    ----------
    K_index : int
        The parameter index of the Hill constant K.
    rate_index : int
        The parameter index of the rate k.
    n_index : int
        The parameter index of the cooperativity n.
    s1_index : int
        The species index of x1.
    d_index : int
        The species index of d.
    """
    # constructor
    def __init__(self):
        """See class docstring."""
        self.propensity_type = PropensityType.proportional_hill_positive

    cdef double get_propensity(self, double* state, double* params, double time):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        cdef double d = state[self.d_index]
        if X < 0 and abs(X) < float_epsilon:
            X = 0.        # Fix rounding errors
        return rate * d *  (X / K) ** n / (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double d = state[self.d_index]
        cdef double rate = params[self.rate_index]
        if X < 0 and abs(X) < float_epsilon:
            X = 0.        # Fix rounding errors
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


    def __str__(self):
        return " {s" + str(self.d_index) + "} * {p" + str(self.rate_index) + "} * ({s" + str(self.s1_index) + "} / {p" + str(self.K_index) + \
                "}) ^ {p" + str(self.n_index) + "} / (1 + ({s" + str(self.s1_index) + "}/{p" + str(self.K_index) + \
                "})^{p" + str(self.n_index) + "})"

cdef class NegativeHillPropensity(Propensity):
    r"""
    Repressing Hill propensity, rate = k / (1 + (x1/K)^n).

    Constructed via `~Model.create_reaction` with
    ``propensity_type='hillnegative'``. Recognized
    `propensity_param_dict` keys: 'k' (rate), 'K' (Hill constant),
    'n' (cooperativity/Hill coefficient), and 's1' (the repressor
    species x1).

    Attributes
    ----------
    K_index : int
        The parameter index of the Hill constant K.
    rate_index : int
        The parameter index of the rate k.
    n_index : int
        The parameter index of the cooperativity n.
    s1_index : int
        The species index of x1.
    """
    # constructor
    def __init__(self):
        """See class docstring."""
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

    def __str__(self):
        return " {p" + str(self.rate_index) + "} / (1 + ({s" + str(self.s1_index) + "}/{p" + str(self.K_index) + \
                "})^{p" + str(self.n_index) + "})"


cdef class NegativeProportionalHillPropensity(Propensity):
    r"""
    Repressing Hill propensity proportional to a second species.

    rate = k*d / (1 + (x1/K)^n). Constructed via
    `~Model.create_reaction` with
    ``propensity_type='proportionalhillnegative'``. Recognized
    `propensity_param_dict` keys: 'k' (rate), 'K' (Hill constant),
    'n' (cooperativity/Hill coefficient), 's1' (the repressor species
    x1), and 'd' (the proportional species d).

    Attributes
    ----------
    K_index : int
        The parameter index of the Hill constant K.
    rate_index : int
        The parameter index of the rate k.
    n_index : int
        The parameter index of the cooperativity n.
    s1_index : int
        The species index of x1.
    d_index : int
        The species index of d.
    """
    # constructor
    def __init__(self):
        """See class docstring."""
        self.propensity_type = PropensityType.proportional_hill_negative

    cdef double get_propensity(self, double* state, double* params, double time):
        cdef double X = state[self.s1_index]
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double rate = params[self.rate_index]
        cdef double d = state[self.d_index]
        if X < 0 and abs(X) < float_epsilon:
            X = 0.        # Fix rounding errors
        return rate * d *  1/ (1 + (X/K)**n)

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        cdef double X = state[self.s1_index] / volume
        cdef double K = params[self.K_index]
        cdef double n = params[self.n_index]
        cdef double d = state[self.d_index]
        cdef double rate = params[self.rate_index]
        if X < 0 and abs(X) < float_epsilon:
            X = 0.        # Fix rounding errors
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
        """
        Bind the `s1` and `d` species names to state-vector indices.

        An alternative to `~Propensity.initialize` for setting just
        the species fields.

        Parameters
        ----------
        species : dict
            Maps 's1' and/or 'd' to species names.
        species_indices : dict
            Maps species names to their index in the state vector.
        """
        for key in species:
            if key == 's1':
                self.s1_index = species_indices[species['s1']]
            elif key == 'd':
                self.d_index = species_indices[species['d']]
            else:
                logging.info('Warning! Useless field for NegativeProportionalHillPropensity '+str(key))
    def set_parameters(self, parameters, parameter_indices):
        """
        Bind the `K`, `n`, and `k` parameter names to indices.

        An alternative to `~Propensity.initialize` for setting just
        the parameter fields.

        Parameters
        ----------
        parameters : dict
            Maps 'K', 'n', and/or 'k' to parameter names.
        parameter_indices : dict
            Maps parameter names to their index in the parameter
            vector.
        """
        for key in parameters:
            if key == 'K':
                self.K_index = parameter_indices[parameters[key]]
            elif key == 'n':
                self.n_index = parameter_indices[parameters[key]]
            elif key == 'k':
                self.rate_index = parameter_indices[parameters[key]]
            else:
                logging.info('Warning! Useless field for NegativeProportionalHillPropensity '+str(key))

    def __str__(self):
        return " {s" + str(self.d_index) + "} * {p" + str(self.rate_index) + "} / (1 + ({s" + str(self.s1_index) + "}/{p" + str(self.K_index) + \
                "})^{p" + str(self.n_index) + "})"


cdef class MassActionPropensity(Propensity):
    r"""
    General mass-action propensity, rate = k times the reactant product.

    Selected automatically by `~Model.create_reaction` for
    ``propensity_type='massaction'`` reactions with three or more
    reactant species (`ConstitutivePropensity`,
    `UnimolecularPropensity`, and `BimolecularPropensity` are used
    instead for zero, one, and two reactant species respectively);
    not normally constructed directly. Recognized
    `propensity_param_dict` keys: 'k' (the rate parameter name) and
    'species' (the reactant species names, e.g. 'A*A*B' for a
    reaction where species 'A' appears with multiplicity 2). The
    deterministic propensity multiplies each distinct species once;
    the stochastic propensity additionally accounts for repeated
    species via falling-factorial terms (e.g. x*(x-1) for a species
    with multiplicity 2), matching combinatorial reaction kinetics.

    Attributes
    ----------
    sp_inds : list of int
        The species index of each distinct reactant species.
    sp_counts : list of int
        The multiplicity of each species in `sp_inds`.
    k_index : int
        The parameter index of the rate k.
    """
    def __init__(self):
        """See class docstring."""
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

    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        return self.get_propensity(state, params, time) / (volume ** (self.num_species - 1) )

    cdef double get_stochastic_volume_propensity(self, double *state, double *params, double volume, double time):
        return self.get_stochastic_propensity(state, params, time) / (volume ** (self.num_species - 1))


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

    def __str__(self):
        string_rep = " {p" + str(self.k_index) + "} "
        for i in range(len(self.sp_inds)):
            string_rep += "* ({s" + str(self.sp_inds[i]) + "}^" + str(self.sp_counts) + ") "
        return string_rep


##################################################                ####################################################
######################################              PARSING                             ##############################
#################################################                     ################################################

@cython.auto_pickle(True)
cdef class Term:
    """
    Base class for nodes in a parsed symbolic expression tree.

    Used internally by `GeneralPropensity` and
    `~bioscrape.types.GeneralODERule` to evaluate arbitrary symbolic rate/rule
    expressions (parsed from a string); not normally constructed directly.
    Subclasses implement `evaluate` and `volume_evaluate` (both `cdef`),
    called from Python via `py_evaluate`/`py_volume_evaluate`.
    """
    cdef double evaluate(self, double *species, double *params, double time):
        raise SyntaxError('Cannot make Term base object')

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        raise SyntaxError('Cannot make Term base object')


    def py_evaluate(self, np.ndarray species, np.ndarray params, double time=0.0):
        """
        Evaluate this expression node.

        Parameters
        ----------
        species : numpy.ndarray
            The state (species) vector.
        params : numpy.ndarray
            The parameter vector.
        time : float, optional
            The current time (default 0.0).

        Returns
        -------
        float
            The value of this expression node.
        """
        return self.evaluate(<double*> species.data, <double*> params.data, time)

    def py_volume_evaluate(self, np.ndarray species, np.ndarray params,
                           double vol, double time=0.0):
        """
        Evaluate this expression node with a cell volume.

        Parameters
        ----------
        species : numpy.ndarray
            The state (species) vector.
        params : numpy.ndarray
            The parameter vector.
        vol : float
            The cell volume.
        time : float, optional
            The current time (default 0.0).

        Returns
        -------
        float
            The value of this expression node.
        """
        return self.volume_evaluate(<double*> species.data, <double*> params.data,
                                    vol, time)

# Base building blocks
@cython.auto_pickle(True)
cdef class ConstantTerm(Term):
    """
    A constant-valued expression node.

    Parameters
    ----------
    val : float
        The constant value.
    """

    def __init__(self, double val):
        """See class docstring."""
        self.value = val

    cdef double evaluate(self, double *species, double *params, double time):
        return self.value
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return self.value

    def __str__(self):
        return f" {self.value} "

@cython.auto_pickle(True)
cdef class SpeciesTerm(Term):
    """
    An expression node referencing a species value by index.

    Parameters
    ----------
    ind : int
        The species index.
    """
    def __init__(self, unsigned ind):
        """See class docstring."""
        self.index = ind

    cdef double evaluate(self, double *species, double *params, double time):
        return species[self.index]
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return species[self.index]

    def __str__(self):
        return " {s" + str(self.index) + "} "

@cython.auto_pickle(True)
cdef class ParameterTerm(Term):
    """
    An expression node referencing a parameter value by index.

    Parameters
    ----------
    ind : int
        The parameter index.
    """
    def __init__(self, unsigned ind):
        """See class docstring."""
        self.index = ind

    cdef double evaluate(self, double *species, double *params, double time):
        return params[self.index]
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return params[self.index]

    def __str__(self):
        return " {p" + str(self.index) + "} "

cdef class VolumeTerm(Term):
    """
    An expression node representing the cell volume.

    Evaluates to 1.0 via `evaluate` (no volume given) and to the cell
    volume via `volume_evaluate`.
    """
    cdef double evaluate(self, double *species, double *params, double time):
        return 1.0
    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return vol

    def __str__(self):
        return " {volume} "

# Putting stuff together

#   To make BinaryTerms picklable, we need to use a __reduce__ 
#   function, which needs to return a function that can reconstruct a BinaryTerm
#   from a state. Internally-defined functions can't be pickled, so that 
#   function has to be external to the class. Thus, we need this function.
# 
#   I also tried to make BinaryTerms picklable using __getstate__ and 
#   __setstate__, which would obviate the need for this ugly external function,
#   but that resulted in mysterious errors about calling unsafe __new__() 
#   methods, which I couldn't resolve. Sorry. =(
def restore_binary_term(state, ClassName):
    new_term = ClassName()
    if state is not None:
        for i, x in enumerate(state):
            new_term.py_add_term(x)
    return new_term

cdef class BinaryTerm(Term):
    """
    Base class for expression nodes that combine a list of child terms.

    Subclasses (`SumTerm`, `ProductTerm`, `MaxTerm`, `MinTerm`)
    implement `evaluate`/`volume_evaluate` by reducing over the child
    terms added via `py_add_term`.
    """
    def __init__(self):
        """See class docstring."""
        self.terms_list = []

    def __reduce__(self):
        return (restore_binary_term, (self.terms_list, self.__class__))

    cdef void add_term(self,Term trm):
        self.terms.push_back(<void*> trm)
        self.terms_list.append(trm)

    def py_add_term(self, Term trm):
        """
        Append a child expression node.

        Parameters
        ----------
        trm : Term
            The child term to add.
        """
        self.add_term(trm)

cdef class SumTerm(BinaryTerm):
    """An expression node summing its child terms."""
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

    def __str__(self):
        string_rep = " (" + str(self.terms_list[0])
        for i in range(1, len(self.terms_list)):
            string_rep += " + " + str(self.terms_list[i])
        string_rep += ") "
        return string_rep

cdef class ProductTerm(BinaryTerm):
    """An expression node multiplying its child terms."""
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

    def __str__(self):
        string_rep = " (" + str(self.terms_list[0])
        for i in range(1, len(self.terms_list)):
            string_rep += " * " + str(self.terms_list[i])
        string_rep += ") "
        return string_rep

cdef class MaxTerm(BinaryTerm):
    """An expression node taking the maximum of its child terms."""
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

    def __str__(self):
        string_rep = " max(" + str(self.terms_list[0])
        for i in range(1, len(self.terms_list)):
            string_rep += ", " + str(self.terms_list[i])
        string_rep += ") "
        return string_rep

cdef class MinTerm(BinaryTerm):
    """An expression node taking the minimum of its child terms."""
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


    def __str__(self):
        string_rep = " min(" + str(self.terms_list[0])
        for i in range(1, len(self.terms_list)):
            string_rep += ", " + str(self.terms_list[i])
        string_rep += ") "
        return string_rep

@cython.auto_pickle(True)
cdef class PowerTerm(Term):
    """An expression node raising a base term to an exponent term."""
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

    def __str__(self):
        return " (" + str(self.base) + " ^ " + str(self.exponent) + ") "

@cython.auto_pickle(True)
cdef class ExpTerm(Term):
    """An expression node computing the natural exponential of its argument term."""
    cdef void set_arg(self, Term arg):
        self.arg = arg

    cdef double evaluate(self, double *species, double *params, double time):
        return exp(self.arg.evaluate(species,params,time))

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return exp(self.arg.volume_evaluate(species,params,vol,time))

    def __str__(self):
        return " (e ^ " + str(self.exponent) + ") "

@cython.auto_pickle(True)
cdef class LogTerm(Term):
    """An expression node computing the natural log of its argument."""
    cdef void set_arg(self, Term arg):
        self.arg = arg

    cdef double evaluate(self, double *species, double *params, double time):
        return log(self.arg.evaluate(species,params,time))

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return log(self.arg.volume_evaluate(species,params,vol,time))

    def __str__(self):
        return " log(" + str(self.arg) + ") "

@cython.auto_pickle(True)
cdef class StepTerm(Term):
    """A Heaviside step node: 1 if its argument is >= 0, else 0."""
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

    def __str__(self):
        return " Heaviside(" + str(self.arg) + ") "

@cython.auto_pickle(True)
cdef class AbsTerm(Term):
    """An expression node computing the absolute value of its argument."""
    cdef void set_arg(self, Term arg):
        self.arg = arg

    cdef double evaluate(self, double *species, double *params, double time):
        return fabs( self.arg.evaluate(species,params,time) )

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return fabs( self.arg.volume_evaluate(species,params,vol,time) )

    def __str__(self):
        return " abs(" + str(self.arg) + ") "

@cython.auto_pickle(True)
cdef class TimeTerm(Term):
    """An expression node representing the current simulation time."""
    cdef double evaluate(self, double *species, double *params, double time):
        return time

    cdef double volume_evaluate(self, double *species, double *params, double vol, double time):
        return time

    def __str__(self):
        return " {time} "

def sympy_species_and_parameters(instring, species2index = None, params2index = None):
    instring = instring.replace('^','**')
    instring = instring.replace('|','_')
    root = sympy.sympify(instring, _clash1)
    nodes = [root]
    index = 0
    while index < len(nodes):
        node = nodes[index]
        index += 1
        nodes.extend(node.args)

    names = [str(n) for n in nodes if type(n) == sympy.Symbol]\
            +[str(n)[1:] for n in nodes if type(n) == sympy.Symbol if str(n)[0] == "_"]

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
    """
    Propensity defined by an arbitrary symbolic rate expression.

    Constructed via `~Model.create_reaction` with
    ``propensity_type='general'``. Recognized
    `propensity_param_dict` key: 'rate', a string expression in
    species and parameter names (e.g. 'k*A*B'), parsed into a
    `Term` expression tree via `parse_expression`.
    """

    cdef double get_propensity(self, double* state, double* params, double time):
        return self.term.evaluate(state,params,time)
    cdef double get_volume_propensity(self, double *state, double *params, double volume, double time):
        return self.term.volume_evaluate(state,params,volume,time)

    def __init__(self):
        """See class docstring."""
        self.propensity_type = PropensityType.general

    def initialize(self, dict dictionary, dict species2index, dict params2index):
        """
        Parse the `rate` field into an expression tree.

        Parameters
        ----------
        dictionary : dict
            Must contain a `rate` key with the symbolic rate
            expression string.
        species2index : dict
            Maps species names to their index in the state vector.
        params2index : dict
            Maps parameter names to their index in the parameter
            vector.
        """
        instring = dictionary['rate']
        self.term = parse_expression(instring, species2index, params2index)

    def get_species_and_parameters(self, dict fields, dict species2index, dict params2index):
        return sympy_species_and_parameters(fields['rate'].strip(), species2index, params2index)

    def py_get_term(self):
        """
        The parsed expression tree for this propensity's rate.

        Returns
        -------
        Term
            The root node of the parsed `rate` expression.
        """
        return self.term

    def __str__(self):
        return str(self.term)


##################################################                ####################################################
######################################              DELAY TYPES                        ##############################
#################################################                     ################################################

cdef class Delay:
    """
    Base class for stochastic reaction delay distributions.

    Contains one function, `get_delay` (`cdef`, called from Python via
    `py_get_delay`), that returns a delay as a function of the state
    and parameter vectors; must be overridden by subclasses. Attached
    to a reaction via `~Model.create_reaction`'s `delay_type` argument
    (e.g. 'fixed', 'gaussian', 'gamma') rather than constructed
    directly.

    Attributes
    ----------
    delay_type : DelayType
        The type of delay.
    """
    def __init__(self):
        """See class docstring."""
        self.delay_type = DelayType.unset_delay

    def py_get_delay(self, np.ndarray[np.double_t,ndim=1] state,
                     np.ndarray[np.double_t,ndim=1] params):
        """
        Sample/compute the delay at a given state and parameter vector.

        This function should NOT be overridden by subclasses. It is
        just a Python wrapper of the `cdef` `get_delay` function,
        which subclasses should override instead.

        Parameters
        ----------
        state : numpy.ndarray
            The state vector.
        params : numpy.ndarray
            The parameter vector.

        Returns
        -------
        float
            The computed delay.
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
        Bind this delay's field names to state/parameter indices.

        Called once by `~bioscrape.types.Model` when the model is
        initialized; not normally called directly. The recognized
        keys in `param_dictionary` (e.g. 'delay', 'mean', 'std', 'k',
        'theta') depend on the delay subclass -- see each subclass's
        docstring.

        Parameters
        ----------
        param_dictionary : dict
            Maps field names to species or parameter names.
        species_indices : dict
            Maps species names to their index in the state vector.
        parameter_indices : dict
            Maps parameter names to their index in the parameter
            vector.
        """
        pass

    def get_species_and_parameters(self, dict fields, **keywords):
        """
        The species and parameter names referenced by this delay.

        Parameters
        ----------
        fields : dict
            The delay's field dictionary (as passed to `initialize`).
        **keywords
            Ignored here; `~Model` calls every delay's
            `get_species_and_parameters` uniformly with
            `species2index`/`params2index` (dicts mapping species/
            parameter names to their state/parameter vector index),
            not used by any built-in `Delay` subclass but available
            to custom subclasses that need them.

        Returns
        -------
        tuple of list of str
            `(species_names, parameter_names)`.
        """
        return [],[]


cdef class NoDelay(Delay):
    """A delay of 0.0 (no delay); the default when no delay is specified."""
    def __init__(self):
        """See class docstring."""
        self.delay_type = DelayType.none

    cdef double get_delay(self, double* state, double* params):
        return 0.0

cdef class FixedDelay(Delay):
    """
    A fixed (deterministic) delay.

    Constructed via `~Model.create_reaction` with
    ``delay_type='fixed'``. Recognized `delay_param_dict` key: 'delay'
    (the name of the parameter giving the fixed delay value).

    Attributes
    ----------
    delay_index : int
        The parameter index of the delay value.
    """

    def __init__(self):
        """See class docstring."""
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
    """
    A Gaussian-distributed delay, resampled at each reaction firing.

    Constructed via `~Model.create_reaction` with
    ``delay_type='gaussian'``. Recognized `delay_param_dict` keys:
    'mean' and 'std' (parameter names giving the distribution's mean
    and standard deviation).

    Attributes
    ----------
    mean_index : int
        The parameter index of the mean.
    std_index : int
        The parameter index of the standard deviation.
    """

    def __init__(self):
        """See class docstring."""
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
    r"""
    A gamma-distributed delay with shape k and scale theta.

    Resampled at each reaction firing; equivalent to the sum of k
    independent exponential random variables with mean theta; has
    mean k*theta and variance k*theta^2.

    Constructed via `~Model.create_reaction` with
    ``delay_type='gamma'``. Recognized `delay_param_dict` keys: 'k'
    and 'theta' (parameter names giving the shape and scale).

    Attributes
    ----------
    k_index : int
        The parameter index of the shape k.
    theta_index : int
        The parameter index of the scale theta.
    """

    def __init__(self):
        """See class docstring."""
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
    Base class for model rules.

    Rules apply updates outside the normal reaction network, either
    once at the start of a simulation or repeatedly during it.
    Subclasses implement `rule_operation` and `rule_volume_operation`
    (both `cdef`); `execute_rule`/`execute_volume_rule` (also `cdef`,
    called from Python via `py_execute_rule`/`py_execute_volume_rule`)
    apply them according to the rule's frequency, set via
    `set_frequency_flag`.
    """
    cdef void rule_operation(self, double *state, double *params, double time, double dt):
        raise NotImplementedError('Creating base Rule class. This should be subclassed.')

    cdef void rule_volume_operation(self, double *state, double *params, double volume, double time, double dt):
        self.rule_operation(state, params, time, dt)

    cdef void execute_rule(self, double *state, double *params, double time, double dt, unsigned rule_step):
        if self.frequency_flag == -1 or self.frequency_flag == time or (rule_step and self.frequency_flag == -2):
            self.rule_operation(state, params, time, dt)

    cdef void execute_volume_rule(self, double *state, double *params, double volume, double time, double dt, unsigned rule_step):
        if self.frequency_flag == -1 or self.frequency_flag == time or (rule_step and self.frequency_flag == -2):
            self.rule_volume_operation(state, params, volume, time, dt)

    def py_execute_rule(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                        double time = 0.0, double dt = .01, rule_step = True):
        """
        Apply this rule if it is due to fire.

        Parameters
        ----------
        state : numpy.ndarray
            The state vector, updated in place if the rule fires.
        params : numpy.ndarray
            The parameter vector, updated in place if the rule targets
            a parameter.
        time : float, optional
            The current time (default 0.0).
        dt : float, optional
            The simulation time step, used by 'dt'-frequency rules
            (default 0.01).
        rule_step : bool, optional
            Whether this call corresponds to a full simulation time
            step; relevant for 'dt'-frequency rules (default True).
        """
        self.execute_rule(<double*> state.data, <double*> params.data,time, dt, rule_step)

    def py_execute_volume_rule(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                               double volume, double time=0.0, double dt = .01, rule_step = True):
        """
        Apply this rule (with a cell volume) if it is due to fire.

        Parameters
        ----------
        state : numpy.ndarray
            The state vector, updated in place if the rule fires.
        params : numpy.ndarray
            The parameter vector, updated in place if the rule targets
            a parameter.
        volume : float
            The cell volume.
        time : float, optional
            The current time (default 0.0).
        dt : float, optional
            The simulation time step, used by 'dt'-frequency rules
            (default 0.01).
        rule_step : bool, optional
            Whether this call corresponds to a full simulation time
            step; relevant for 'dt'-frequency rules (default True).
        """
        self.execute_volume_rule(<double*> state.data, <double*> params.data, volume,time, dt, rule_step)

    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices, rule_frequency = "repeat"):
        """
        Bind this rule's field names to state/parameter indices.

        Called once by `~bioscrape.types.Model` when the model is
        initialized; not normally called directly. The recognized
        keys in `dictionary` depend on the rule subclass -- see each
        subclass's docstring.

        Parameters
        ----------
        dictionary : dict
            Maps field names to species/parameter names or
            expressions.
        species_indices : dict
            Maps species names to their index in the state vector.
        parameter_indices : dict
            Maps parameter names to their index in the parameter
            vector.
        rule_frequency : str or float, optional
            How often the rule fires; see `set_frequency_flag`
            (default 'repeat').
        """
        pass

    def get_species_and_parameters(self, dict fields, **keywords):
        """
        The species and parameter names referenced by this rule.

        Parameters
        ----------
        fields : dict
            The rule's field dictionary (as passed to `initialize`).

        Returns
        -------
        tuple of list of str
            `(species_names, parameter_names)`.
        """
        return (None,None)

    def set_frequency_flag(self, rule_frequency):
        """
        Set when this rule fires.

        Parameters
        ----------
        rule_frequency : str or float
            'start' fires the rule once, at simulation start. 'repeat'
            (or 'repeated') fires it at every reaction/timepoint. 'dt'
            fires it once per simulation time step (used by
            `GeneralODERule`). A non-negative number fires it only at
            that specific time.

        Raises
        ------
        ValueError
            If `rule_frequency` is not one of the recognized strings
            and cannot be interpreted as a non-negative float.
        """
        if rule_frequency == "start":
            self.frequency_flag = 0.0

        elif rule_frequency == "repeat" or rule_frequency == "repeated":
            self.frequency_flag = -1.0

        elif rule_frequency == "dt":
            self.frequency_flag = -2.0

        elif float(rule_frequency) >= 0:
            self.frequency_flag = float(rule_frequency)

        else:
            raise ValueError(f"Invalid rule frequency: {rule_frequency} for {self}.")


cdef class AdditiveAssignmentRule(Rule):
    """
    Assigns a species to the sum of a list of other species.

    Constructed via `~Model.create_rule` with
    ``rule_type='additive'``. `rule_attributes['equation']` must be of
    the form 'dest = src1 + src2 + ...', naming the destination
    species and the source species to sum (no other expression forms
    are supported; see `GeneralAssignmentRule` for arbitrary
    expressions).
    """

    cdef void rule_operation(self, double *state, double *params, double time, double dt):
        cdef unsigned i = 0
        cdef double answer = 0.0
        for i in range(self.species_source_indices.size()):
            answer += state[ self.species_source_indices[i] ]

        state[self.dest_index] = answer

    def initialize(self, dict dictionary, dict species_indices, dict parameter_indices, rule_frequency = "repeat"):
        self.set_frequency_flag(rule_frequency)
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
    Assigns a species or parameter to an arbitrary expression's value.

    Constructed via `~Model.create_rule` with
    ``rule_type='assignment'``. `rule_attributes['equation']` must be
    of the form 'dest = expression', where `expression` is parsed
    (via `parse_expression`) into a `Term` tree and re-evaluated each
    time the rule fires; `dest` may be a species or a parameter name.
    """
    cdef void rule_operation(self, double *state, double *params, double time, double dt):
        if self.param_flag > 0:
            params[self.dest_index] = self.rhs.evaluate(state,params,time)
        else:
            state[self.dest_index] = self.rhs.evaluate(state,params,time)

    cdef void rule_volume_operation(self, double *state, double *params, double volume, double time, double dt):
        if self.param_flag > 0:
            params[self.dest_index] = self.rhs.volume_evaluate(state,params,volume, time)
        else:
            state[self.dest_index] = self.rhs.volume_evaluate(state,params,volume, time)

    def initialize(self, dict fields, species2index, params2index, rule_frequency = "repeat"):
        self.set_frequency_flag(rule_frequency)
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


cdef class GeneralODERule(Rule):
    """
    Numerically integrates a target via Euler's method.

    Constructed via `~Model.create_rule` with ``rule_type='ode'``
    (frequency is forced to fire once per simulation time step).
    `rule_attributes` must contain 'target' (the destination species
    or parameter name) and 'equation' (a symbolic expression `f`,
    parsed via `parse_expression`); each firing adds
    f(state, params, time) * dt to the target's current value (a
    first-order Euler step).
    """
    cdef void rule_operation(self, double *state, double *params, double time, double dt):
        if self.param_flag > 0:
            params[self.dest_index] = params[self.dest_index] + self.rhs.evaluate(state,params,time)*dt
        else:
            state[self.dest_index] = state[self.dest_index] + self.rhs.evaluate(state,params,time)*dt

    cdef void rule_volume_operation(self, double *state, double *params, double volume, double time, double dt):
        if self.param_flag > 0:
            params[self.dest_index] = params[self.dest_index] + self.rhs.volume_evaluate(state,params,volume, time)*dt
        else:
            state[self.dest_index] = state[self.dest_index] + self.rhs.volume_evaluate(state,params,volume, time)*dt

    def initialize(self, dict fields, species2index, params2index, rule_frequency = "dt"):
        print("Initializing ODE Rule")
        self.set_frequency_flag(rule_frequency)
        self.rhs = parse_expression(fields['equation'], species2index, params2index)

        dest_name = fields['target'].strip()

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
        dest_name = fields['target'].strip()
        instring = fields['equation']

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
    """
    Base class for cell volume models.

    Used by volume-aware simulators
    (`~bioscrape.simulator.VolumeSSASimulator`,
    `~bioscrape.simulator.DelayVolumeSSASimulator`) to track a cell's
    volume and determine when it divides. Subclasses implement
    `get_volume_step`, `initialize`, `cell_divided`, and `copy` (all
    `cdef`).

    Attributes
    ----------
    current_volume : float
        The current volume of the cell (should be positive).
    """
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
        A deep copy of this volume object.

        Returns
        -------
        Volume
            The copy.
        """
        return self.copy()

    def py_get_volume_step(self, np.ndarray[np.double_t,ndim=1] state, np.ndarray[np.double_t,ndim=1] params,
                           double time, double volume, double dt):
        """
        Compute the volume change over a time step.

        Parameters
        ----------
        state : numpy.ndarray
            The state vector.
        params : numpy.ndarray
            The parameter vector.
        time : float
            The ending time after the volume step has occurred.
        volume : float
            The volume before the volume step occurs.
        dt : float
            The width of the time step.

        Returns
        -------
        float
            The change in cell volume from `time - dt` to `time`.
        """
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
        """
        Initialize the volume at the start of a simulation.

        Required to handle non-memoryless properties (e.g. a
        pre-sampled division time); called once before simulation
        begins.

        Parameters
        ----------
        state : numpy.ndarray
            The state vector.
        params : numpy.ndarray
            The parameter vector.
        time : float
            The initial time.
        volume : float
            The initial volume.
        """
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
        """
        Whether the cell divided during a time step.

        Parameters
        ----------
        state : numpy.ndarray
            The state vector.
        params : numpy.ndarray
            The parameter vector.
        time : float
            The ending time of the time step.
        volume : float
            The volume after the time step occurred.
        dt : float
            The width of the time step.

        Returns
        -------
        bool
            True if the cell divided in `[time - dt, time]`.
        """
        return self.cell_divided(<double*> state.data, <double*> params.data, time, volume, dt)


    def py_set_volume(self, double v):
        """
        Set the current volume.

        Parameters
        ----------
        v : float
            The volume to set (should be positive).
        """
        self.set_volume(v)
    def py_get_volume(self):
        """
        The current volume.

        Returns
        -------
        float
            The current volume.
        """
        return self.get_volume()




cdef class StochasticTimeThresholdVolume(Volume):
    r"""
    A cell that grows exponentially and divides at a random time.

    Growth is deterministic; the division time is sampled once, at
    `initialize`, as T times a normally-distributed random variable
    with mean 1 and standard deviation `division_noise`, where T is
    the deterministic time remaining to reach
    `average_division_volume` at the current growth rate.

    Parameters
    ----------
    cell_cycle_time : float
        The average cell cycle time; sets the growth rate to
        log(2) / cell_cycle_time.
    average_division_volume : float
        The average volume at which the cell divides.
    division_noise : float
        The relative noise (coefficient of variation, much less than 1)
        in the division time.

    Attributes
    ----------
    division_time : float
        The (pre-sampled) time at which the cell will divide.
    growth_rate : float
        The volume growth rate g, where dV/dt = g*V.
    """
    def __init__(self, double cell_cycle_time, double average_division_volume, double division_noise):
        """See class docstring."""
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
    r"""
    A cell with a state-dependent growth rate and a random division volume.

    Must be configured via `setup` after construction (`__init__`
    takes no arguments).

    Attributes
    ----------
    division_volume : float
        The volume at which the cell will divide.
    average_division_volume : float
        The average volume at which to divide.
    division_noise : float
        The relative noise (coefficient of variation, much less than 1)
        in the division volume.
    growth_rate : Term
        The growth rate expression, evaluated based on the state.
    """

    def __init__(self):
        """See class docstring."""
        pass

    def setup(self, double average_division_volume, double division_noise, growth_rate, Model m):
        """
        Configure the growth rate and division volume distribution.

        Parameters
        ----------
        average_division_volume : float
            The average volume at which to divide.
        division_noise : float
            The relative noise (coefficient of variation) in the
            division volume.
        growth_rate : str
            A symbolic expression for the growth rate as a function of
            model species/parameters, parsed via
            `~Model.parse_general_expression`.
        m : Model
            The model providing the species/parameter context used to
            parse `growth_rate`.
        """
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
    """A chemical reaction network model with optional delays.

    A `Model` can be built from an SBML file, a (deprecated) bioscrape
    XML file, or programmatically from the `species`, `reactions`,
    `parameters`, and `rules` arguments. These can be combined, and
    the model API (`create_reaction`, `create_rule`, `set_parameter`,
    etc.) can be used to extend the model further after construction.
    See the bioscrape wiki for full API documentation:
    https://github.com/biocircuits/bioscrape/wiki

    Volumes are not represented internally; they must be supplied
    externally, e.g. by a `Volume` object driving a simulator.

    Parameters
    ----------
    sbml_filename : str, optional
        Path to an SBML file to import. Cannot be combined with any
        other constructor argument except `input_printout` and
        `**kwargs`; import the SBML model first, then edit it with
        the API. If neither `sbml_filename` nor `filename` is given,
        no file is loaded and the model is built entirely from
        `species`, `reactions`, `parameters`, and `rules`.
    filename : str, optional
        Path to a (deprecated) bioscrape XML file to import. Prefer
        `sbml_filename` for new models.
    species : list of str, default []
        Species names to add to the model.
    reactions : list of tuple, default []
        Reaction tuples, each either
        `(reactants, products, propensity_type, propensity_param_dict)`
        or the 8-tuple form that also specifies a delay. See the
        bioscrape wiki for the full reaction tuple format.
    parameters : list of tuple or dict, default []
        `(name, value)` pairs, or a dict, giving initial parameter
        values.
    rules : list of tuple, default []
        Rule tuples, each `(rule_type, rule_attributes)` or
        `(rule_type, rule_attributes, rule_frequency)`. See the
        bioscrape wiki for the rule tuple format.
    initial_condition_dict : dict, optional
        Maps species names to initial values.
    input_printout : bool, default False
        If True, print verbose status messages while constructing
        the model.
    initialize_model : bool, default True
        If True, validate the model and build its internal arrays
        (see `py_initialize`) after construction.
    **kwargs
        Additional keyword arguments forwarded to `import_sbml` when
        loading from `sbml_filename`.

    Attributes
    ----------
    species2index : dict
        Maps species names to their index in `species_values` and in
        the model's stoichiometric matrices.
    params2index : dict
        Maps parameter names to their index in `params_values`.
    species_values : numpy.ndarray
        Current species values, indexed via `species2index`.
    params_values : numpy.ndarray
        Parameter values, indexed via `params2index`.
    propensities : list of Propensity
        Propensity object for each reaction, in the order the
        reactions were added.
    delays : list of Delay
        Delay object for each reaction, in the order the reactions
        were added.
    update_array : numpy.ndarray
        Stoichiometric matrix (num_species x num_reactions) for
        changes that occur immediately.
    delay_update_array : numpy.ndarray
        Stoichiometric matrix (num_species x num_reactions) for
        changes that occur after a delay.
    """
    def __init__(self, sbml_filename = None, filename = None, species = [], reactions = [],
                 parameters = [], rules = [], initial_condition_dict = None,
                 input_printout = False, initialize_model = True, **kwargs):
        """See class docstring."""

        ########################################################################
        # DEVELOPER WARNING 
        # 
        # In order to be copiable and usable with multiprocessing, Model must be 
        # picklable. To do that, Model implements a __getstate__ method and a 
        # __setstate__ method, which respectively compress all of the Model's 
        # state variables into a picklable tuple and use those tuples to make a 
        # new Model identical to the old one. 
        # 
        # IF YOU ADD, REMOVE, OR CHANGE ANY VARIABLES HERE, YOU MUST REFLECT 
        # THOSE CHANGES IN THE __getstate__ AND __setstate__ METHODS.
        #
        # This is especially important for newly-added variables. If you add 
        # variables but don't update the pickling methods, then you will 
        # introduce SILENT bugs whenever a user makes a copy of a Model or tries
        # to use a Model in multiple threads/processes with multiprocessing. 
        ########################################################################
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
        """Validate the model and build its internal C-level arrays.

        Creates the propensity/delay C-vectors and the stoichiometric
        matrices, and checks that all parameters and species have
        values (see `check_parameters` and `check_species`). Called
        automatically by the constructor unless
        `initialize_model=False` was passed.
        """
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
        if species not in self.species2index and species is not None and species != '':
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
            #Now anything not declared as a Species will be interpreted as a parameter
            pass
        for param_name in param_names:
            self._add_param(param_name)

        #Moved to Model._initialize
        self.delay_reaction_updates.append(delay_reaction_update_dict)
        delay_object.initialize(delay_param_dict, self.species2index, self.params2index)
        self.reaction_list.append((propensity_object, delay_object, reaction_update_dict, delay_reaction_update_dict))


    def create_propensity(self, propensity_type, propensity_param_dict, input_printout = False):
        """Construct a `Propensity` object of the given type.

        Automatically picks the fastest matching subclass for
        `propensity_type='massaction'`, based on the number of
        reactant species (`ConstitutivePropensity`,
        `UnimolecularPropensity`, `BimolecularPropensity`, or
        `MassActionPropensity`).

        Parameters
        ----------
        propensity_type : str
            One of 'massaction', 'hillpositive',
            'proportionalhillpositive', 'hillnegative',
            'proportionalhillnegative', or 'general'.
        propensity_param_dict : dict
            Parameters for the given propensity type. Numeric values
            are automatically converted to dummy parameters.
        input_printout : bool, default False
            If True, print the propensity type and parameters.

        Returns
        -------
        Propensity
            The constructed propensity object.

        Raises
        ------
        SyntaxError
            If `propensity_type` is not a recognized type.
        """
        if input_printout:
            print("Creating Propensity: prop_type="+str(propensity_type)+" params="+str(propensity_param_dict))
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

    def create_reaction(self, reactants, products, propensity_type,
                        propensity_param_dict, delay_type = None,
                        delay_reactants = None, delay_products = None,
                        delay_param_dict = None, input_printout = False):
        """Create a reaction and add it to the model.

        Supports all native propensity types and delay types.

        Parameters
        ----------
        reactants : list of str
            Reactant species names.
        products : list of str
            Product species names.
        propensity_type : str
            One of 'massaction', 'hillpositive',
            'proportionalhillpositive', 'hillnegative',
            'proportionalhillnegative', or 'general'.
        propensity_param_dict : dict
            Parameters for the given propensity type.
        delay_type : str, optional
            One of 'none', 'fixed', 'gaussian', or
            'gamma'. Defaults to no delay.
        delay_reactants : list of str, optional
            Reactant species names for the delayed part of the
            reaction.
        delay_products : list of str, optional
            Product species names for the delayed part of the
            reaction.
        delay_param_dict : dict, optional
            Parameters for the delay distribution.
        input_printout : bool, default False
            If True, print the reaction being created.
        """
        if input_printout:
            print("creating reaction with:"+
                "\n\tPropensity_type="+str(propensity_type)+" Inputs="+str(reactants)+" Outputs="+str(products)+
                "\n\tpropensity_param_dict="+str(propensity_param_dict)+
                "\n\tdelay type="+str(delay_type)+" delay inputs="+str(delay_reactants)+" delay outputs="+str(delay_products)+
                "\n\tdelay parameters="+str(delay_param_dict))
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

        prop_object = self.create_propensity(propensity_type, propensity_param_dict, input_printout = input_printout)

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
        self.reaction_definitions.append((reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict))



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

    def create_rule(self, rule_type, rule_attributes, rule_frequency = "repeated", input_printout = False):
        """Create a rule and add it to the model.

        Parameters
        ----------
        rule_type : str
            One of 'additive' (`AdditiveAssignmentRule`),
            'assignment' (`GeneralAssignmentRule`), or 'ode'
            (`GeneralODERule`, always run at every timestep
            regardless of `rule_frequency`).
        rule_attributes : dict
            Rule parameters/attributes. For additive/assignment
            rules, the only attribute used is 'equation'.
        rule_frequency : str, default "repeated"
            When the rule fires; see `Rule.set_frequency_flag` for
            the supported values.
        input_printout : bool, default False
            If True, print the rule being created.

        Raises
        ------
        SyntaxError
            If `rule_type` is not a recognized type.
        """
        if input_printout:
            print("Rule Created with \n\trule_type = "+str(rule_type)+"\n\trule_attributes="+str(rule_attributes)+"\n\trule_frequency="+str(rule_frequency))

        self.initialized = False

        # Parse the rule by rule type
        if rule_type == 'additive':
            rule_object = AdditiveAssignmentRule()
        elif rule_type == 'assignment':
            rule_object = GeneralAssignmentRule()
        elif rule_type == "ode":
            rule_frequency = "dt"
            rule_object = GeneralODERule()
        else:
            raise SyntaxError('Invalid type of Rule: ' + rule_type)

        # Add species and params to model
        species_names, params_names = rule_object.get_species_and_parameters(rule_attributes, species2index = self.species2index, params2index=self.params2index)
        for s in species_names: pass #self._add_species(s) No species should be added here
        for p in params_names: self._add_param(p)

        # initialize the rule
        if 'type' in rule_attributes:
            rule_attributes.pop('type')
        rule_object.initialize(rule_attributes,self.species2index,self.params2index, rule_frequency = rule_frequency)
        # Add the rule to the right place
        self.repeat_rules.append(rule_object)
        self.rule_definitions.append((rule_type, rule_attributes, rule_frequency))


    def set_parameter(self, param_name, param_value):
        """Set the value of a parameter in the model.

        If the parameter does not already exist in the model, it is
        added (and a warning is logged).

        Parameters
        ----------
        param_name : str
            The parameter name.
        param_value : float
            The value to set.
        """
        if param_name not in self.params2index:
            logging.info('Warning! parameter '+ param_name+" does not show up in any currently defined reactions or rules.")
            self._add_param(param_name)

        param_index = self.params2index[param_name]
        self.params_values[param_index] = param_value

    def create_parameter(self, param_name, param_value):
        """Add a new parameter to the model and set its value.

        Parameters
        ----------
        param_name : str
            The parameter name.
        param_value : float
            The value to set.
        """
        self._add_param(param_name)
        self.set_parameter(param_name, param_value)

    def check_parameters(self):
        """Check that every parameter in the model has a value.

        Raises
        ------
        ValueError
            If any parameter has no value (is NaN).
        """
        error_string = "Unspecified Parameters: "
        unspecified_parameters = False
        for p in self.params2index:
            i = self.params2index[p]
            if np.isnan(self.params_values[i]):
                unspecified_parameters = True
                error_string += p+"="+str(self.params_values[i])+', '

        if unspecified_parameters:
            error_string += f" (params2index is {self.params2index}; param_values is {self.params_values})"
            raise ValueError(error_string[:-2])

    def check_species(self):
        """Check that every species in the model has a value.

        Any species without an initial condition is defaulted to 0
        and a warning is logged.
        """
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
                if sp != "":
                    self.update_array[self.species2index[sp],reaction_index] = reaction_update_dict[sp]

            for sp in delay_reaction_update_dict:
                if sp != "":
                    self.delay_update_array[self.species2index[sp],reaction_index] = delay_reaction_update_dict[sp]

        return self.update_array, self.delay_update_array

    def parse_model(self, filename, input_printout = False):
        """Parse the model from a (deprecated) bioscrape XML file.

        Fills in all local variables (propensities, delays, update
        arrays) and maps species/parameters to indices.

        .. deprecated:: 1.0.2.2
            Bioscrape XML is being replaced by SBML and will no
            longer be supported in a future version of the software.

        Parameters
        ----------
        filename : str or file
            The model file. If a string, the file is opened;
            otherwise it is assumed to already be a file handle.
        input_printout : bool, default False
            If True, print verbose status messages while parsing.
        """
        # open XML file from the filename and use BeautifulSoup to parse it
        warnings.warn("Deprecated Warning: Bioscrape XML is being replaced by SBML and will no longer be supported in a future version of the software.")

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

            self.create_reaction(reactants = reactants, products = products, propensity_type = propensity['type'],
                                 propensity_param_dict = propensity_param_dict, delay_reactants=delay_reactants, 
                                 delay_products=delay_products, delay_param_dict = delay_param_dict, 
                                 input_printout = input_printout)
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


    def has_delays(self):
        """Return True if any reaction in the model has a delay."""
        if self.has_delay:
            return True
        else:
            return False

    def get_params2index(self):
        """Return the dict mapping parameter names to indices."""
        return self.params2index

    def get_species2index(self):
        """Return the dict mapping species names to indices."""
        return self.species2index

    def get_species_list(self):
        """Return species names as a list, ordered by index."""
        l = [None] * self.get_number_of_species()
        for s in self.species2index:
            l[self.species2index[s]] = s
        return l

    def get_species_array(self):
        """Return the current species values as a `numpy.ndarray`."""
        return np.array(self.species_values)

    def get_param_list(self):
        """Return parameter names as a list, ordered by index."""
        l = [None] * self.get_number_of_params()
        for p in self.params2index:
            l[self.params2index[p]] = p
        return l

    def get_params(self):
        """Return the set of parameter names.

        Returns
        -------
        dict_keys of str
            The parameter names.
        """
        return self.params2index.keys()

    def get_number_of_params(self):
        """Return the number of parameters in the model."""
        return len(self.params2index.keys())

    def get_parameter_values(self):
        """Return the parameter values as a `numpy.ndarray`."""
        return self.params_values

    def get_parameter_dictionary(self):
        """Return a dict mapping parameter names to their values."""
        param_dict = {}
        keys = self.get_params()
        values = self.get_parameter_values()
        for (key, value) in zip(keys, values):
            param_dict[key] = value
        return param_dict

    def get_species(self):
        """Return the set of species names.

        Returns
        -------
        dict_keys of str
            The species names.
        """
        return self.species2index.keys()

    def get_species_dictionary(self):
        """Return a dict mapping species names to their values."""
        A = self.get_species_array()
        species_dict = {}
        for s in self.species2index:
            species_dict[s] = A[self.species2index[s]]
        return species_dict
        # return {(s, A[self.species2index[s]]) for s in self.species2index}

    def get_number_of_species(self):
        """Return the number of species in the model."""
        return len(self.species2index.keys())

    def set_params(self, param_dict):
        """Set parameter values.

        Parameters
        ----------
        param_dict : dict of str to float
            Maps parameter names to the values to set. Names not
            already in the model trigger a warning and are ignored.
        """
        param_names = set(self.params2index.keys())
        for p in param_dict:
            if p in param_names:
                self.params_values[self.params2index[p]] = param_dict[p]
            else:
                warnings.warn('Trying to set parameter that is not in model: %s'  % p)


    def set_species(self, species_dict):
        """Set initial species values.

        Parameters
        ----------
        species_dict : dict of str to float
            Maps species names to the values to set. Names not
            already in the model trigger a warning and are ignored.
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
        """Return the list of `Propensity` objects, one per reaction."""
        return self.propensities

    def get_delays(self):
        """Return the list of `Delay` objects, one per reaction."""
        return self.delays

    def get_reactions(self):
        """Return the internal list of reaction tuples.

        Each entry is
        `(propensity_object, delay_object, reaction_update_dict,
        delay_reaction_update_dict)`.
        """
        return self.reaction_list

    def get_rules(self):
        """Return the list of rule definition tuples.

        Each entry is `(rule_type, rule_attributes, rule_frequency)`.
        """
        return self.rule_definitions

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
        """Return the stoichiometric matrix for immediate changes.

        Returns
        -------
        numpy.ndarray
            A 2-D array with 1 row per species, 1 column per
            reaction.
        """
        return self.update_array

    cdef np.ndarray get_delay_update_array(self):
        """
        Get the stoichiometric matrix for changes that occur after a delay.
        :return: (np.ndarray) A 2-D array with 1 row per species, 1 column for each reaction.
        """
        return self.delay_update_array

    def py_get_delay_update_array(self):
        """Return the stoichiometric matrix for delayed changes.

        Returns
        -------
        numpy.ndarray
            A 2-D array with 1 row per species, 1 column per
            reaction.
        """
        return self.delay_update_array

    def get_param_index(self, param_name):
        """Return the index of a parameter, or -1 if not found.

        Parameters
        ----------
        param_name : str
            The parameter name.
        """
        if param_name in self.params2index:
            return self.params2index[param_name]
        return -1

    def get_species_index(self, species_name):
        """Return the index of a species, or -1 if not found.

        Parameters
        ----------
        species_name : str
            The species name.
        """
        if species_name in self.species2index:
            return self.species2index[species_name]
        return -1

    def get_param_value(self, param_name):
        """Return the value of a parameter.

        Parameters
        ----------
        param_name : str
            The parameter name.

        Raises
        ------
        LookupError
            If no parameter with that name exists in the model.
        """
        if param_name in self.params2index:
            return self.params_values[self.params2index[param_name]]
        else:
            raise LookupError('No parameter with name '+ param_name)

    def get_species_value(self, species_name):
        """Return the value of a species.

        Parameters
        ----------
        species_name : str
            The species name.

        Raises
        ------
        LookupError
            If no species with that name exists in the model.
        """
        if species_name in self.species2index:
            return self.species_values[self.species2index[species_name]]
        else:
            raise LookupError('No species with name '+ species_name)

    def parse_general_expression(self, instring):
        """Parse a symbolic expression string into a `Term` tree.

        Parameters
        ----------
        instring : str
            The expression to parse, in terms of this model's
            species and parameter names.

        Returns
        -------
        Term
            The root of the parsed expression tree.
        """
        return parse_expression(instring,self.species2index,self.params2index)


    def generate_sbml_model(self, stochastic_model = False, **keywords):
        """Build an SBML representation of this model.

        Parameters
        ----------
        stochastic_model : bool, default False
            If True, annotate reactions for stochastic simulation.
        **keywords
            Additional keyword arguments forwarded to
            `create_sbml_model`.

        Returns
        -------
        document : libsbml.SBMLDocument
            The generated SBML document.
        model : libsbml.Model
            The generated SBML model.
        """
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
            if delay_type != None:
                delay_dict = {'type':delay_type, 'reactants':delay_reactants, 
                            'products':delay_products, 'parameters':delay_param_dict}
            else:
                delay_dict = None
            add_reaction(model, reactants, products, rxn_id, propensity_type,
                         propensity_param_dict, stochastic = stochastic_model,
                         delay_annotation_dict = delay_dict)
            rxn_count += 1

        rule_count = 0
        for rule_tuple in self.rule_definitions:
            rule_id = "rule" + str(rule_count)
            # Syntax of rule_tuple = (rule_type, rule_dict, rule_frequency)
            (rule_type, rule_dict, rule_frequency) = rule_tuple
            # Extract the rule variable id from rule_dict:
            
            if rule_type in ["ode", "ODE", 'GeneralODERule']:
                rule_formula = rule_dict['equation']
                rule_variable = rule_dict['target']
            else:
                equation = rule_dict['equation']
                split_eqn = [s.strip() for s in equation.split('=') ]
                try:
                    assert(len(split_eqn) == 2) # Checking rule_dict equation structure.
                except AssertionError as e:
                    e.args += ('rule equation', equation, 'not of the form VARIABLE = F(X).')
                    raise

                # Extract the rule formula for the variable above from rule_dict:
                rule_formula = split_eqn[1]
                rule_variable = split_eqn[0]
            add_rule(model, rule_id, rule_type, rule_variable, rule_formula, rule_frequency)
            rule_count += 1

        if document.getNumErrors():
            warnings.warn('The generated SBML model has errors:')
            err_message = document.getErrorLog().toString()
            print(err_message)
        return document, model

    def write_sbml_model(self, file_name, stochastic_model = False, **keywords):
        """Write this model to a file in SBML format.

        Parameters
        ----------
        file_name : str
            Path to write the SBML file to.
        stochastic_model : bool, default False
            If True, annotate reactions for stochastic simulation.
        **keywords
            Additional keyword arguments forwarded to
            `generate_sbml_model`.

        Returns
        -------
        bool
            True on success.
        """
        document, _ = self.generate_sbml_model(stochastic_model = stochastic_model, **keywords)
        sbml_string = libsbml.writeSBMLToString(document)
        with open(file_name, 'w') as f:
            f.write(sbml_string)
        return True

    def get_reaction_strings(self):
        """Return a list of strings describing the model's reactions.

        Not delay-compatible.

        Returns
        -------
        list of str
            One string per reaction, formatted as
            'reactants ->[propensity] products'.
        """
        reaction_strings = []

        # Propensity objects can only report their string descriptions in terms of 
        # species and parameter indices in this model. We'll need to fill those in with
        # concrete species and parameter names using this dictionary.
        decoding_args = {}
        for sp, idx in self.species2index.items():
            decoding_args[f's{idx}'] = sp
        for param, idx in self.params2index.items():
            decoding_args[f'p{idx}'] = param

        for reaction_index in range(len(self.reaction_list)):
            prop_object, delay_object, reaction_update_dict, delay_reaction_update_dict \
                = self.reaction_list[reaction_index]
            reactant_string = ""
            product_string = ""
            for sp, stoich in reaction_update_dict.items():
                if stoich > 0:
                    if len(reactant_string) > 0:
                        reactant_string += " + "
                    if stoich > 1:
                        reactant_string += f" {stoich} *"
                    reactant_string += f" {sp}"
                elif stoich < 0:
                    if len(product_string) > 0:
                        product_string += " + "
                    if stoich < -1:
                        product_string += f" {stoich} *"
                    product_string += f" {sp}"
            propensity_string = str(prop_object).format(**decoding_args)
            reaction_strings.append(f"{reactant_string} ->[{propensity_string}] {product_string}")
        return reaction_strings

    def get_derivative_strings(self):
        """Return each species' ODE dynamics as a string.

        Not delay-compatible.

        Returns
        -------
        dict of str to str
            Maps each species name to a string description of its
            ODE right-hand side.
        """
        # Propensity objects can only report their string descriptions in terms of
        # species and parameter indices in this model. We'll need to fill those in with
        # concrete species and parameter names using this dictionary.
        decoding_args = {}
        for sp, idx in self.species2index.items():
            decoding_args[f's{idx}'] = sp
        for param, idx in self.params2index.items():
            decoding_args[f'p{idx}'] = param

        deriv_strings = {sp: "" for sp in self.species2index.keys()}
        for reaction_index in range(len(self.reaction_list)):
            prop_object, delay_object, reaction_update_dict, delay_reaction_update_dict \
                = self.reaction_list[reaction_index]
            term_string = str(prop_object).format(**decoding_args)
            for sp, stoich in reaction_update_dict.items():
                if stoich > 0:
                    if len(deriv_strings[sp]) > 0:
                        deriv_strings[sp] += " + " 
                    if stoich > 1:
                        deriv_strings[sp] += f"{stoich} * "
                    deriv_strings[sp] += term_string
                elif stoich < 0:
                    deriv_strings[sp] += " - "
                if stoich < -1:
                    deriv_strings[sp] += f" {stoich} * "
                deriv_strings[sp] += term_string

        return deriv_strings

    # Update this if you change any of Model's member variables!
    def __getstate__(self):
        '''Returns the Model's state as a tuple of picklable Python objects. 

        Note that c_propensities, c_delays, and c_repeat_rules are just 
        pointers to the objects in propensities, delays, and repeat_rules,
        respectively, so only the latter are needed to fully represent the 
        Model's state.
        '''
        return (self._next_species_index,
                self._next_params_index,
                self._dummy_param_counter,
                self.has_delay,
                self.propensities,
                self.delays,
                self.repeat_rules,
                self.species2index,
                self.params2index,
                self.species_values,
                self.params_values,
                self.update_array,
                self.delay_update_array,
                self.reaction_list,
                self.reaction_updates,
                self.delay_reaction_updates,
                self.initialized,
                self.reaction_definitions,
                self.rule_definitions)

    # Update this if you change any of Model's member variables!
    def __setstate__(self, state):
        '''Sets this Model's state to that of another, using a tuple generated
        by the other Model's __getstate__ method. 

        Note that c_propensities, c_delays, and c_repeat_rules are just 
        pointers to the objects in propensities, delays, and repeat_rules,
        respectively, so only the latter are needed to fully reconstruct a 
        Model's state.
        '''
        self._next_species_index = state[0]
        self._next_params_index = state[1]
        self._dummy_param_counter = state[2]
        self.has_delay = state[3]

        self.propensities = state[4]
        self.c_propensities.clear()
        if state[4] is not None:
            for x in state[4]:
                self.c_propensities.push_back(<void *> x)
        self.delays = state[5]
        self.c_delays.clear()
        if state[5] is not None:
            for x in state[5]:
                self.c_delays.push_back(<void *> x)
        self.repeat_rules = state[6]
        self.c_repeat_rules.clear()
        if state[6] is not None:
            for x in state[6]:
                self.c_repeat_rules.push_back(<void *> x)

        self.species2index = state[7]
        self.params2index = state[8]
        self.species_values = state[9]
        self.params_values = state[10]
        self.update_array = state[11]
        self.delay_update_array = state[12]
        self.reaction_list = state[13]
        self.reaction_updates = state[14]
        self.delay_reaction_updates = state[15]
        self.initialized = state[16]
        self.reaction_definitions = state[17]
        self.rule_definitions = state[18]



##################################################                ####################################################
######################################              DATA    TYPES                       ##############################
#################################################                     ################################################

cdef class Schnitz:
    """The data acquired from a single cell trajectory.

    Parents and daughters are left as None and must be set later if
    required.

    Parameters
    ----------
    time : numpy.ndarray
        1-D array with time points.
    data : numpy.ndarray
        2-D array with one row for each time point, one column for
        each measured output.
    volume : numpy.ndarray
        1-D array with volume at each time point.

    Attributes
    ----------
    data : numpy.ndarray
        2-D array with one row for each time point, one column for
        each measured output.
    time : numpy.ndarray
        1-D array with the time points.
    volume : numpy.ndarray
        1-D array with the volume at each time point.
    parent : Schnitz
        The parent cell's Schnitz, or None if it has no parent.
    daughter1, daughter2 : Schnitz
        The daughter cells' Schnitzes, or None if they don't exist.
    """
    def __init__(self, time, data, volume):
        """See class docstring."""
        self.parent = None
        self.daughter1 = None
        self.daughter2 = None
        self.time = time
        self.volume = volume
        self.data = data

    def py_get_data(self):
        """Return the 2-D array of measured outputs."""
        return self.data

    def py_set_data(self, data):
        """Set the 2-D array of measured outputs."""
        self.data = data

    def py_get_time(self):
        """Return the 1-D array of time points."""
        return self.time

    def py_get_volume(self):
        """Return the 1-D array of volumes."""
        return self.volume

    def py_get_dataframe(self, Model = None):
        """Return this Schnitz's data as a pandas DataFrame.

        Parameters
        ----------
        Model : Model, optional
            If given, used to label the data columns with species
            names. Otherwise the columns are left unnamed.

        Returns
        -------
        pandas.DataFrame or numpy.ndarray
            The data, time, and volume as a DataFrame, or (if
            pandas is not installed) the raw data array.
        """
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
            return self.py_get_data()

    def py_get_parent(self):
        """Return the parent cell's Schnitz, or None."""
        return self.parent

    def py_get_daughters(self):
        """Return the `(daughter1, daughter2)` Schnitz tuple."""
        return (self.daughter1, self.daughter2)

    def py_set_parent(self, Schnitz p):
        """Set the parent cell's Schnitz."""
        self.set_parent(p)

    def py_set_daughters(self,Schnitz d1, Schnitz d2):
        """Set the daughter cells' Schnitzes."""
        self.set_daughters(d1,d2)

    def get_sub_lineage(self, dict species_dict = None):
        """Return a `Lineage` rooted at this Schnitz.

        Includes this Schnitz and all of its descendants.

        Parameters
        ----------
        species_dict : dict, optional
            If given, a species-name-to-index dict, and the result
            is an `ExperimentalLineage` using it. Otherwise a plain
            `Lineage` is returned.

        Returns
        -------
        Lineage
            The sub-lineage rooted at this Schnitz.
        """
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
        """Return a copy of this Schnitz.

        `time`, `data`, and `volume` are deep-copied; `parent` and
        the daughters are shared with the original.
        """
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
    """A cell lineage consisting of many `Schnitz` objects.

    Attributes
    ----------
    schnitzes : list of Schnitz
        All of the Schnitzes in the lineage.
    """
    def __init__(self):
        """See class docstring."""
        self.schnitzes = []

    def py_size(self):
        """Return the total number of Schnitzes in the lineage."""
        return self.c_schnitzes.size()

    def py_get_schnitz(self, unsigned index):
        """Return a specific Schnitz from the lineage.

        Parameters
        ----------
        index : unsigned
            The Schnitz to retrieve, 0 <= index < py_size().

        Returns
        -------
        Schnitz
            The requested Schnitz.

        Raises
        ------
        IndexError
            If `index` is out of range.
        """
        if index >= self.py_size():
            raise IndexError(f"index {index} > lineage.py_size() = {self.py_size()}")

        return (<Schnitz> (self.c_schnitzes[index]))

    def py_add_schnitz(self, Schnitz s):
        """Add a Schnitz to the lineage."""
        self.add_schnitz(s)

    def __setstate__(self, state):
        self.schnitzes = []
        self.c_schnitzes.clear()
        for s in state[0]:
            self.add_schnitz(s)

    def __getstate__(self):
        return (self.schnitzes,)

    def truncate_lineage(self,double start_time, double end_time):
        """Return a copy of this lineage restricted to a time window.

        Schnitzes entirely outside `[start_time, end_time]` are
        dropped; remaining Schnitzes have their `time`/`data`/
        `volume` arrays sliced to the window, with parent/daughter
        links preserved between the kept Schnitzes.

        Parameters
        ----------
        start_time : float
            Start of the time window to keep.
        end_time : float
            End of the time window to keep.

        Returns
        -------
        Lineage
            The truncated lineage.
        """
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

    def get_schnitzes_by_generation(self):
        """Group this lineage's Schnitzes by generation.

        Returns
        -------
        list of list of Schnitz
            One list per generation; generation 0 holds the
            Schnitzes with no parent, generation 1 their daughters,
            and so on.
        """
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
    """A `Lineage` whose Schnitz data columns map to species names.

    Parameters
    ----------
    species_indices : dict of str to int, default {}
        Maps species names to their column index in each Schnitz's
        data array.
    """
    def __init__(self, dict species_indices={}):
        """See class docstring."""
        super().__init__()
        self.species_dict = species_indices

    def py_set_species_indices(self, dict species_indices):
        """Set the species-name-to-data-column-index mapping."""
        self.species_dict = species_indices.copy()

    def py_get_species_index(self, str species):
        """Return the data column index for a species.

        Parameters
        ----------
        species : str
            The species name.

        Returns
        -------
        int
            The column index, or -1 (with a warning) if `species`
            is not in the mapping.
        """
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
        """Return a copy of this lineage restricted to a time window.

        Same as `Lineage.truncate_lineage`, but preserves the
        species-index mapping on the returned `ExperimentalLineage`.

        Parameters
        ----------
        start_time : float
            Start of the time window to keep.
        end_time : float
            End of the time window to keep.

        Returns
        -------
        ExperimentalLineage
            The truncated lineage.
        """
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










