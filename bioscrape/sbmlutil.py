import numpy as np
import sympy
from sympy.abc import _clash1
import warnings
from bioscrape.types import Model

def read_model_from_sbml(sbml_file):
    return import_sbml(sbml_file)

def import_sbml(sbml_file, bioscrape_model = None, input_printout = False):
    """
    Convert SBML document to bioscrape Model object. Note that events, compartments, non-standard function definitions,
    and some kinds of rules will be ignored. 
    Adds mass action kinetics based reactions with the appropriate mass action propensity in bioscrape. 
    Propensities with the correct annotation are added as compiled propensity types.
    All other propensities are added as general propensity.
    Local parameters are renamed if there is a conflict since bioscrape does not have a local environment.
    """
    # Attempt to import libsbml and read the SBML model.
    try:
        import libsbml
    except:
        raise ImportError("libsbml not found. See sbml.org for installation help!\n" +
                        'If you are using anaconda you can run the following:\n' +
                            'conda install -c SBMLTeam python-libsbml\n\n\n')

    reader = libsbml.SBMLReader()
    doc = reader.readSBML(sbml_file)

    if doc.getNumErrors() > 1:
        raise SyntaxError("SBML File {0} cannot be read without errors".format(sbml_file))

    model = doc.getModel()

    # Parse through species and parameters and keep a set of both along with their values.
    allspecies = {}
    allparams = {}
    allreactions = []
    for s in model.getListOfSpecies():
        sid = s.getId()
        if sid == "volume" or sid == "t":
            warnings.warn("You have defined a species called '" + sid +
                            ". This species is being ignored and treated as a keyword.")
            continue
        allspecies[sid] = 0.0
        if np.isfinite(s.getInitialAmount()):
            allspecies[sid] = s.getInitialAmount()
        if np.isfinite(s.getInitialConcentration()) and allspecies[sid] == 0:
            allspecies[sid] = s.getInitialConcentration()

    for p in model.getListOfParameters():
        pid = p.getId()
        allparams[pid] = 0.0
        if np.isfinite(p.getValue()):
            allparams[pid] = p.getValue()

    # Now go through reactions one at a time to get stoich and rates, then append to reaction_list.
    reaction_list = []
    for reaction in model.getListOfReactions():
        # get the propensity 
        kl = reaction.getKineticLaw()
        # capture any local parameters
        for p in kl.getListOfParameters():
            pid = p.getId()
            if pid in allparams:
                # If local parameter ID already exists in allparams due to another local/global parameter with same ID
                oldid = pid
                newid = oldid + '_' + reaction.getId()
                # Rename the ID everywhere it's used (such as in the Kinetic Law)
                kl.renameSIdRefs(oldid, newid)
                p.setId(newid)
                # Rename its usages
                for element in reaction.getListOfAllElements():
                    element.renameSIdRefs(oldid, newid)
                pid = newid

            allparams[pid] = 0.0
            if np.isfinite(p.getValue()):
                allparams[pid] = p.getValue()

        # get the formula as a string and then add
        # a leading _ to parameter names
        kl_formula = libsbml.formulaToL3String(kl.getMath())
        rate_string = _add_underscore_to_parameters(kl_formula, allparams)

        if reaction.getReversible():
            warnings.warn('SBML model contains reversible reaction!\n' +
                            'Please check rate expressions and ensure they are non-negative before doing '+
                            'stochastic simulations.') 

        #Get Reactants and Products    
        reactant_list = []
        product_list = []
        for reactant in reaction.getListOfReactants():
            reactantspecies = model.getSpecies(reactant.getSpecies())
            reactantspecies_id = reactantspecies.getId()
            if reactantspecies_id in allspecies:
                reactant_list.append(reactantspecies_id)
            else:
                warnings.warn('Reactant in reaction {0} not found in the list of species in the SBML model.' + 
                            ' The species will be added with zero initial amount'.format(reaction.getId()))
                allspecies[reactantspecies_id] = 0.0
                reactant_list.append(reactantspecies_id)
            
        for product in reaction.getListOfProducts():
            productspecies = model.getSpecies(product.getSpecies())
            productspecies_id = productspecies.getId()
            if productspecies_id in allspecies:
                product_list.append(productspecies_id)
            else:
                warnings.warn('Reactant in reaction {0} not found in the list of species in the SBML model.' + 
                            ' The species will be added with zero initial amount'.format(reaction.getId()))
                allspecies[productspecies_id] = 0.0
                product_list.append(productspecies_id)

        #Identify propensities based upon annotations
        annotation_string = reaction.getAnnotationString()
        if "PropensityType" in annotation_string:
            ind0 = annotation_string.find("<PropensityType>")
            ind1 = annotation_string.find("</PropensityType>")
            propensity_definition = {}
            annotation_list = annotation_string[ind0:ind1].split(" ")
            key_vals = [(i.split("=")[0], i.split("=")[1]) for i in annotation_list if "=" in i]
            propensity_params = {}
            for (k, v) in key_vals:
                try:
                    propensity_params[k] = float(v)
                except ValueError:
                    propensity_params[k] = v
            if input_printout:
                print("Reaction found:", reactant_list, "->", product_list)
                print("Annotated propensity found with params:", propensity_params)
            rxn  = (reactant_list, product_list, propensity_params['type'], propensity_params)

        else: #No annotation found
            propensity_type = 'general'
            general_kl_formula = {}
            general_kl_formula['rate'] = rate_string 
            rxn = (reactant_list, product_list, propensity_type, general_kl_formula)
            if input_printout:
                print("Reaction found:", reactant_list, "->", product_list)
                print("Unnotated propensity found with ratestring:", rate_string)
        allreactions.append(rxn)
            
    # Go through rules one at a time
    allrules = []
    #"Rules must be a tuple: (rule_type (string), rule_attributes (dict), rule_frequency (optional))")
    for rule in model.getListOfRules():
        rule_formula = libsbml.formulaToL3String(rule.getMath())
        rulevariable = rule.getVariable()
        if rulevariable in allspecies:
            rule_string = rulevariable + '=' + _add_underscore_to_parameters(rule_formula,allparams)
        elif rulevariable in allparams:
            rule_string = '_' + rulevariable + '=' + _add_underscore_to_parameters(rule_formula,allparams)
        else:
            warnings.warn('SBML: Attempting to assign something that is not a parameter or species %s'
                            % rulevariable)
            continue
        if rule.getElementName() == 'algebraicRule':
            warnings.warn('Unsupported rule type: %s' % rule.getElementName())
            continue
        elif rule.getElementName() == 'assignmentRule':
            rule_type = 'assignment'
        elif rule.getElementName() == 'rateRule':
            rate_rule_formula = _add_underscore_to_parameters(rule_formula, allparams)
            rule_rxn = ([''], [rulevariable], 'general', rate_rule_formula) # Create --> X type reaction to model rate rules.
            allreactions.append(rule_rxn)
            continue
        else:
            raise ValueError('Invalid SBML Rule type.')
        rule_dict = {}
        rule_dict['equation'] = rule_string
        rule_frequency = 'repeated'
        rule_tuple = (rule_type, rule_dict, rule_frequency)
        allrules.append(rule_tuple)
    
    #print('allparams = {0}'.format(allparams))
    #print('allspecies = {0}'.format(allspecies))
    #print('allreactions = {0}'.format(allreactions))
    #print(allrules)

    # Check and warn if there are any unrecognized components (function definitions, packages, etc.)
    if len(model.getListOfCompartments()) > 0 or len(model.getListOfUnitDefinitions()) > 0  or len(model.getListOfEvents()) > 0: 
        warnings.warn('Compartments, UnitDefintions, Events, and some other SBML model components are not recognized by bioscrape. ' + 
                        'Refer to the bioscrape wiki for more information.')
    
    #If no Model is passed into the function, a Model is returned
    if bioscrape_model == None:
        bioscrape_model = Model()
    #If a Model is passed into the function, that Model is modified
    if isinstance(bioscrape_model, Model):
        for species in allspecies.keys():
            bioscrape_model._add_species(species)

        for (param, val) in allparams.items():
            bioscrape_model._add_param(param)
            bioscrape_model.set_parameter(param, val)
            if input_printout:
                print("Adding Parameter:", param, "=", val)

        for rxn in allreactions:
            if len(rxn) == 4:
                reactants, products, propensity_type, propensity_param_dict = rxn
                delay_type, delay_reactants, delay_products, delay_param_dict = None, None,  None, None
            elif len(rxn) == 8:
                reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict = rxn
            bioscrape_model.create_reaction(reactants, products, propensity_type, propensity_param_dict, delay_type, delay_reactants, delay_products, delay_param_dict, input_printout = input_printout)
        
        for rule in allrules:
            if len(rule) == 2:
                rule_type, rule_attributes = rule
                bioscrape_model.create_rule(rule_type, rule_attributes, input_printout = input_printout)
            elif len(rule) == 3:
                rule_type, rule_attributes, rule_frequency = rule
                bioscrape_model.create_rule(rule_type, rule_attributes, rule_frequency = rule_frequency, input_printout = input_printout)
        bioscrape_model.set_species(allspecies)
        bioscrape_model.py_initialize()
        return bioscrape_model
    else:
        raise ValueError("bioscrape_model keyword must be a Bioscrape Model object or None (in which case a Model object is returned).")

     

# Helpful utility functions start here 
def _add_underscore_to_parameters(formula, parameters):
    sympy_rate = sympy.sympify(formula, _clash1)
    nodes = [sympy_rate]
    index = 0
    while index < len(nodes):
        node = nodes[index]
        index += 1
        nodes.extend(node.args)

    for node in nodes:
        if type(node) == sympy.Symbol:
            if node.name in parameters:
                node.name = '_' + node.name

    return str(sympy_rate)


# Renames lists of SIds in an SBML Document
def renameSIds(document, oldSIds, newSIds, debug = False):
    '''
    Updates the SId from oldSId to newSId for any component of the Subsystem.
    Returns the SBMLDocument of the updated Subsystem
    '''

    #
    # @file    renameSId.py
    # @brief   Utility program, renaming a specific SId
    #          while updating all references to it.
    # @author  Frank T. Bergmann
    #
    # <!--------------------------------------------------------------------------
    # This sample program is distributed under a different license than the rest
    # of libSBML.  This program uses the open-source MIT license, as follows:
    #
    # Copyright (c) 2013-2018 by the California Institute of Technology
    # (California, USA), the European Bioinformatics Institute (EMBL-EBI, UK)
    # and the University of Heidelberg (Germany), with support from the National
    # Institutes of Health (USA) under grant R01GM070923.  All rights reserved.
    #
    # Permission is hereby granted, free of charge, to any person obtaining a
    # copy of this software and associated documentation files (the "Software"),
    # to deal in the Software without restriction, including without limitation
    # the rights to use, copy, modify, merge, publish, distribute, sublicense,
    # and/or sell copies of the Software, and to permit persons to whom the
    # Software is furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included in
    # all copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    # THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    # FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    # DEALINGS IN THE SOFTWARE.
    #
    # Neither the name of the California Institute of Technology (Caltech), nor
    # of the European Bioinformatics Institute (EMBL-EBI), nor of the University
    # of Heidelberg, nor the names of any contributors, may be used to endorse
    # or promote products derived from this software without specific prior
    # written permission.
    # ------------------------------------------------------------------------ -->
    #

    try:
        import libsbml
    except:
        raise ImportError("libsbml not found. See sbml.org for installation help!\n" +
                            'If you are using anaconda you can run the following:\n' +
                            'conda install -c SBMLTeam python-libsbml\n\n\n')


    if len(oldSIds) != len(newSIds):
        raise ValueError("Length oldSIds != length newSIds")

    for ind in range(len(oldSIds)):
        oldSId = oldSIds[ind]
        newSId = newSIds[ind]

        if oldSId == newSId:
            warnings.warn("The Ids are identical: " +str(oldSId)+". SId skipped.")

        if not libsbml.SyntaxChecker.isValidInternalSId(newSId):
            warnings.warn("The new SId '{0}' does not represent a valid SId.".format(newSId))
            

        element = document.getElementBySId(oldSId)

        if element == None:
            if debug:
                warnings.warn("Found no element with SId '{0}' in subsystem {1}".format(oldSId,document.getModel().getId()))

        # update all references to this element
        allElements = document.getListOfAllElements()
        for i in range(allElements.getSize()):
            current = allElements.get(i)
            current.renameSIdRefs(oldSId, newSId)
    return document