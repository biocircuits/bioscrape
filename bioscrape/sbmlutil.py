import numpy as np
import sympy
from sympy.abc import _clash1
import warnings
from bioscrape.types import Model
import libsbml
from collections import OrderedDict # Need this to remove duplicates from lists

def read_model_from_sbml(sbml_file):
    return import_sbml(sbml_file)

def import_sbml(sbml_file, bioscrape_model = None, input_printout = False, **kwargs):
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
    if model is None:
        raise ValueError("SBML File {0} not found. Model could not be read.".format(sbml_file))
    if 'sbml_warnings' in kwargs:
        sbml_warnings = kwargs.get('sbml_warnings')
    elif 'bioscrape' in model.getId() or 'biocrnpyler' in model.getId():
        sbml_warnings = False
    else:
        sbml_warnings = True
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
    for reaction in model.getListOfReactions():
        # get the propensity
        kl = reaction.getKineticLaw()


        # capture any local parameters
        # also must save renamed local parameters to rename annotations later
        renamed_params = {}
        for p in kl.getListOfParameters():
            pid = p.getId()
            if pid in allparams:
                # If local parameter ID already exists in allparams due to another local/global parameter with same ID
                oldid = pid
                newid = oldid + '_' + reaction.getId()
                # Rename the ID everywhere it's used (such as in the Kinetic Law)
                kl.renameSIdRefs(oldid, newid)
                p.setId(newid)
                renamed_params[oldid] = newid #save the oldid-->newid mapping
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
        #We should no longer add underscores to parameters
        #rate_string = _add_underscore_to_parameters(kl_formula, allparams)
        rate_string = kl_formula

        if reaction.getReversible() and sbml_warnings:
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
                if np.isfinite(reactant.getStoichiometry()):
                    for i in range(int(reactant.getStoichiometry())):
                        reactant_list.append(reactantspecies_id)
                else:
                    reactant_list.append(reactantspecies_id)
            else:
                warnings.warn('Reactant in reaction {0} not found in the list of species in the SBML model.'
                 + ' The species will be added with zero initial amount'.format(reaction.getId()))
                allspecies[reactantspecies_id] = 0.0
                if np.isfinite(reactant.getStoichiometry()):
                    for i in range(int(reactant.getStoichiometry())):
                        reactant_list.append(reactantspecies_id)
                else:
                    reactant_list.append(reactantspecies_id)
        for product in reaction.getListOfProducts():
            productspecies = model.getSpecies(product.getSpecies())
            productspecies_id = productspecies.getId()
            if productspecies_id in allspecies:
                if np.isfinite(product.getStoichiometry()):
                    for i in range(int(product.getStoichiometry())):
                        product_list.append(productspecies_id)
                else:
                    product_list.append(productspecies_id)
            else:
                warnings.warn('Reactant in reaction {0} not found in the list of species in the SBML model.' +
                            ' The species will be added with zero initial amount'.format(reaction.getId()))
                allspecies[productspecies_id] = 0.0
                if np.isfinite(product.getStoichiometry()):
                    for i in range(int(product.getStoichiometry())):
                        product_list.append(productspecies_id)
                else:
                    product_list.append(productspecies_id)

        #Identify propensities based upon annotations
        annotation_string = reaction.getAnnotationString()
        if "PropensityType" in annotation_string:
            ind0 = annotation_string.find("<PropensityType>")
            ind1 = annotation_string.find("</PropensityType>")
            if ind0 == -1 or ind1 == -1:
                # Annotation could not be read
                if input_printout:
                    print('Annotation could not be read properly, adding reaction with general propensity.')
                propensity_type = 'general'
                general_kl_formula = {}
                general_kl_formula['rate'] = rate_string
                rxn = (reactant_list, product_list, propensity_type, general_kl_formula)
            else:
                # propensity_definition = {}
                annotation_list = annotation_string[ind0:ind1].split(" ")
                key_vals = [(i.split("=")[0], i.split("=")[1]) for i in annotation_list if "=" in i]
                propensity_params = {}
                for (k, v) in key_vals:
                    #Change the name of a parameter if it was renamed earlier
                    if v in renamed_params:
                        v = renamed_params[v]
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
                print("Propensity found with general ratestring:", rate_string)

        allreactions.append(rxn)

    # Go through rules one at a time
    allrules = []
    #"Rules must be a tuple: (rule_type (string), rule_attributes (dict), rule_frequency (optional))")
    for rule in model.getListOfRules():
        rule_formula = libsbml.formulaToL3String(rule.getMath())
        rulevariable = rule.getVariable()
        if rulevariable in allspecies:
            #rule_string = rulevariable + '=' + _add_underscore_to_parameters(rule_formula,allparams)
            rule_string = rulevariable + '=' + rule_formula
        elif rulevariable in allparams:
            #rule_string = '_' + rulevariable + '=' + _add_underscore_to_parameters(rule_formula,allparams)
            rule_string = rulevariable + '=' + rule_formula
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
            #rate_rule_formula = _add_underscore_to_parameters(rule_formula, allparams)
            rate_rule_formula = rule_formula
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
    
    # Check and warn if there are any unrecognized components (function definitions, packages, etc.)
    if len(model.getListOfCompartments()) > 0 or len(model.getListOfUnitDefinitions()) > 0  or len(model.getListOfEvents()) > 0: 
        if sbml_warnings:
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

def _get_species_list_in_formula(formula, species):
    sympy_rate = sympy.sympify(formula, _clash1)
    nodes = [sympy_rate]
    index = 0
    while index < len(nodes):
        node = nodes[index]
        index += 1
        nodes.extend(node.args)
    species_return = []
    for node in nodes:
        if type(node) == sympy.Symbol:
            if node.name in species:
                species_return.append(node.name)
    return species_return

def _remove_underscore_from_parameters(formula, parameters):
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
                node.name = node.name.replace('_','',1)
    return str(sympy_rate).replace('**','^')

def create_sbml_model(compartment_id="default", time_units='second', extent_units='mole', substance_units='mole',
                      length_units='metre', area_units='square_metre', volume_units='litre', volume = 1e-6):
    # Create an empty SBML Document of Level 3 Version 2 of SBML
    document = libsbml.SBMLDocument(3, 2) 
    model = document.createModel()
    model.setId('bioscrape_generated_model_' + str(np.random.randint(1e6)))
    # Define units for area (not used, but keeps COPASI from complaining)
    unitdef = model.createUnitDefinition()
    unitdef.setId('square_metre')
    unitdef.setName('square_metre')
    unit = unitdef.createUnit()
    unit.setKind(libsbml.UNIT_KIND_METRE)
    unit.setExponent(2)
    unit.setScale(0)
    unit.setMultiplier(1)

    # Set up required units and containers
    model.setTimeUnits(time_units)  # set model-wide time units
    model.setExtentUnits(extent_units)  # set model units of extent
    model.setSubstanceUnits(substance_units)  # set model substance units
    model.setLengthUnits(length_units)  # area units (never used?)
    model.setAreaUnits(area_units)  # area units (never used?)
    model.setVolumeUnits(volume_units)  # default volume unit

    # Define the default compartment
    compartment = model.createCompartment()
    compartment.setId(compartment_id)
    compartment.setName(compartment_id)
    compartment.setConstant(True)  # keep compartment size constant
    compartment.setSpatialDimensions(3)  # 3 dimensional compartment
    compartment.setVolume(volume)  # 1 microliter

    # Returning document is enough. document.getModel() gives the model, and model.getCompartment(0) gives the compartment.
    return document, model


# Creates an SBML id from a chemical_reaction_network.species object
def species_sbml_id(species_name, document=None):
    # Construct the species ID
    all_ids = []
    if document:
        all_ids = getAllIds(document.getListOfAllElements())

    trans = SetIdFromNames(all_ids)
    species_id = trans.getValidIdForName(species_name)
    return species_id


# Helper function to add a species to the model
# species must be chemical_reaction_network.species objects
def add_species(model, compartment, species, debug=False, initial_concentration=None):
    model = model  # Get the model where we will store results

    # Construct the species name
    species_name = species

    # Construct the species ID
    species_id = species_sbml_id(species_name, model.getSBMLDocument())
    if species_name != species_id:
        raise ValueError('Species names used are invalid strings to write into an SBML file.' + 
                        'Remove colons, semicolons, and other special characters.' + 
                        'Duplicate species names are also not allowed.' + 
                        'Starting species names with numbers is also not allowed')

    if debug: print("Adding species", species_name, species_id)
    sbml_species = model.createSpecies()
    sbml_species.setName(species_name)
    sbml_species.setId(species_id)
    sbml_species.setName(species_id)
    sbml_species.setCompartment(compartment.getId())
    sbml_species.setConstant(False)
    sbml_species.setBoundaryCondition(False)
    sbml_species.setHasOnlySubstanceUnits(False)
    if initial_concentration is None:
        initial_concentration = 0
    sbml_species.setInitialConcentration(initial_concentration)
    return sbml_species


# Helper function to add a parameter to the model
def add_parameter(model, param_name, param_value, debug=False):
    # Check to see if this parameter is already present
    parameter = model.createParameter()
    # all_ids = getAllIds(model.getSBMLDocument().getListOfAllElements())
    # trans = SetIdFromNames(all_ids)
    # parameter.setId(trans.getValidIdForName(param_name))
    if debug: print("Adding parameter", param_name)
    # param_name might be an invalid SBML identifier. But, if we change the name here
    # then we need to make sure the changes are propagated to KineticLaws etc. TODO.
    if param_name[0] == '_':
        param_name = param_name.replace('_','',1)
    parameter.setId(param_name)
    parameter.setName(param_name)
    # Set the value of the parameter
    parameter.setValue(float(param_value))
    parameter.setConstant(True) # Is set to False while creating rules (if some parameter is changed using Rules)

    return parameter

# Helper function to add a rule to an sbml model
# rule must be a chemical_reaction_network.reaction object
# propensity params is a dictionary of the parameters for non-massaction propensities.
def add_rule(model, rule_id, rule_type, rule_variable, rule_formula, **kwargs):
    # Create SBML equivalent of bioscrape rule:
    if rule_type == 'algebraic':
        raise NotImplementedError
    if rule_type == 'assignment' or rule_type == 'additive':
        # Simply create SBML assignment rule type. For additive rule type as well,
        # AssignmentRule type of SBML will work as $s_0$ is the artificial species that
        # exists in the bioscrape model.
        if rule_variable[0] == '_':
            rule_variable = rule_variable.replace('_','',1)
        for param in model.getListOfParameters():
            if rule_variable == param.getId():
                param.setConstant(False)
        allparams = {}
        for p in model.getListOfParameters():
            pid = p.getId()
            pid = '_' + pid
            allparams[pid] = p.getValue()
        rule_formula = _remove_underscore_from_parameters(rule_formula, allparams)
        rule = model.createAssignmentRule()
        rule.setId(rule_id)
        rule.setName(rule_id)
        rule.setVariable(rule_variable)
        rule.setFormula(rule_formula)
    return rule


# Helper function to add a reaction to an sbml model
# propensity params is a dictionary of the parameters for non-massaction propensities.
# propensity_params is a dictionary with keyword 'rate' for general propensity
def add_reaction(model, inputs_list, outputs_list,
                 reaction_id, propensity_type, propensity_params,
                 stochastic = False, propensity_annotation = True):

    # Create the reaction
    # We cast to an OrderedDict and back to remove duplicates.
    # We could cast to a regular dict instead, but only in Python 3.7 or higher.
    inputs = list(OrderedDict.fromkeys(inputs_list))
    #inputs.sort()
    input_coefs = [inputs_list.count(i) for i in inputs]
    outputs = list(OrderedDict.fromkeys(outputs_list))
    output_coefs = [outputs_list.count(o) for o in outputs]

    reaction = model.createReaction()
    reaction.setReversible(False)
    # reaction.setFast(False) # Deprecated in SBML
    all_ids = getAllIds(model.getSBMLDocument().getListOfAllElements())
    trans = SetIdFromNames(all_ids)
    reaction.setId(trans.getValidIdForName(reaction_id))
    reaction.setName(reaction.getId())
    ratestring = "" #Stores the string representing the rate function
    annotation_dict = {"type":propensity_type}
    allspecies = []
    for s in model.getListOfSpecies():
        sid = s.getId()
        allspecies.append(sid)

    allparams = {}
    for p in model.getListOfParameters():
        pid = p.getId()
        pid = '_' + pid 
        allparams[pid] = p.getValue()
    ratelaw = reaction.createKineticLaw()
    #Create Local Propensity Parameters
    if propensity_type=="massaction":
        annotation_dict["k"] = propensity_params['k']
        ratestring = propensity_params['k']

    #Hill Function Propensities
    elif propensity_type in ["hillpositive", "hillnegative", "proportionalhillpositive", "proportionalhillnegative"]:
        ratestring = propensity_params['k']
        annotation_dict["k"] = propensity_params['k']
        annotation_dict["K"] = propensity_params['K']
        annotation_dict["n"] = propensity_params['n']

    elif propensity_type == "general":
        pass
        #annotation_dict["rate"] = propensity_params['rate']
    else:
        raise ValueError(propensity_type+" is not a supported propensity_type")

    # Create the reactants
    reactants_list = []
    for i in range(len(inputs)):
        species = str(inputs[i]).replace("'", "")
        stoichiometry = input_coefs[i]
        # Multiple species with same name should be an invalid bioscrape construct.
        species_id = getSpeciesByName(model,species).getId()
        reactant = reaction.createReactant()
        reactant.setSpecies(species_id)  # ! TODO: add error checking
        reactants_list.append(species_id)
        reactant.setConstant(False)
        if stoichiometry is None or stoichiometry is np.nan:
            stoichiometry = 1.0
        reactant.setStoichiometry(stoichiometry)

        #Create Rate-strings for massaction propensities
        if propensity_type=="massaction" and stochastic:
            for i in range(stoichiometry):
                if i > 0:
                    ratestring += f" * ( {species_id} - {i} )"
                else:
                    ratestring += f" * {species_id}"

        elif propensity_type=="massaction" and not stochastic:
            if stoichiometry > 1:
                ratestring += f" * {species_id}^{stoichiometry}"
            else:
                ratestring += f" * {species_id}"

    # Create the products
    products_list = []
    for i in range(len(outputs)):
        species = str(outputs[i]).replace("'", "")
        stoichiometry = output_coefs[i]
        product = reaction.createProduct()
        species_id = getSpeciesByName(model, species).getId()
        product.setSpecies(species_id)
        products_list.append(species_id)
        if stoichiometry is None or stoichiometry is np.nan:
            stoichiometry = 1.0
        product.setStoichiometry(stoichiometry)
        product.setConstant(False)


    #Create ratestring for non-massaction propensities
    if propensity_type == "hillpositive":
        if not ("s1" in propensity_params):
            raise ValueError("hillpositive propensities, p(s1; k, K, n) "
                    "= k*s1^n/(s1^n + K), require the following key in the propensity_params dictionary:"
                    "'s1':species (chemical_reaction_network.species)")

        s = str(propensity_params['s1']).replace("'", "")
        s_species_id = getSpeciesByName(model, s).getId()
        if s_species_id not in reactants_list and s_species_id not in products_list:
            modifier = reaction.createModifier()
            modifier.setSpecies(s_species_id)
        n = propensity_params['n']
        K = propensity_params['K']
        ratestring+=f"*{s_species_id}^n/({s_species_id}^{n}+{K})"

        annotation_dict["s1"] = s_species_id

    elif propensity_type == "hillnegative":
        if not ("s1" in propensity_params):
            raise ValueError("hillnegative propensities, "
                    "p(s1; k, K, n) = k*1/(s1^n + K), require the following key in the propensity_params dictionary:"
                    "'s1':species (chemical_reaction_network.species)")
        s = str(propensity_params['s1']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()
        if s_species_id not in reactants_list and s_species_id not in products_list:
            modifier = reaction.createModifier()
            modifier.setSpecies(s_species_id)
        n = propensity_params['n']
        K = propensity_params['K']
        ratestring+=f"/({s_species_id}^{n}+{K})"
        annotation_dict["s1"] = s_species_id

    elif propensity_type == "proportionalhillpositive":
        if not ("s1" in propensity_params and "d" in propensity_params):
            raise ValueError("proportionalhillpositive propensities, "
                "p(s1, d; k, K, n) = k*d*s1^n/(s1^n + K), require the following key in the propensity_params dictionary:"
                "'s1':species (chemical_reaction_network.species)"
                "'d':species (chemical_reaction_network.species), ")

        s = str(propensity_params['s1']).replace("'", "")
        d = str(propensity_params['d']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()
        if s_species_id not in reactants_list and s_species_id not in products_list:
            modifier = reaction.createModifier()
            modifier.setSpecies(s_species_id)
        d_species_id = getSpeciesByName(model,d).getId()
        if d_species_id not in reactants_list and d_species_id not in products_list:
            modifier = reaction.createModifier()
            modifier.setSpecies(d_species_id)
        n = propensity_params['n']
        K = propensity_params['K']

        ratestring+=f"*{d_species_id}*{s_species_id}^n/({s_species_id}^{n} + {K})"

        annotation_dict["s1"] = s_species_id
        annotation_dict["d"] = d_species_id

    elif propensity_type == "proportionalhillnegative":
        if not ("s1" in propensity_params and "d" in propensity_params):
            raise ValueError("proportionalhillnegative propensities, "
                "p(s1, d; k, K, n) = k*d/(s1^n + K), require the following key in the propensity_params dictionary:"
                "'s1':species (chemical_reaction_network.species)"
                "'d':species (chemical_reaction_network.species), ")

        s = str(propensity_params['s1']).replace("'", "")
        d = str(propensity_params['d']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()
        if s_species_id not in reactants_list and s_species_id not in products_list:
            modifier = reaction.createModifier()
            modifier.setSpecies(s_species_id)
        d_species_id = getSpeciesByName(model,d).getId()
        if d_species_id not in reactants_list and d_species_id not in products_list:
            modifier = reaction.createModifier()
            modifier.setSpecies(d_species_id)
        n = propensity_params['n']
        K = propensity_params['K']

        ratestring+=f"*{d_species_id}/({s_species_id}^{n}+{K})"

        annotation_dict["s1"] = s_species_id
        annotation_dict["d"] = d_species_id
    elif propensity_type == "general":
        ratestring = propensity_params['rate']

        species_list = _get_species_list_in_formula(ratestring, allspecies)
        for s in species_list:
            if s not in reactants_list and s not in products_list:
                modifier = reaction.createModifier()
                modifier.setSpecies(s)

    ratestring = _remove_underscore_from_parameters(ratestring, allparams)
    # Set the ratelaw to the ratestring
    math_ast = libsbml.parseL3Formula(ratestring)
    ratelaw.setMath(math_ast)

    #Add propensity annotation
    if propensity_annotation and propensity_type != "general":
        annotation_string = "<PropensityType>"
        for k in annotation_dict:
            annotation_string += " "+k + "=" + str(annotation_dict[k])
        annotation_string += "</PropensityType>"
        reaction.appendAnnotation(annotation_string)

    return reaction

#  # Returns a list of all ids from the given list of elements
#  #
def getAllIds(allElements):
    result = []
    if (allElements == None or allElements.getSize() == 0):
        return result

    for i in range(0, allElements.getSize()):
        current = allElements.get(i)
        if (current.isSetId() \
           and current.getTypeCode() != libsbml.SBML_LOCAL_PARAMETER):
            result.append(current.getId())
    return result

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


# !/usr/bin/env python
##
## @file    setIdFromNames.py
## @brief   Utility program, renaming all SIds that also has
##          names specified. The new id will be derived from
##          the name, with all invalid characters removed.
##
## @author  Frank T. Bergmann
##
##
## <!--------------------------------------------------------------------------
## This sample program is distributed under a different license than the rest
## of libSBML.  This program uses the open-source MIT license, as follows:
##
## Copyright (c) 2013-2017 by the California Institute of Technology
## (California, USA), the European Bioinformatics Institute (EMBL-EBI, UK)
## and the University of Heidelberg (Germany), with support from the National
## Institutes of Health (USA) under grant R01GM070923.  All rights reserved.
##
## Permission is hereby granted, free of charge, to any person obtaining a
## copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
##
## Neither the name of the California Institute of Technology (Caltech), nor
## of the European Bioinformatics Institute (EMBL-EBI), nor of the University
## of Heidelberg, nor the names of any contributors, may be used to endorse
## or promote products derived from this software without specific prior
## written permission.
## ------------------------------------------------------------------------ -->
##
##

import sys
import os.path
import time


# This class implements an identifier transformer, that means it can be used
# to rename all sbase elements.
class SetIdFromNames(libsbml.IdentifierTransformer):
    def __init__(self, ids):
        # call the constructor of the base class
        libsbml.IdentifierTransformer.__init__(self)
        # remember existing ids ...
        self.existingIds = ids

        # The function actually doing the transforming. This function is called

    # once for each SBase element in the model.
    def transform(self, element):
        # return in case we don't have a valid element
        if (element == None \
           or element.getTypeCode() == libsbml.SBML_LOCAL_PARAMETER):
            return libsbml.LIBSBML_OPERATION_SUCCESS

            # or if there is nothing to do
        if (element.isSetName() == False \
           or element.getId() == element.getName()):
            return libsbml.LIBSBML_OPERATION_SUCCESS

            # find the new id
        newId = self.getValidIdForName(element.getName())

        # set it
        element.setId(newId)

        # remember it
        self.existingIds.append(newId)

        return libsbml.LIBSBML_OPERATION_SUCCESS

    def nameToSbmlId(self, name):
        IdStream = []
        count = 0
        end = len(name)

        if '0' <= name[count] and name[count] <= '9':
            IdStream.append('x_')
        if '*' in name:
            IdStream.append('xx')
        for count in range(0, end):
            if (('0' <= name[count] and name[count] <= '9') or
                    ('a' <= name[count] and name[count] <= 'z') or
                    ('A' <= name[count] and name[count] <= 'Z')):
                IdStream.append(name[count])
            else:
                IdStream.append('_')
        Id = ''.join(IdStream)
        if (Id[len(Id) - 1] != '_'):
            return Id

        return Id[:-1]

    #
    # Generates the id out of the name, and ensures it is unique.
    # It does so by appending numbers to the original name.
    #
    def getValidIdForName(self, name):
        baseString = self.nameToSbmlId(name)
        id = baseString
        count = 1
        while (self.existingIds.count(id) != 0):
            id = "{0}_{1}".format(baseString, count)
            count = count + 1
        return id

    #  #


def getSpeciesByName(model, name, compartment=''):
    '''
    Returns a list of species in the Model with the given name
    compartment : (Optional) argument to specify the compartment name in which
    to look for the species.
    '''
    if type(name) is not str:
        raise ValueError('"name" must be a string.')
    species_found = []
    for species in model.getListOfSpecies():
        if species.getName() == name:
            if compartment != '':
                comp_elem = species.getCompartment()
                comp_name = model.getElementBySId(comp_elem).getName()
                if comp_name == compartment:
                    species_found.append(species)
                else:
                    continue
            else:
                species_found.append(species)

    if len(species_found) == 1:
        return species_found[0]
    elif not species_found:
        raise ValueError('The species ' + name + ' not found.')
    else:
        warnings.warn('Multiple species with name ' + name + ' found. Returning a list')
        return species_found