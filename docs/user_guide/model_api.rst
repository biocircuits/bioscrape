Model API
=========

This page is an overview of an internal model API for bioscrape. This
API allows for models to be constructed programmatically as well as
through SBML files.

Creating Models
---------------

Models can be created from SBML files, lists of species, reactions,
parameters, and rules, or a combination of the above.

.. _bioscrapetypesmodelinit-constructor:

bioscrape.types.Model.\ **init** [Constructor]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inputs:
'''''''

-  filename = None: the filename of a bioscrape xml file.
-  species = []: A list of species (str) to add to the model. Species
   will also be inferred from reactions, rules, and the
   initial_condition_dict parameter.
-  reactions = []: A list of reactions to add to the model. Reactions
   without delay are represented as tuples: (reactants (list of str),
   products (list of str), propensity_type (str), propensity_parameters
   (dictionary)) Reactions with delay are represented by larger tuples:
   (reactants (list of str), products (list of str), propensity_type
   (str), propensity_param_dict delay_type (str), delay_reactants (list
   of str), delay_products (list of str), delay_parameters
   (dictionary)). See create_reaction for more details.
-  parameters = []: A list of parameters to add to the model. Parameters
   are represented as a tuples: (param_name(str), param_value(float)).
-  rules = []: A list of rules to add to the model. Rules are
   represented as tuples: (rule_type (str), rule_attributes
   (dictionary), rule_frequency (str)). See create_rule for more
   details.
-  initial_condition_dict = {}: a dictionary {Specie(str) :
   initial_amount(int/double)}. Species not passed in this way default
   to a value of 0.
-  sbml_filename = (str), A string for the SBML file. This argument
   cannot be used together with the above model construction arguments
   in the Model constructor call. However, once you have created a model
   (either using an SBML or using species, reactions), you can edit the
   model however you like.

Creating Reactions
------------------

Reactions can be created automatically if they use the default supported
types with Model.create_reaction. Reactions with custom delay types and
propensity types (such as user defined delays and propensities) can be
added to the model with Model._add_reaction.

.. _modelcreate_reaction:

Model.create_reaction
~~~~~~~~~~~~~~~~~~~~~

Required Inputs:
''''''''''''''''

-  reactants (list): a list of reactant specie names (strs)
-  products (list): a list of product specie names (strs)
-  propensity_type (str): the type of propensity

   -  Supported types: "massaction", "hillpositive",
      "proportionalhillpositive", "hillnegative",
      "proportionalhillnegative", "general"

-  propensity_param_dict: a dictionary of parameters for the given
   propensity type

   -  Specifications of the attributes required by each propensity type
      can be found in the :doc:`propensity <propensities>` section.
   -  Note: parameter dictionaries can be of the form {key (str) ->
      var_name (str)} or {key (str) --> var_val (float or a string
      representation of a float)}. In the first case, the parameter
      var_name must be added with a value explicitly. In the second
      case, dummy variables are automatically created in the model and
      set to the value var_val.

Optional Inputs:
''''''''''''''''

-  delay_type = None: a str indicating the type of delay

   -  Supported types: "fixed", "gaussian", and "gamma"

-  delay_reactants = []: a list of delay reaction reactant specie names
   (str)
-  delay_products: a list of delay reaction products specie names (str)
-  delay_param_dict: a dictionary of the parameters for the delay
   distribution

   -  Specifications of the attributes required by each delay type can
      be found in the :doc:`delay <delays>` section.
   -  Note: parameter dictionaries can be of the form {key (str) ->
      var_name (str)} or {key (str) --> var_val (float or a string
      representation of a float)}. In the first case, the parameter
      var_name must be added with a value explicitly. In the second
      case, dummy variables are automatically created in the model and
      set to the value var_val.

.. _model_add_reaction:

Model._add_reaction
~~~~~~~~~~~~~~~~~~~

.. _inputs-1:

Inputs:
'''''''

-  reaction_update_dict: A dictionary species_index (accessible via
   Model.species2index) --> stochiometric change (a positive integer for
   species created by the reaction, a negative integer for species
   consumed by the reaction, and 0 for catalysts). Only species involved
   in the reaction as inputs or outputs should be in the dictionary.
-  propensity_object: An instance of a propensity object
-  propensity_param_dict: a dictionary of propensity object parameters
   used by the propensity_object.initialize function.

   -  Specifications of the attributes required by each propensity type
      can be found in the :doc:`propensity <propensities>` section.
   -  Note: parameter dictionaries can be of the form {key (str) ->
      var_name (str)} or {key (str) --> var_val (float or a string
      representation of a float)}. In the first case, the parameter
      var_name must be added with a value explicitly. In the second
      case, dummy variables are automatically created in the model and
      set to the value var_val.

.. _optional-inputs-1:

Optional Inputs:
''''''''''''''''

-  delay_reaction_update_dict = {}: same as reaction_update_dict but
   occurs with a delay.
-  delay_object = None: An instance of a delay object.
-  delay_param_dict = {}: a dictionary of delay object parameters used
   by the delay_object.initialize function.

   -  Specifications of the attributes required by each delay type can
      be found in the :doc:`delay <delays>` section.
   -  Note: parameter dictionaries can be of the form {key (str) ->
      var_name (str)} or {key (str) --> var_val (float or a string
      representation of a float)}. In the first case, the parameter
      var_name must be added with a value explicitly. In the second
      case, dummy variables are automatically created in the model and
      set to the value var_val.

Adding Rules
------------

Rules take the form of algebraic expressions involving species and
parameters which are updated periodically during simulation.

.. _modelcreate_rule:

Model.create_rule
~~~~~~~~~~~~~~~~~

Creates a rule using one of the supported rule types.

.. _inputs-2:

Inputs:
'''''''

-  rule_type (str): the type of rule

   -  Supported Types: 'additive' and 'assignment'

-  rule_attributes (dict): a dictionary of rule attributes and
   parameters

   -  Required attributes: 'Equation' --> Algebraic Rule Equation
      written using BioSCRAPE conventions

-  rule_frequency (str): frequency of how often the rule is updated.

   -  Supported Frequencies: "repeated" (updated every simulation step)

Helpful utility functions for species and parameters
----------------------------------------------------

Species are usually added automatically from reactions and rules.
However they can also be added directly and have their initial
conditions changed (the default is always 0).

.. _modelset_species:

Model.set_species
~~~~~~~~~~~~~~~~~

Sets the value of multiple species stored in a dictionary. Species not
already in the model will not be added this way. Recommended function
for changing initial conditions.

.. _inputs-3:

Inputs:
'''''''

-  species_dict: specie_name (str) --> species value (double)

.. _model_set_species_value:

Model._set_species_value
~~~~~~~~~~~~~~~~~~~~~~~~

Sets the value of one species and adds it to the model if it has not
already been added.

.. _inputs-4:

Inputs:
'''''''

-  specie (str): specie name
-  value (double): specie value

.. _model_add_species:

Model._add_species
~~~~~~~~~~~~~~~~~~

Adds a specie to the model with default value of 0.

.. _inputs-5:

Inputs:
'''''''

-  species (str): species name

Adding and Setting Parameters
-----------------------------

.. _modelset_params:

Model.set_params
~~~~~~~~~~~~~~~~

Sets the value of multiple parameters stored in a dictionary. Parameters
not already in the model will not be added this way.

.. _inputs-6:

Inputs:
'''''''

-  param_dict: param_name (str) --> param value (double)

.. _model_add_param:

Model._add_param
~~~~~~~~~~~~~~~~

Adds a parameter to the model, with the default value nan.

.. _inputs-7:

Inputs:
'''''''

-  parameter_name (str): The name of the parameter.

.. _modelset_parameter:

Model.set_parameter
~~~~~~~~~~~~~~~~~~~

Sets the value of a parameter and adds it to the model (if its not
already in the model).

Inputs:

-  param_name (str): the name of the parameter
-  param_value (double): the value of the parameter

.. _modelget_parameter_dictionary:

Model.get_parameter_dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the dictionary of parameter names and values in the model.

Output:

-  params_dict (dict): the dictionary with key as parameter names and
   values for the parameter values

.. _modelget_species_dictionary:

Model.get_species_dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the dictionary of species names and initial values in the model.

Output:

-  species_dict (dict): the dictionary with key as species names and
   values for the species initial values
