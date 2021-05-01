import warnings, os, pickle

# We don't want warnings in dependencies to show up in bioscrape's tests.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	import numpy as np
	import pylab as plt
	import random
	import pytest

import test_utils
import pandas as pd
from bioscrape.simulator import *
from bioscrape.types import *
from bioscrape.lineage import *

# Seed RNG value. All tests use this value.
seed = 54173

# Set true to get more diagnostic prints
debug = False


def test_random_complex_model():
	'''
	Creates a two-reaction Model with randomish (but replicable) 
	propensities, delays, and rules, then tests whether that model produces
	the same results after serialization (pickling) and deserialization 
	(unpickling). 

	This test currently does NOT distinguish between pickling and depickling
	errors -- it just checks that if you pickle and depickle a Model, it still
	produces the same results when simulated.
	'''

	test_utils.set_seed(seed)

	################
	# CREATE MODEL #
	################

	# Start with a basic reaction A+B-->C
	inputs = ["A", "B"]
	outputs = ["C"]
	all_species = inputs + outputs
	param_min, param_max = -4, 4

	if debug:
		print('simulating propensity type ', prop_type)

	
	# We will use a random(ish) rational function for propensity
	param_dict = {}
	rate_str = "(1+"
	numerator_terms = np.random.randint(0, 5)
	denominator_terms = np.random.randint(0, 5)
	for i in range(numerator_terms):
		coef = str(round(np.exp(np.random.uniform(low = param_min,
												 high = param_max)), 3))
		exp = str(round(np.random.uniform(low = 0,high = param_max), 3))
		species = all_species[np.random.randint(len(all_species))]
		rate_str += coef + "*" + species + "^" + exp + "+"
	rate_str = rate_str[:-1] + ")"
	rate_str += "/(1+"
	for i in range(denominator_terms):
		coef =str(round(np.exp(np.random.uniform(low = param_min,
												 high = param_max)), 3))
		exp = str(round(np.random.uniform(low = 0,high = param_max), 3))
		species = all_species[np.random.randint(len(all_species))]
		rate_str += coef + "*" + species + "^" + exp + "+"
	rate_str =  rate_str[:-1] + ")"
	param_dict['rate'] = rate_str
	if debug:
		print('\t params =', param_dict)

	consolidation_rxn = (inputs, outputs, "general", param_dict)

	# We'll add an A-->B reaction with delay.
	delay_k = round(np.exp(np.random.uniform(low = -1,
                                             high = 1)), 3)
	delay_params = dict()
	for p in ["mean", "std"]:
		delay_params[p] = round(np.exp(np.random.uniform(low = param_min, 
                                                         high = param_max)), 3)
	conversion_rxn = (["A"], [], "massaction", {"k": delay_k}, "gaussian", 
                      [], ["B"], delay_params)

	x0 = {"A":25, "B": 25, "C":0, "I": 0}
	M = Model(species = ["A", "B", "C", "I"], 
			  reactions = [consolidation_rxn, conversion_rxn], 
			  initial_condition_dict = x0)
	M.set_species(x0)


    # Finally, to test assignment rules, we'll add a rule adding an "inducer" at 
    # time 75.
	M.create_parameter("I0", 10) #Inducer concentration
	M.create_parameter("T_I0", 2) #Initial time inducer is added
	M.create_rule("assignment", {"equation":"I = _I0*Heaviside(t-_T_I0)"})


	###############################
	# TEST SERIALIZATION OF MODEL #
	###############################
	timepoints = np.arange(0, 25, 0.01)
	for is_stochastic in [False, True]:
		prepickling_results = py_simulate_model(timepoints, 
									Model = M, 
									stochastic = is_stochastic, 
									return_dataframe = False).py_get_result()
		pickled_M = pickle.dumps(M)
		postpickled_M = pickle.loads(pickled_M)

		test_utils.set_seed(seed)
		postpickling_results = py_simulate_model(timepoints, 
									Model = postpickled_M, 
									stochastic = is_stochastic, 
									return_dataframe = False).py_get_result()

		assert np.allclose(prepickling_results, postpickling_results, equal_nan = True), \
			   	f"Complex Model changes in " + \
			   	f"{'stochastic' if is_stochastic else 'deterministic'} " + \
			   	"simulation when pickled and unpickled."


def test_random_complex_lineagemodel():
	'''
	Creates a two-reaction LineageModel with more or less the same settings as
	the Model in test_random_complex_model, plus linage_propensities, a volume 
	event, a division rule, and a death rule.

	This test currently does NOT distinguish between pickling and depickling
	errors -- it just checks that if you pickle and depickle a Model, it still
	produces the same results when simulated.
	'''

	test_utils.set_seed(seed)

	################
	# CREATE MODEL #
	################

	# Start with a basic reaction A+B-->C
	inputs = ["A", "B"]
	outputs = ["C"]
	all_species = inputs + outputs
	param_min, param_max = -4, 4

	if debug:
		print('simulating propensity type ', prop_type)

	
	# We will use a random(ish) rational function for propensity
	param_dict = {}
	rate_str = "(1+"
	numerator_terms = np.random.randint(0, 5)
	denominator_terms = np.random.randint(0, 5)
	for i in range(numerator_terms):
		coef = str(round(np.exp(np.random.uniform(low = param_min,
												 high = param_max)), 3))
		exp = str(round(np.random.uniform(low = 0,high = param_max), 3))
		species = all_species[np.random.randint(len(all_species))]
		rate_str += coef + "*" + species + "^" + exp + "+"
	rate_str = rate_str[:-1] + ")"
	rate_str += "/(1+"
	for i in range(denominator_terms):
		coef =str(round(np.exp(np.random.uniform(low = param_min,
												 high = param_max)), 3))
		exp = str(round(np.random.uniform(low = 0,high = param_max), 3))
		species = all_species[np.random.randint(len(all_species))]
		rate_str += coef + "*" + species + "^" + exp + "+"
	rate_str =  rate_str[:-1] + ")"
	param_dict['rate'] = rate_str
	if debug:
		print('\t params =', param_dict)


	### DELETEME
	param_dict["rate"] = "A*3"

	consolidation_rxn = (inputs, outputs, "general", param_dict)

	# We'll add an A-->B reaction with delay.
	delay_k = round(np.exp(np.random.uniform(low = -1,
                                             high = 1)), 3)
	delay_params = dict()
	for p in ["mean", "std"]:
		delay_params[p] = round(np.exp(np.random.uniform(low = param_min, 
                                                         high = param_max)), 3)
	conversion_rxn = (["A"], [], "massaction", {"k": delay_k}, "gaussian", 
                      [], ["B"], delay_params)

	x0 = {"A":25, "B": 25, "C":0, "I": 0}
	M = LineageModel(species = ["A", "B", "C", "I"], 
			  reactions = [consolidation_rxn, conversion_rxn], 
			  initial_condition_dict = x0)
	M.set_species(x0)


    # Finally, to test assignment rules, we'll add a rule adding an "inducer" at 
    # time 100.
	M.create_parameter("I0", 10) #Inducer concentration
	M.create_parameter("T_I0", 2.5) #Initial time inducer is added
	M.create_rule("assignment", {"equation":"I = _I0*Heaviside(t-_T_I0)"})

	# Add lineage-specific features.
	M.create_volume_event("linear volume", {"growth_rate":1}, 
						  "massaction", {"k":1.2, "species":""})

	vsplit_options = {
					"default":"binomial",
					"G":"duplicate",
					}
	vsplit = LineageVolumeSplitter(M, options = vsplit_options)
	M.create_division_rule("deltaV", {"threshold":1}, vsplit)

	M.create_death_event("death", {}, "hillpositive", {"k":1, "s1":"I", "n":2, "K":5})



	###############################
	# TEST SERIALIZATION OF MODEL #
	###############################
	timepoints = np.arange(0, 5, 0.01)
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		prepickling_lineage = py_SimulateCellLineage(timepoints, Model = M)
	schnitz_tree = prepickling_lineage.get_schnitzes_by_generation()
	prepickling_results = [s.py_get_dataframe(Model=M) for gen in schnitz_tree \
													  for s in gen]
	

	pickled_M = pickle.dumps(M)
	postpickled_M = pickle.loads(pickled_M)

	# test_utils.set_seed(seed)
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		postpickling_lineage = py_SimulateCellLineage(timepoints, 
													  Model = postpickled_M)
	schnitz_tree = prepickling_lineage.get_schnitzes_by_generation()
	postpickling_results = [s.py_get_dataframe(Model=M) for gen in schnitz_tree\
													  for s in gen]

	assert len(prepickling_results) == len(postpickling_results), \
    	"LineageModel simulation produces different number of schnitzes "+\
    	"after pickling."
	for i in range(len(prepickling_results)):
		pd.testing.assert_frame_equal(prepickling_results[i], 
									  postpickling_results[i],
									  check_exact=False)
