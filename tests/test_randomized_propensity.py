import warnings, os

# We don't want warnings in dependencies to show up in bioscrape's tests.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	import numpy as np
	import pylab as plt
	import random
	import pytest

import test_utils
from bioscrape.simulator import *
from bioscrape.types import *

# Seed RNG value. All tests use this value.
seed = 54173

# Set true to get more diagnostic prints
debug = False

# Parameter ranges to randomly choose parameters (on a log scale)
param_min = -4
param_max = 4

# parameter names required for each propensity (general will be treated by
# itself)
propensity_param_requirements = {
	'massaction':['k'],
	'hillpositive':['k', 'K', 'n'],
	'hillnegative':['k', 'K', 'n'],
	'proportionalhillpositive':["k", "K", "n"],
	'proportionalhillnegative':["k", "K", "n"]
}
# species (passed in as parameters) requires for each propensity (general
# will be treated by itself)
propensity_species_requirements = {
	'hillpositive':['s1'],
	'hillnegative':['s1'],
	'proportionalhillpositive':['s1', 'd'],
	'proportionalhillnegative':['s1', 'd'],
	"massaction":[]
}

all_prop_types = ['hillpositive',
				  'proportionalhillpositive',
				  'hillnegative',
				  'proportionalhillnegative',
				  'massaction', 'general']

TEST_NAME = "random_propensities"


def random_prop_model(prop_type):
	'''
	Returns a randomish model with a specified propensity type. Set to always
	return the same model, for any particular propensity type.

	WARNING: To produce consistent Models, this function resets the random seeds
	used during Model construction. This may have unexpected effects on random
	number generation outside this function as a side-effect.
	'''

	test_utils.set_seed(seed)

	#Will always consider the reaction: A+B-->C
	inputs = ["A", "B"]
	outputs = ["C"]
	all_species = inputs + outputs
	x0 = {"A":25, "B": 25, "C":0}

	if debug:
		print('simulating propensity type ', prop_type)

	param_dict = {}
	# Here we will use a random(ish) rational function
	if prop_type == 'general':
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
	else:
		required_params = propensity_param_requirements[prop_type]
		required_species = propensity_species_requirements[prop_type]
		param_dict = {}
		for p in required_params:
			param_dict[p] = \
					round(np.exp(np.random.uniform(low = param_min,
												   high = param_max)), 3)
		for i in range(len(required_species)):
			k = required_species[i]
			param_dict[k] = inputs[i]
	if debug:
		print('\t params =', param_dict)

	rxn = (inputs, outputs, prop_type, param_dict)
	M = Model(reactions = [rxn], initial_condition_dict = x0)
	M.set_species(x0)

	return M

# def test_debug():
# 	import bioscrape.sbmlutil
# 	bioscrape.sbmlutil.import_sbml("frozen_sbml_outputs/random_propensities/hillnegative.sbml.tmp")

@pytest.mark.parametrize('prop_type', all_prop_types)
def test_random_propensity_outputs(prop_type):
	test_results = dict()
	model = random_prop_model(prop_type)

	timepoints = np.arange(0, 50, .01)

	results_d = py_simulate_model(timepoints, Model = model, stochastic = False, return_dataframe = False).py_get_result()
	results_s = py_simulate_model(timepoints, Model = model, stochastic = True, return_dataframe = False).py_get_result()

	test_results[prop_type + "_deterministic"] = results_d
	test_results[prop_type + "_stochastic"]    = results_s

	test_utils.check_sim_results(TEST_NAME, test_results)

@pytest.mark.parametrize('prop_type', all_prop_types)
def test_random_propensity_sbml(prop_type):
	model_dict = dict()
	model_dict[prop_type] = random_prop_model(prop_type)

	test_utils.check_sbml_IO(TEST_NAME, model_dict)

def debug_random_prop_tests():
	'''
	This is not a test.

	Plot frozen results for debugging purposes.
	'''
	propensity_types = ['hillpositive', 'proportionalhillpositive',
					    'hillnegative', 'proportionalhillnegative',
					    'massaction', 'general']
	colors = {
		'massaction':'blue',
		'hillpositive': 'cyan',
		'hillnegative': 'red',
		'proportionalhillpositive': 'orange',
		'proportionalhillnegative': 'purple',
		'general': 'black'
		}
	test_loc = os.path.join(test_utils.frozen_results_loc, TEST_NAME)

	plt.figure()
	for prop_type in propensity_types:
		results_d = np.load(os.path.join(test_loc, 
										 prop_type + "_deterministic.npy"))
		plt.plot(results_d[:,0], results_d[:,3],
				label = "deterministic "+str(prop_type),
						# +"params = "+str(param_dict),
				color = colors[prop_type])

		results_s = np.load(os.path.join(test_loc, 
										 prop_type + "_stochastic.npy"))
		plt.plot(results_s[:,0], results_s[:,3], ":",
			       label = "stochastic "+str(prop_type),
			       		   # +"params = "+str(param_dict),
			       color = colors[prop_type])

	# plt.legend()
	plt.xlabel("time")
	plt.ylabel("C")
	plt.legend()
	plt.show()