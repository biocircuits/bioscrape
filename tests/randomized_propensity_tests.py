import warnings

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


# Set to True to get in-line prints and plotting
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

TEST_NAME = "random_propensities"

@pytest.mark.parametrize('prop_type', 
					   ['hillpositive', 'proportionalhillpositive', 
					    'hillnegative', 'proportionalhillnegative', 
					    'massaction', 'general'])
def test_random_propensities(prop_type):
	# RNG seed
	seed = 54173
	test_utils.set_seed(seed)

	#Will always consider the reaction: A+B-->C
	inputs = ["A", "B"]
	outputs = ["C"]
	all_species = inputs + outputs
	timepoints = np.arange(0, 50, .01)
	x0 = {"A":25, "B": 25, "C":0}

	test_results = dict()

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
	results_d = py_simulate_model(timepoints, Model = M)
	results_s = py_simulate_model(timepoints, Model = M, stochastic = True)

	test_results[prop_type + "_deterministic"] = results_d
	test_results[prop_type + "_stochastic"]    = results_s

	test_utils.check_sim_results(TEST_NAME, test_results)

def debug_random_prop_tests():
	'''
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
		'proportionalhillnegative': 'purple'
		}				   					    
	test_loc = os.path.join(test_utils.frozen_results_loc, TEST_NAME)

	plt.figure()
	for prop_type in propensity_types:

		results_d = np.load(os.path.join(test_loc, prop_type))

		plt.plot(timepoints, results_d["C"], 
				label = "deterministic "+str(prop_type)\
						+"params = "+str(param_dict), 
				color = colors[prop_type])
		plt.plot(timepoints, results_s["C"], ":", 
			       label = "stochastic "+str(prop_type)+"params = "\
			       		   +str(param_dict),
			       color = colors[prop_type])

	# plt.legend()
	plt.xlabel("time")
	plt.ylabel("C")
	plt.show()