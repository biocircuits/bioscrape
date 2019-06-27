import numpy as np
import pylab as plt
from bioscrape.simulator import *
from bioscrape.types import *
import warnings




#Parameter ranges to randomly choose parameters (on a log scale)
param_min = -1
param_max = 1
include_delays = False

#Names of different supported propensities
propensity_types = ['massaction']#, 'hillpositive', 'proportionalhillpositive', 'hillnegative', 'proportionalhillnegative', 'massaction']#, 'general']

#parameter names required for each propensity (general will be treated by itself)
propensity_param_requirements = {
	'massaction':['k'], 'hillpositive':['k', 'K', 'n'], 'hillnegative':['k', 'K', 'n'],
	'proportionalhillpositive':["k", "K", "n"], 'proportionalhillnegative':["k", "K", "n"]
}
#species (passed in as parameters) requires for each propensity (general will be treated by itself)
propensity_specie_requirements = {
	'hillpositive':['s1'], 'hillnegative':['s1'], 'proportionalhillpositive':['s1', 'd'], 'proportionalhillnegative':['s1', 'd'], "massaction":[]
}
delay_types = [None, 'fixed', 'gaussian', 'gamma']
delay_required_params = {
	"fixed":["delay"], "gaussian":["mean", "std"], "gamma":["k", "theta"]
}


species = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
n_species = len(species)

tests = 10000
for test in range(tests):
	n_reactions = np.random.randint(1, 3)

	reactions = []
	for r in range(n_reactions):

		try_again = True
		while try_again:#Loop ensures no positive feedback which leads to long simulations
			inputs = []
			outputs = []
			while(len(inputs) == 0 and len(outputs) == 0):
			
				n_inputs = np.random.randint(0, 5)
				for i in range(n_inputs):
					inputs.append(species[np.random.randint(0, len(species))])

				n_outputs = np.random.randint(0, 5)
				for i in range(n_outputs):
					outputs.append(species[np.random.randint(0, len(species))])

			delay_inputs = []
			delay_outputs = []
			while(len(delay_inputs) == 0 and len(delay_outputs)==0):
				n_delay_inputs = np.random.randint(0, 5)
				
				for i in range(n_delay_inputs):
					delay_inputs.append(species[np.random.randint(0, len(species))])

				n_delay_outputs = np.random.randint(0, 5)
				for i in range(n_delay_outputs):
					delay_outputs.append(species[np.random.randint(0, len(species))])

			inputs_in_outputs = len([i for i in inputs if i in outputs or i in delay_outputs])
			if inputs_in_outputs >= len(inputs):
				try_again = True
			else:
				try_again = False


		prop_type = propensity_types[np.random.randint(0, len(propensity_types))]
		param_dict = {}
		if prop_type != 'general':
			required_params = propensity_param_requirements[prop_type]
			required_species = propensity_specie_requirements[prop_type]
			param_dict = {}
			for p in required_params:
				param_dict[p] = round(np.exp(np.random.uniform(low = param_min, high = param_max)), 3)
			for i in range(len(required_species)):
				k = required_species[i]
				param_dict[k] = species[np.random.randint(0, len(species))]

		elif prop_type == 'general': #Here we will use a random(ish) rational function
			rate_str = "(1+"
			numerator_terms = np.random.randint(0, 5)
			denominator_terms = np.random.randint(0, 5)
			for i in range(numerator_terms):
				coef = str(round(np.exp(np.random.uniform(low = param_min, high = param_max)), 3))
				exp = str(round(np.random.uniform(low = 0, high = param_max), 3))
				specie = species[np.random.randint(0, len(species))]
				rate_str += coef+"*"+specie+"^"+exp+"+"
			rate_str =  rate_str[:-1] + ")"
			rate_str += "/(1+"
			for i in range(denominator_terms):
				coef =str(round(np.exp(np.random.uniform(low = param_min, high = param_max)), 3))
				exp = str(round(np.random.uniform(low = 0, high = param_max), 3))
				specie = species[np.random.randint(0, len(species))]
				rate_str += coef+"*"+specie+"^"+exp+"+"
			rate_str =  rate_str[:-1] + ")"
			param_dict['rate'] = rate_str

		delay_type = delay_types[np.random.randint(0, len(delay_types))]
		if delay_type == None or not include_delays:
			delay_inputs = None
			delay_outputs = None
			delay_params = None
			rxn = (inputs, outputs, prop_type, param_dict)
		else:
			delay_params = {}
			for p in delay_required_params[delay_type]:
				delay_params[p] = round(np.exp(np.random.uniform(low = param_min, high = param_max)), 3)

			rxn = (inputs, outputs, prop_type, param_dict, delay_type, delay_inputs, delay_outputs, delay_params)
		reactions.append(rxn)

	timepoints = np.arange(0, 0.001, .0001)
	print("Simulating model", test, "#rxns=", len(reactions), "rxns=", reactions)
	M = Model(reactions = reactions, initial_condition_dict = {s:np.random.randint(10, 100) for s in species})
	results_s = py_simulate_model(timepoints, Model = M, stochastic = True, delay = include_delays, safe = True)		
	print("simulation successful")
