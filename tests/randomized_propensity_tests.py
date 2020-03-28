import numpy as np
import pylab as plt

from bioscrape.simulator import *
from bioscrape.types import *


#Parameter ranges to randomly choose parameters (on a log scale)
param_min = -4
param_max = 4

#Names of different supported propensities
propensity_types = ['hillpositive', 'proportionalhillpositive', 'hillnegative', 'proportionalhillnegative', 'massaction', 'general']

#parameter names required for each propensity (general will be treated by itself)
propensity_param_requirements = {
	'massaction':['k'], 'hillpositive':['k', 'K', 'n'], 'hillnegative':['k', 'K', 'n'],
	'proportionalhillpositive':["k", "K", "n"], 'proportionalhillnegative':["k", "K", "n"]
}
#species (passed in as parameters) requires for each propensity (general will be treated by itself)
propensity_specie_requirements = {
	'hillpositive':['s1'], 'hillnegative':['s1'], 'proportionalhillpositive':['s1', 'd'], 'proportionalhillnegative':['s1', 'd'], "massaction":[]
}

#Will always consider the reaction: A+B-->C
inputs = ["A", "B"]
outputs = ["C"]
timepoints = np.arange(0, 10, .01)
x0 = {"A":25, "B": 25, "C":0}

plt.figure()
color_list = ['blue', 'cyan', 'red', 'orange', 'purple', 'green', 'yellow']
for ind in range(len(propensity_types)):
	prop_type = propensity_types[ind]
	print('simulating propensity type', prop_type)

	param_dict = {}
	if prop_type != 'general':
		required_params = propensity_param_requirements[prop_type]
		required_species = propensity_specie_requirements[prop_type]
		param_dict = {}
		for p in required_params:
			param_dict[p] = round(np.exp(np.random.uniform(low = param_min, high = param_max)), 3)
		for i in range(len(required_species)):
			k = required_species[i]
			param_dict[k] = inputs[i]

	elif prop_type == 'general': #Here we will use a random(ish) rational function
		rate_str = "(1+"
		numerator_terms = np.random.randint(0, 5)
		denominator_terms = np.random.randint(0, 5)
		for i in range(numerator_terms):
			coef = str(round(np.exp(np.random.uniform(low = param_min, high = param_max)), 3))
			exp = str(round(np.random.uniform(low = 0, high = param_max), 3))
			specie = (inputs+outputs)[np.random.randint(0, len(inputs)+len(outputs))]
			rate_str += coef+"*"+specie+"^"+exp+"+"
		rate_str =  rate_str[:-1] + ")"
		rate_str += "/(1+"
		for i in range(denominator_terms):
			coef =str(round(np.exp(np.random.uniform(low = param_min, high = param_max)), 3))
			exp = str(round(np.random.uniform(low = 0, high = param_max), 3))
			specie = (inputs+outputs)[np.random.randint(0, len(inputs)+len(outputs))]
			rate_str += coef+"*"+specie+"^"+exp+"+"
		rate_str =  rate_str[:-1] + ")"
		param_dict['rate'] = rate_str
	print('\t params =', param_dict)

	rxn = (inputs, outputs, prop_type, param_dict)
	M = Model(reactions = [rxn])
	M.set_species(x0)
	results_d = py_simulate_model(timepoints, Model = M)
	results_s = py_simulate_model(timepoints, Model = M, stochastic = True)

	plt.plot(timepoints, results_d["C"], label = "deterministic "+str(prop_type)+"params = "+str(param_dict), color = color_list[ind])
	plt.plot(timepoints, results_s["C"], ":", label = "stochastic "+str(prop_type)+"params = "+str(param_dict), color = color_list[ind])
plt.legend()
plt.xlabel("time")
plt.ylabel("C")
plt.show()




