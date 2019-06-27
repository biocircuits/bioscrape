import numpy as np
import pylab as plt
from bioscrape.simulator import *
from bioscrape.types import *

delay_types = ['fixed', 'gaussian', 'gamma']
delay_required_params = {
	"fixed":["delay"], "gaussian":["mean", "std"], "gamma":["k", "theta"]
}

#Will always consider the reaction: A+B-->C @ rate k=1 (but occurs with a delay)
inputs = ["A", "B"]
outputs = []
delay_inputs = []
delay_outputs = ["C"]
timepoints = np.arange(0, 10, .01)
x0 = {"A":25, "B": 25, "C":0}
prop_type = "massaction"
prop_params = {"k":1.0}

#Parameter ranges to randomly choose parameters (on a log scale)
param_min = -4
param_max = 4


plt.figure()
color_list = ['blue', 'cyan', 'red', 'orange', 'purple', 'green', 'yellow']
for ind in range(len(delay_types)):
	delay_type = delay_types[ind]
	delay_params = {}
	for p in delay_required_params[delay_type]:
		delay_params[p] = round(np.exp(np.random.uniform(low = param_min, high = param_max)), 3)
	rxn = (inputs, outputs, prop_type, prop_params, delay_type, delay_inputs, delay_outputs, delay_params)

	M = Model(reactions = [rxn])
	M.set_species(x0)
	print("Simulating Delay =", delay_type, "params=", delay_params)
	results_s = py_simulate_model(timepoints, Model = M, stochastic = True, delay = True)
	plt.plot(timepoints, results_s["C"], label = "stochastic "+str(delay_type)+"params = "+str(delay_params), color = color_list[ind])

plt.legend()
plt.title("Simulating A+B-->C as a delay reaction with random parameters")
plt.xlabel("time")
plt.ylabel("C")
plt.show()

