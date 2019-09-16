from bioscrape.types import import_sbml
from bioscrape.simulator import py_simulate_model

import numpy as np
import matplotlib.pyplot as plt

M_sbml = import_sbml('models/sbml_test.xml')
timepoints = np.linspace(0,100,1000)

result = py_simulate_model(timepoints, Model = M_sbml)

plt.figure()
for s in M_sbml.get_species_list():
    plt.plot(timepoints, result[s], label = s)

plt.legend()