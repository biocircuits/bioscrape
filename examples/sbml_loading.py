from bioscrape.types import import_sbml, read_model_from_sbml, Model
from bioscrape.simulator import py_simulate_model

import numpy as np
import matplotlib.pyplot as plt

# Import SBML
M_sbml = import_sbml('models/sbml_test.xml')

# Load SBML using the old method that works with strings - 
# M_sbml = read_model_from_sbml('models/sbml_test.xml')

# Write it to bioscrape
M_sbml.write_bioscrape_xml('models/sbml_test_bioscrape.xml')

timepoints = np.linspace(0,100,1000)

result = py_simulate_model(timepoints, Model = M_sbml)

plt.figure()
for s in M_sbml.get_species_list():
    plt.plot(timepoints, result[s], label = s)

plt.legend()
plt.show()
