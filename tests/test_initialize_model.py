from bioscrape.simulator import *
from bioscrape.types import *
import numpy as np
import pytest

timepoints = np.arange(0, 100, .1)

#This test confirms that calling model.py_initialize twice before simulation does not cause bugs
#  This test was written due to a bug found where the c_propensities vector was not cleared between initializations.
def test_initialize():
    rxn1 = (["A", "A"], ["B"], "massaction", {"k":1.0})
    rxn2 = (["B"], [], "massaction", {"k":.1})
    rxn3 = (["A", "B"], ["C"], "massaction", {"k":2.0})

    #create a model which will be automatically initialized
    M1 = Model(species = ["A", "B", "C"], reactions = [rxn1, rxn2, rxn3], initial_condition_dict = {"A":100})
    R1 = py_simulate_model(Model = M1, timepoints = timepoints)

    #Create another model which will be initialied twice
    M2 = Model(species = ["A", "B", "C"], reactions = [rxn2, rxn3], initial_condition_dict = {"A":100})
    M2.create_reaction(rxn1[0], rxn1[1], rxn1[2], rxn1[3]) #this will uninitialize the model
    #Simulating it will cause it to initialize again
    R2 = py_simulate_model(Model = M2, timepoints = timepoints)

    #the models are identical, so R1 and R2 should have the same output
    assert np.allclose(R1["A"], R2["A"])
    assert np.allclose(R1["B"], R2["B"])
    assert np.allclose(R1["C"], R2["C"])

