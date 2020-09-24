import pytest
import test_utils
from bioscrape.simulator import *
from bioscrape.types import *
import numpy as np

#The same random seeds will be used to compare simulation results that should be identical
seed = 1234
timepoints = np.arange(0, 30.0, .01)


# Test ModelCSimInterface and SafeModelCSimInterface relative to eachother
# Their results should be the same until the curve goes negative
def test_safe_modelsiminterface_deterministic():
    rxn = (["A"], [""], "hillpositive", {"k":1, "K":10, "s1":"B", "n":2})
    m = Model(species = ["A", "B"], reactions = [rxn], initial_condition_dict = {"A":10, "B":10})

    sim_det = DeterministicSimulator()

    i_fast = ModelCSimInterface(m)
    i_fast.py_prep_deterministic_simulation()
    R_det = sim_det.py_simulate(i_fast, timepoints).py_get_result()

    i_safe = SafeModelCSimInterface(m)
    i_safe.py_prep_deterministic_simulation()
    R_det_s = sim_det.py_simulate(i_safe, timepoints).py_get_result()

    #The beginnings of the simulations should be the same
    assert np.allclose(R_det_s[:100, :], R_det[:100, :])
    #The ends of the simulations should be different
    assert not np.allclose(R_det_s[100:, :], R_det[100:, :])
    #the non safe simulation should go below 0
    assert (R_det[-100:, 0] < 0).all()
    #the safe simulation should stop at 0
    #the non safe simulation should go below 0
    assert (R_det_s[-100:, 0] == 0).all()


    