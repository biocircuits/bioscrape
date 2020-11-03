import pytest
import test_utils
from bioscrape.simulator import *
from bioscrape.types import *
import bioscrape
import numpy as np

#The same random seeds will be used to compare simulation results that should be identical
seed = 1234
timepoints = np.arange(0, 100.0, .01)


# Test ModelCSimInterface and SafeModelCSimInterface relative to eachother
# Their results should be the same until the non-safe curve goes negative

#Deterministic test
def test_safe_modelsiminterface_deterministic():
    rxn = (["A", "A"], [""], "hillpositive", {"k":1, "K":10, "s1":"B", "n":2})
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
    assert (R_det[-100:, 0] <= 0).all()
    #the safe simulation should stop at 0
    #the non safe simulation should go below 0
    assert (np.round(R_det_s[-100:, 0], 6) == 0).all()

#Stochastic test
def test_safe_modelsiminterface_stochastic():
    rxn = (["A", "A"], [""], "hillpositive", {"k":1, "K":10, "s1":"B", "n":2})
    m = Model(species = ["A", "B"], reactions = [rxn], initial_condition_dict = {"A":10, "B":10})
    sim_stoch = SSASimulator()

    i_fast = ModelCSimInterface(m)
    np.random.seed(seed)
    bioscrape.random.py_seed_random(seed)
    R_stoch_f = sim_stoch.py_simulate(i_fast, timepoints).py_get_result()

    i_safe = SafeModelCSimInterface(m)
    np.random.seed(seed)
    bioscrape.random.py_seed_random(seed)
    R_stoch_s = sim_stoch.py_simulate(i_safe, timepoints).py_get_result()

    #The beginnings of the simulations should be the same
    assert np.allclose(R_stoch_f[:100, :], R_stoch_s[:100, :])
    #The ends of the simulations should be different
    assert not np.allclose(R_stoch_f[100:, :], R_stoch_s[100:, :])
    #the non safe simulation should go below 0
    assert (R_stoch_f[-100:, 0] <= 0).any()
    #the safe simulation should stop at 0
    #the non safe simulation should go below 0
    assert (np.round(R_stoch_s[-100:, 0], 6) == 0).any()

#Stochastic with Volume test
def test_safe_modelsiminterface_stochastic_volume():
    rxn = (["A", "A"], [""], "hillpositive", {"k":1, "K":10, "s1":"B", "n":2})
    m = Model(species = ["A", "B"], reactions = [rxn], initial_condition_dict = {"A":10, "B":10})
    v = Volume()
    v.py_set_volume(1.5)
    sim_stoch_v = VolumeSSASimulator()

    i_fast = ModelCSimInterface(m)
    np.random.seed(seed)
    bioscrape.random.py_seed_random(seed)
    R_stochv_f = sim_stoch_v.py_volume_simulate(i_fast, v, timepoints).py_get_result()

    i_safe = SafeModelCSimInterface(m)
    np.random.seed(seed)
    bioscrape.random.py_seed_random(seed)
    R_stochv_s = sim_stoch_v.py_volume_simulate(i_safe, v, timepoints).py_get_result()


    #The beginnings of the simulations should be the same
    assert np.allclose(R_stochv_f[:100, :], R_stochv_s[:100, :])
    #The ends of the simulations should be different
    assert not np.allclose(R_stochv_f[100:, :], R_stochv_s[100:, :])
    #the non safe simulation should go below 0
    assert (R_stochv_f[-100:, 0] <= 0).any()
    #the safe simulation should stop at 0
    #the non safe simulation should go below 0
    assert (np.round(R_stochv_s[-100:, 0], 6) == 0).any()


    