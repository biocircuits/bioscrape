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

# Seed RNG value. All tests use this value.
seed = 54173

all_delay_types = ["fixed", "gaussian", "gamma"]
delay_classes   = {
                    "no": Delay,
                    "fixed": FixedDelay,
                    "gaussian": GaussianDelay,
                    "gamma": GammaDelay
                  }
delay_required_params = {
                    "fixed":["delay"], 
                    "gaussian":["mean", "std"], 
                    "gamma":["k", "theta"]
                    }

param_min = -4
param_max = 4

TEST_NAME = "random_delays"

def random_delay_model(delay_type):
    '''
    Returns a randomish model with a non-delayed reaction and a delay reaction. 
    Set to always return the same model, for any particular delay type. 

    WARNING: To produce consistent Models, this function resets the random seeds
    used during Model construction. This may have unexpected effects on random
    number generation outside this function as a side-effect.
    '''

    test_utils.set_seed(seed)

    # Will always consider the reactions A-->B and A+B-->C, where only the 
    # first reaction has delay.
    all_species = ["A", "B", "C"]
    x0 = {"A": 25, "B": 5, "C": 0}

    # A+B-->C
    consol_k = round(np.exp(np.random.uniform(low = param_min,
                                              high = param_max)), 3)
    consolidation_rxn = (["A", "B"], ["C"], "massaction", {'k': consol_k})

    # A-->B reaction, with delay.
    delay_class = delay_classes[delay_type]
    delay_k = round(np.exp(np.random.uniform(low = -1,
                                             high = 1)), 3)
    delay_params = dict()
    for p in delay_required_params[delay_type]:
        delay_params[p] = round(np.exp(np.random.uniform(low = param_min, 
                                                         high = param_max)), 3)
    conversion_rxn = (["A"], [], "massaction", {"k": delay_k}, delay_type, 
                      [], ["B"], delay_params)

    M = Model(reactions = [consolidation_rxn, conversion_rxn])
    M.set_species(x0)

    return M

# print("Simulating Delay =", delay_type, "params=", delay_params)
# results_s = py_simulate_model(timepoints, Model = M, stochastic = True, delay = True)
# plt.plot(timepoints, results_s["C"], label = "stochastic "+str(delay_type)+"params = "+str(delay_params), color = color_list[ind])

@pytest.mark.parametrize('delay_type', all_delay_types)
def test_delay_model(delay_type):
    test_results = dict()
    model = random_delay_model(delay_type)

    timepoints = np.arange(0, 50, 0.01)

    results_s = py_simulate_model(timepoints, Model = model, stochastic = True, return_dataframe = False).py_get_result()

    test_results[delay_type + "_stochastic"]    = results_s

    test_utils.check_sim_results(TEST_NAME, test_results)


## Uncomment this when SBML-writing has been implemented for delay reactions.

# @pytest.mark.parametrize('delay_type', all_delay_types)
# def test_random_propensity_sbml(delay_type):
#     model_dict = dict()
#     model_dict[delay_type] = random_delay_model(delay_type)

#     test_utils.check_sbml_IO(TEST_NAME, model_dict)