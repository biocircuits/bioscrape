import numpy as np
import os
import warnings

import bioscrape.random


frozen_results_loc = os.path.join(os.path.dirname(__file__), 
                                  "frozen_sim_results")

def set_seed(seed):
    '''
    Sets any RNG seeds relevant to bioscrape.
    '''
    if seed == 0:
        raise Warning("Testing with seed = 0! This will not generate ",
                      "consistent results and tests will not pass.")
    np.random.seed(seed)
    bioscrape.random.py_seed_random(seed)

def check_sim_results(test_name, results_dict):
    '''
    Checks whether a dictionary of simulation results matches a "frozen" version
    of those results, within numeric error. 

    If this test has not been run yet (i.e., there is no folder associated with 
    the test), then whatever results are contained in results_dict will be saved 
    as the new frozen output of that test. In this case, the test will always 
    fail, and this function will raise a Warning.

    Parameters:
        test_name - a string descriptor of the test. This will be used to name
                    the folder of frozen results.
        results_dict - a dictionary mapping names of results (usually a string)
                        to numpy arrays of simulation results. 

    Returns: None

    Effects: 
             If frozen results do not already exist, creates a directory of 
                results freezing the data in results_dict, arranged as follows:

            |- frozen_results_loc
            |  - test_name
            |    - results_dict<key#1>.npy
            |    - results_dict<key#2>.npy
            |    ...

            Finally, each result for which there was already a frozen result 
                is asserted to be equal to that result.
    '''
    test_folder = os.path.join(frozen_results_loc, test_name)
    if not os.path.exists(test_folder):
        os.mkdir(test_folder)
        warnings.warn(Warning("No result folder for test ", test_name, 
                              "; creating."))

    all_exist = True

    for sim_name, sim_data in results_dict.items():
        result_file = os.path.join(frozen_results_loc, 
                                   test_name,sim_name + ".npy")
    
        if not os.path.exists(result_file):
            np.save(result_file, sim_data)
            warnings.warn(Warning("Simulation " + test_name +  ":" +  sim_name + 
                          " has no saved result; freezing this result."))
            continue

        frozen_data = np.load(result_file)
        assert np.allclose(sim_data, frozen_data), test_name + ":" + sim_name +\
                                                " doesn't match frozen results"

