import numpy as np
import os
import warnings

import bioscrape.random
import bioscrape.sbmlutil


frozen_results_loc = os.path.join(os.path.dirname(__file__),
                                  "frozen_sim_results")
frozen_sbml_loc    = os.path.join(os.path.dirname(__file__),
                                  "frozen_sbml_outputs")

def set_seed(seed):
    '''
    Sets any RNG seeds relevant to bioscrape.
    '''
    if seed == 0:
        raise Warning("Testing with seed = 0! This will not generate ",
                      "consistent results and tests will not pass.")
    np.random.seed(seed)
    bioscrape.random.py_seed_random(seed)

def check_sbml_IO(test_name, model_dict):
    '''
    Checks that a model's SBML output is the same as a "frozen" version, and
    that the model can be loaded correctly.

    If this test has not been run yet (i.e., there is no folder associated with
    the test), then whatever models are contained in model_dict will be saved
    as SBML as the new frozen output of that test. In this case, the test SHOULD
    always succeed, but this function will raise a Warning.

    Parmeters:
        test_name - a string descriptor of the test. This will be used to name
                    the folder of frozen SBML outputs.
        model_dict - a dictionary mapping names of models (usually a string)
                        to Model objects.
    Returns: None

    Effects:
             If frozen SBML does not already exist, creates a directory of
                results freezing the SBML outputs of the models in model_dict,
                arranged as follows:

            |- frozen_sbml_loc
            |  - test_name
            |    - model_dict<key#1>.sbml
            |    - model_dict<key#2>.sbml
            |    ...

            Finally, each model for which there was already a frozen SBML output
                writes to a temporary SBML file, that file is checked against
                the frozen version, and that temporary file is loaded again
                and the resulting model is asserted to be equal to the original.
    '''
    test_folder = os.path.join(frozen_sbml_loc, test_name)
    if not os.path.exists(test_folder):
        os.mkdir(test_folder)
        warnings.warn(Warning("No SBML output folder for test ", test_name,
                              "; creating."))

    all_exist = True

    for model_name, model in model_dict.items():
        frozen_sbml_file = os.path.join(test_folder, model_name + ".sbml")
        temp_sbml_file   = os.path.join(test_folder, model_name + ".sbml.tmp")

        if not os.path.exists(frozen_sbml_file):
    #         # SAVE MODEL AS SBML: Currently not implemented
            model.write_sbml_model(frozen_sbml_file)
            warnings.warn(Warning("Simulation " + test_name +  ":" + model_name+
                          " has no saved result; freezing this result."))
            continue

        # Test SBML writing.
        # Note: As written, this test is extremely sensitive to formatting.
        # Should probably make it so it ignores e.g. whitespace.
        model.write_sbml_model(temp_sbml_file)
        with open(frozen_sbml_file, 'r') as frozen_text:
            with open(temp_sbml_file, 'r') as temp_text:
                for frozen_line in frozen_text:
                    temp_line = next(temp_text)
                    assert temp_line == frozen_line, \
                            f"{test_name}:{model_name} Model's SBML write " + \
                             "does not match frozen results."
        #We do not know how to test model equality.
        #reloaded_model = bioscrape.sbmlutil.import_sbml(temp_sbml_file)
        #assert reloaded_model == model, f"{test_name}:{model_name} changes " + \
        #                                 "when saved as SBML and reloaded."

def check_sim_results(test_name, results_dict):
    '''
    Checks whether a dictionary of simulation results matches a "frozen" version
    of those results, within numeric error.

    If this test has not been run yet (i.e., there is no folder associated with
    the test), then whatever results are contained in results_dict will be saved
    as the new frozen output of that test. In this case, the test SHOULD always
    succeed, but this function will raise a Warning.

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
                                   test_name, sim_name + ".npy")

        if not os.path.exists(result_file):
            np.save(result_file, sim_data)
            warnings.warn(Warning("Simulation " + test_name +  ":" +  sim_name +
                          " has no saved result; freezing this result."))
            continue

        frozen_data = np.load(result_file, allow_pickle=False)
        assert np.allclose(sim_data, frozen_data), test_name + ":" + sim_name +" doesn't match frozen results"
