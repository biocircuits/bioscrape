from bioscrape.simulator import *
from bioscrape.types import *
import numpy as np
import pytest
import os
import test_utils


#This test file compares all propensity types 
#by loading SBML models with and without bioscrape annotations. 
#This will help ensure that specific propensity types match
#general propensities. Deterministic and SSA without Volume are tested.
#general propensities will not necessarily be volume dependent,
#so volume SSA requires a different set of tests.
#Models were generated via BioCRNpyler
all_prop_types = ['hill_positive',
                  'prop_hill_positive',
                  'hill_negative',
                  'prop_hill_negative',
                  'massaction_1',
                  'massaction_2']

model_path = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs", "models")
timepoints = np.arange(0, 10, .1)
seed = 54173

@pytest.mark.parametrize('prop_type', all_prop_types)
def test_props_via_annotations_all(prop_type):
    sbml_bs = os.path.join(model_path, f"{prop_type}_crn_bs.xml")
    sbml_general = os.path.join(model_path, f"{prop_type}_crn.xml")
    CRN_bs = Model(sbml_filename = sbml_bs, sbml_warnings = False, input_printout = True)
    CRN_general = Model(sbml_filename = sbml_general, sbml_warnings = False)

    #Test Deterministic Simulation
    R_bs_det = py_simulate_model(Model = CRN_bs, timepoints = timepoints)
    R_general_det = py_simulate_model(Model = CRN_general, timepoints = timepoints)

    #test that the models have the same species
    assert set(R_bs_det.columns) == set(R_general_det.columns)
    #test that the output of X is the same deterministically
    assert np.allclose(R_bs_det["protein_X"], R_general_det["protein_X"])

    #Test Stochastic Simulation
    test_utils.set_seed(seed)
    R_bs_stoch = py_simulate_model(Model = CRN_bs, timepoints = timepoints, stochastic = True)
    test_utils.set_seed(seed)
    R_general_stoch = py_simulate_model(Model = CRN_general, timepoints = timepoints, stochastic = True)

    #test that the models have the same species
    assert set(R_general_stoch.columns) == set(R_bs_stoch.columns)
    #test that the output of X is the same deterministically
    assert np.allclose(R_general_stoch["protein_X"], R_bs_stoch["protein_X"])
