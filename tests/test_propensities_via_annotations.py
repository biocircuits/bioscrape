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


model_path = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs", "models")
timepoints = np.arange(0, 10, .1)
seed = 54173

def test_massaction():
    #First model
    sbml_bs_1 = os.path.join(model_path, "massaction_crn_1_bs.xml")
    sbml_general_1 = os.path.join(model_path, "massaction_crn_1.xml")
    CRN_bs_1 = Model(sbml_filename = sbml_bs_1, sbml_warnings = False)
    CRN_general_1 = Model(sbml_filename = sbml_general_1, sbml_warnings = False)

    #Test Deterministic Simulation
    R_bs_1_det = py_simulate_model(Model = CRN_bs_1, timepoints = timepoints)
    R_general_1_det = py_simulate_model(Model = CRN_general_1, timepoints = timepoints)

    #test that the models have the same species
    assert set(R_bs_1_det.columns) == set(R_general_1_det.columns)
    #test that the output of X is the same deterministically
    assert np.allclose(R_bs_1_det["protein_X"], R_general_1_det["protein_X"])

    #Test Stochastic Simulation
    test_utils.set_seed(seed)
    R_bs_1_stoch = py_simulate_model(Model = CRN_bs_1, timepoints = timepoints, stochastic = True)
    test_utils.set_seed(seed)
    R_general_1_stoch = py_simulate_model(Model = CRN_general_1, timepoints = timepoints, stochastic = True)

    #test that the models have the same species
    assert set(R_general_1_stoch.columns) == set(R_bs_1_stoch.columns)
    #test that the output of X is the same deterministically
    assert np.allclose(R_general_1_stoch["protein_X"], R_bs_1_stoch["protein_X"])

    sbml_bs_2 = os.path.join(model_path, "massaction_crn_2_bs.xml")
    sbml_general_2 = os.path.join(model_path, "massaction_crn_2.xml")
    CRN_bs_2 = Model(sbml_filename = sbml_bs_2, sbml_warnings = False)
    CRN_general_2 = Model(sbml_filename = sbml_general_2, sbml_warnings = False)

    #Test Deterministic Simulation
    R_bs_2_det = py_simulate_model(Model = CRN_bs_2, timepoints = timepoints)
    R_general_2_det = py_simulate_model(Model = CRN_general_2, timepoints = timepoints)

    #test that the models have the same species
    assert set(R_bs_2_det.columns) == set(R_general_2_det.columns)
    #test that the output of X is the same deterministically
    assert np.allclose(R_bs_2_det["protein_X"], R_general_2_det["protein_X"])

    #Test Stochastic Simulation
    test_utils.set_seed(seed)
    R_bs_2_stoch = py_simulate_model(Model = CRN_bs_2, timepoints = timepoints, stochastic = True)
    test_utils.set_seed(seed)
    R_general_2_stoch = py_simulate_model(Model = CRN_general_2, timepoints = timepoints, stochastic = True)

    #test that the models have the same species
    assert set(R_general_2_stoch.columns) == set(R_bs_2_stoch.columns)
    #test that the output of X is the same deterministically
    assert np.allclose(R_general_2_stoch["protein_X"], R_bs_2_stoch["protein_X"])

def test_hill_positive():
    sbml_bs = os.path.join(model_path, "hill_positive_crn_bs.xml")
    sbml_general = os.path.join(model_path, "hill_positive_crn.xml")
    CRN_bs = Model(sbml_filename = sbml_bs, sbml_warnings = False)
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

def test_hill_negative():
    sbml_bs = os.path.join(model_path, "hill_negative_crn_bs.xml")
    sbml_general = os.path.join(model_path, "hill_negative_crn.xml")
    CRN_bs = Model(sbml_filename = sbml_bs, sbml_warnings = False)
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

def test_proportional_hill_positive():
    sbml_bs = os.path.join(model_path, "prop_hill_positive_crn_bs.xml")
    sbml_general = os.path.join(model_path, "prop_hill_positive_crn.xml")
    CRN_bs = Model(sbml_filename = sbml_bs, sbml_warnings = False)
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

def test_proportional_hill_negative():
    sbml_bs = os.path.join(model_path, "prop_hill_negative_crn_bs.xml")
    sbml_general = os.path.join(model_path, "prop_hill_negative_crn.xml")
    CRN_bs = Model(sbml_filename = sbml_bs, sbml_warnings = False)
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