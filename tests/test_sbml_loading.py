from bioscrape.simulator import *
from bioscrape.types import *
import numpy as np
import pytest
import os

def test_sbml_with_and_without_annotations():
    #This test loads 4 SBML files which all encode the same model
    #all these files should output the same results when simulated
    timepoints = np.arange(0, 10, .1)

    model_path    = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs")

    #an annotated SBML file with local parameters
    #sbml_annotated_local = model_path+"\\models\\sbml_annotated_local.xml"
    sbml_annotated_local = os.path.join(model_path,"models", "sbml_annotated_local.xml")
    CRN1 = Model(sbml_filename = sbml_annotated_local, sbml_warnings = False)
    R1 = py_simulate_model(Model = CRN1, timepoints = timepoints)
    #Check that the simulation results contain all relevant species
    assert "X" in R1
    assert "Y" in R1
    assert "Z" in R1

    #an unnotated SBML file with local parameters
    #sbml_local = model_path+"\\models\\sbml_local.xml"
    sbml_local =  os.path.join(model_path,"models", "sbml_local.xml")

    CRN2 = Model(sbml_filename = sbml_local, sbml_warnings = False)
    R2 = py_simulate_model(Model = CRN2, timepoints = timepoints)
    #Check that the simulation results contain all relevant species
    assert "X" in R2
    assert "Y" in R2
    assert "Z" in R2

    #an annotated SBML file with global parameters

    
    #sbml_annotated_global= model_path+"\\models\\sbml_annotated_global.xml"
    sbml_annotated_global = os.path.join(model_path,"models", "sbml_annotated_global.xml")
    CRN3 = Model(sbml_filename = sbml_annotated_global)
    R3 = py_simulate_model(Model = CRN3, timepoints = timepoints)
    #Check that the simulation results contain all relevant species
    assert "X" in R3
    assert "Y" in R3
    assert "Z" in R3

    #an unnotated SBML file with global parameters
    sbml_global = os.path.join(model_path,"models", "sbml_global.xml")
    CRN4 = Model(sbml_filename = sbml_global, sbml_warnings = False)
    R4 = py_simulate_model(Model = CRN4, timepoints = timepoints)
    #Check that the simulation results contain all relevant species
    assert "X" in R4
    assert "Y" in R4
    assert "Z" in R4

    #Compare that all the outputs are equal
    assert all(R1["X"] == R2["X"])
    assert all(R1["Y"] == R2["Y"])
    assert all(R1["Z"] == R2["Z"])
    assert set(R1.columns) == set(R2.columns)

    assert all(R1["X"] == R3["X"])
    assert all(R1["Y"] == R3["Y"])
    assert all(R1["Z"] == R3["Z"])
    assert set(R1.columns) == set(R3.columns)

    assert all(R3["X"] == R4["X"])
    assert all(R3["Y"] == R4["Y"])
    assert all(R3["Z"] == R4["Z"])
    assert set(R3.columns) == set(R4.columns)


def test_sbml_resaving_local_params_to_global_with_bs_annotation():
    #loads an SBML file with local parameters with bioscrape annotations
    #resaves it via bioscrape, reloads it, simulates, and comapares results.    
    timepoints = np.arange(0, 10, .1)
    model_path = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs")

    #sbml_annotated_local = model_path+"\\models\\sbml_annotated_local.xml"
    sbml_annotated_local = os.path.join(model_path, "models", "sbml_annotated_local.xml")
    CRN1 = Model(sbml_filename = sbml_annotated_local, sbml_warnings = False)
    R1 = py_simulate_model(Model = CRN1, timepoints = timepoints)

    #sbml_annotated_local_bs = model_path+"\\models\\sbml_annotated_local_bs.xml"
    sbml_annotated_local_bs = os.path.join(model_path, "models", "sbml_annotated_local_bs.xml")
    CRN1.write_sbml_model(sbml_annotated_local_bs)
    CRN1l = Model(sbml_filename = sbml_annotated_local_bs, sbml_warnings = False)
    R1l = py_simulate_model(Model = CRN1l, timepoints = timepoints)

    #Check that the saved and loaded models match
    assert all(R1["X"] == R1l["X"])
    assert all(R1["Y"] == R1l["Y"])
    assert all(R1["Z"] == R1l["Z"])
    assert set(R1.columns) == set(R1l.columns)



def test_sbml_resaving_global_params_with_bs_annotation():
    #loads an SBML file with global parameters and bioscrape annotations
    #resaves it via bioscrape, reloads it, simulates, and comapares results.

    timepoints = np.arange(0, 10, .1)
    model_path = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs")
    #an annotated SBML file with global parameters

    #sbml_annotated_global= model_path+"\\models\\sbml_annotated_global.xml"
    sbml_annotated_global = os.path.join(model_path,"models", "sbml_annotated_global.xml")
    CRN3 = Model(sbml_filename = sbml_annotated_global)
    R3 = py_simulate_model(Model = CRN3, timepoints = timepoints)

    #sbml_annotated_global_bs = model_path+"\\models\\sbml_annotated_global_bs.xml"
    sbml_annotated_global_bs = os.path.join(model_path,"models", "sbml_annotated_global_bs.xml")
    CRN3.write_sbml_model(sbml_annotated_global_bs)
    CRN3l = Model(sbml_filename = sbml_annotated_global_bs, sbml_warnings = False)
    R3l = py_simulate_model(Model = CRN3l, timepoints = timepoints)
    #Check that the saved and loaded models match
    assert all(R3["X"] == R3l["X"])
    assert all(R3["Y"] == R3l["Y"])
    assert all(R3["Z"] == R3l["Z"])
    assert set(R3.columns) == set(R3l.columns)



#THE FOLLOWING TWO TESTS HAVE BEEN TURNED OFF DUE TO AN UNRESOLVED PROBLEM
#The issue appears to be related to inconsistencies in how sympy is parsing
#general SBML propensities. In particular, in SBML written by bioscrape,
#Sympy is removing "_" from the start of parameters, and treating these
#strings as species. The origin of this problem remains elusive, but one fix
#is to remove the "_name" convention for parameters.
def test_sbml_resaving_local_params_to_global_no_bs_annotation():
    #loads an SBML file with local parameters and no bioscrape annotations
    #resaves it via bioscrape, reloads it, simulates, and comapares results.
    
    timepoints = np.arange(0, 10, .1)
    model_path = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs")
    #an unnotated SBML file with local parameters
    #sbml_local = model_path+"\\models\\sbml_local.xml"
    sbml_local = os.path.join(model_path,"models", "sbml_local.xml")
    CRN2 = Model(sbml_filename = sbml_local, sbml_warnings = False)
    R2 = py_simulate_model(Model = CRN2, timepoints = timepoints)

    #sbml_local_bs = model_path+"\\models\\sbml_local_bs.xml"
    sbml_local_bs = os.path.join(model_path,"models", "sbml_local_bs.xml")
    CRN2.write_sbml_model(sbml_local_bs)
    CRN2l = Model(sbml_filename = sbml_local_bs, sbml_warnings = False)
    R2l = py_simulate_model(Model = CRN2l, timepoints = timepoints)
    #Check that the saved and loaded models match
    assert set(R2.columns) == set(R2l.columns)
    assert all(R2["X"] == R2l["X"])
    assert all(R2["Y"] == R2l["Y"])
    assert all(R2["Z"] == R2l["Z"])
    

def test_sbml_resaving_global_params_no_bs_annotation():
    #loads an SBML file with global parameters and no bioscrape annotations
    #resaves it via bioscrape, reloads it, simulates, and comapares results.
    
    timepoints = np.arange(0, 10, .1)
    model_path = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs")

    #an unnotated SBML file with global parameters
    #sbml_global = model_path+"\\models\\sbml_global.xml"
    sbml_global = os.path.join(model_path,"models", "sbml_global.xml")

    CRN4 = Model(sbml_filename = sbml_global, sbml_warnings = False)
    R4 = py_simulate_model(Model = CRN4, timepoints = timepoints)

    #sbml_global_bs = model_path+"\\models\\sbml_global_bs.xml"
    sbml_global_bs = os.path.join(model_path,"models", "sbml_global_bs.xml")
    CRN4.write_sbml_model(sbml_global_bs)
    CRN4l = Model(sbml_filename = sbml_global_bs, sbml_warnings = False)
    R4l = py_simulate_model(Model = CRN4l, timepoints = timepoints)
    #Check that the saved and loaded models match
    assert set(R4.columns) == set(R4l.columns)
    assert all(R4["X"] == R4l["X"])
    assert all(R4["Y"] == R4l["Y"])
    assert all(R4["Z"] == R4l["Z"])

def test_delay_annotation():
    """Tests reaction annotation in SBML file for Bioscrape delays, and load accordingly.
    """
    model_path = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs")
    timepoints = np.arange(0, 10, .1)
    species = ["G", "T", "X", "I", "X_m"]
    params = [("ktx", 1.5), ("ktl", 10.0), ("KI", 10), ("n", 2.0), ("KR", 20), ("delta", .1)]
    rxn1d = (["G"], ["G"], "proportionalhillpositive", {"d":"G", "s1":"I", "k":"ktx", "K":"KI", "n":"n"},
       "gaussian", [], ["T"], {"mean":10.0, "std":1.0})
    rxn2d = (["T"], ["T"], "hillpositive", {"s1":"T", "k":"ktl", "K":"KR", "n":1},
            "gamma", [], ["X"], {"k":10.0, "theta":3.0})
    rxn3d = (["X"], [], "massaction", {"k":0.1},
            "fixed", [], ["X_m"], {"delay":10.0})
    rxns_delay = [rxn1d, rxn2d, rxn3d]
    M_delay = Model(species = species, parameters = params, reactions = rxns_delay)
    # Uncomment to create new delay_model in frozeon_sbml_outputs folder.
    M_delay.write_sbml_model(os.path.join(model_path, "models", "delay_model.xml"))
    sbml_delay = os.path.join(model_path, "models", "delay_model.xml")
    CRN1 = Model(sbml_filename = sbml_delay, sbml_warnings = False, input_printout = False)
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(sbml_delay)
    sbml_model = doc.getModel()
    annotation_rxn1 = "<annotation>\n<BioscrapeAnnotation>\n"\
                      "<PropensityType> type=proportionalhillpositive k=ktx K=KI n=n s1=I d=G</PropensityType>\n"\
                      "<DelayType> type=gaussian reactants= products=T "\
                      "mean=DummyVar_GaussianDelay_mean_0 std=DummyVar_GaussianDelay_std_1</DelayType>"\
                      "\n</BioscrapeAnnotation>\n</annotation>"
    annotation_rxn2 = "<annotation>\n<BioscrapeAnnotation>\n"\
                      "<PropensityType> type=hillpositive k=ktl K=KR n=DummyVar_PositiveHillPropensity_n_2 s1=T</PropensityType>\n"\
                      "<DelayType> type=gamma reactants= products=X "\
                      "k=DummyVar_GammaDelay_k_3 theta=DummyVar_GammaDelay_theta_4</DelayType>"\
                      "\n</BioscrapeAnnotation>\n</annotation>"
    annotation_rxn3 = "<annotation>\n<BioscrapeAnnotation>\n"\
                      "<PropensityType> type=massaction k=DummyVar_UnimolecularPropensity_k_5</PropensityType>\n"\
                      "<DelayType> type=fixed reactants= products=X_m delay=DummyVar_FixedDelay_delay_6</DelayType>"\
                      "\n</BioscrapeAnnotation>\n</annotation>"
    assert sbml_model.getReaction(0).getAnnotationString().replace(" ", "") == annotation_rxn1.replace(" ", "")
    assert sbml_model.getReaction(1).getAnnotationString().replace(" ", "") == annotation_rxn2.replace(" ", "")
    assert sbml_model.getReaction(2).getAnnotationString().replace(" ", "") == annotation_rxn3.replace(" ", "")
    R1 = py_simulate_model(Model = CRN1, timepoints = timepoints)
    R2 = py_simulate_model(Model = M_delay, timepoints = timepoints)
    #Compare that all the outputs are equal
    assert all(R1["G"] == R2["G"])
    assert all(R1["T"] == R2["T"])
    assert all(R1["X"] == R2["X"])
    assert all(R1["X_m"] == R2["X_m"])
    assert all(R1["I"] == R2["I"])
    assert set(R1.columns) == set(R2.columns)