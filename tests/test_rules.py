from bioscrape.simulator import *
from bioscrape.types import *
import numpy as np
import pytest
import os 

timepoints = np.arange(0, 10, .01)
model_path = os.path.join(os.path.dirname(__file__), "frozen_sbml_outputs")

def test_dt_assignment_rule():
    #Tests that this rule is applied to both parameters (B) and species (A) every dt
    species = ["A"]
    M = Model(species = species, parameters = {"B":0})
    M.create_rule("assignment", {'equation':"B=B + 1"}, rule_frequency = "dt")
    M.create_rule("assignment", {'equation':"A=B"}, rule_frequency = "dt")

    #Test stochastic simulation
    R = py_simulate_model(timepoints, M, stochastic = True)
    assert np.allclose(R["A"].to_numpy(), np.array([i+1 for i in range(len(timepoints))]))

    #Test deterministic simulation
    M = Model(species = species, parameters = {"B":0})
    M.create_rule("assignment", {'equation':"B=B + 1"}, rule_frequency = "dt")
    M.create_rule("assignment", {'equation':"A=B"}, rule_frequency = "dt")
    R1 = py_simulate_model(timepoints, M, stochastic = False)
    assert np.allclose(R1["A"].to_numpy(), np.array([i+1 for i in range(len(timepoints))]))
    # SBML test
    M.write_sbml_model(os.path.join(model_path, "models", "rule_model_dt.xml"))
    sbml_rule = os.path.join(model_path, "models", "rule_model_dt.xml")
    CRN1 = Model(sbml_filename = sbml_rule, sbml_warnings = False, input_printout = False)
    for r1,r2 in zip(M.get_rules(), CRN1.get_rules()):
        for i,j in zip(r1, r2):
            assert i==j

def test_repeat_assignment_rule():
    #Tests that this rule is applied to both parameters (B) and species (A) every timestep and everytime a reaction fires
    species = ["A", "X"]
    reactions = [([], ["X"], "massaction", {"k":10.0})]
    M = Model(species = species, reactions = reactions, parameters = {"B":0})

    M.create_rule("assignment", {'equation':"B=B + 1"}, rule_frequency = "repeated")
    M.create_rule("assignment", {'equation':"A=B"}, rule_frequency = "repeated")
    R1 = py_simulate_model(timepoints, M, stochastic = True)

    #Every time a reaction occurs or a timestep happens, the rule should be applied
    assert np.allclose(R1["A"].to_numpy(), R1["X"].to_numpy() + np.array([i+1 for i in range(len(timepoints))]))

    #NOTE: unclear what the correct value is for deterministic simulation because the RHS might be called in order to do bounds checking,
    #which can effect the values in a black-box kind of way.
    # SBML test
    M.write_sbml_model(os.path.join(model_path, "models", "rule_model_repeated.xml"))
    sbml_rule = os.path.join(model_path, "models", "rule_model_repeated.xml")
    CRN1 = Model(sbml_filename = sbml_rule, sbml_warnings = False, input_printout = False)
    for r1,r2 in zip(M.get_rules(), CRN1.get_rules()):
        for i,j in zip(r1, r2):
            assert i==j

def test_ode_rule():
    #tests a rule of the form A = A + f(x)dt
    timepoints = np.arange(0, .01, .0000001) #Very small timesteps are required for Eulers method to be accurate

    species = ["A"]
    M = Model(species = species, parameters = {"B":0, "k":.1}, debug = True)
    M.create_rule("ode", {'equation':"k", "target":"B"}, rule_frequency = "dt")
    M.create_rule("ode", {'equation':"B", "target":"A"}, rule_frequency = "dt")
    R1 = py_simulate_model(timepoints, M, stochastic = True)

    quadratic = np.array([.05*i**2 for i in timepoints])
    assert np.allclose(R1["A"].to_numpy(), quadratic)
    # SBML test
    ### Uncomment the following lines to run this test after fixing SBML writing of ODE rule.
    # M.write_sbml_model(os.path.join(model_path, "models", "rule_model_ode.xml"))
    # sbml_rule = os.path.join(model_path, "models", "rule_model_ode.xml")
    # CRN1 = Model(sbml_filename = sbml_rule, sbml_warnings = False, input_printout = False)
    # for r1,r2 in zip(M.get_rules(), CRN1.get_rules()):
    #     for i,j in zip(r1, r2):
    #         assert i==j

def test_time_rule():
    #tests rules that fire at specific times
    species = ["A", "X"]
    reactions = [([], ["X"], "massaction", {"k":10.0})]
    M = Model(species = species, reactions = reactions, parameters = {"B":0})

    #A is assigned to B so it is output
    M.create_rule("assignment", {'equation':f"A=B"}, rule_frequency = "repeated")

    #B will switch between 0 and 1 based on the parity of the time
    for i in range(10):
        val = i % 2
        M.create_rule("assignment", {'equation':f"B={val}"}, rule_frequency = f"{i}")

    R1 = py_simulate_model(timepoints, M, stochastic = True)
    
    assert np.all(R1["A"].to_numpy() == np.array([int(r%2) for r in R1["A"]]))
    # SBML test
    M.write_sbml_model(os.path.join(model_path, "models", "rule_model_time.xml"))
    sbml_rule = os.path.join(model_path, "models", "rule_model_time.xml")
    CRN1 = Model(sbml_filename = sbml_rule, sbml_warnings = False, input_printout = False)
    for r1,r2 in zip(M.get_rules(), CRN1.get_rules()):
        for i,j in zip(r1, r2):
            assert i==j
