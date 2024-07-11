from bioscrape.simulator import *
from bioscrape.types import *
from bioscrape.analysis import *
import numpy as np
import pytest
import os

@pytest.fixture(scope = "module")
def model_setup():
    """
    Set up the model for sensitivity analysis tests.
    As an example, we take the toggle switch model.
    """
    # Toggle Switch dynamics
    species = ["m_t", "m_l", "p_t", "p_l"]
    params = [("K",100),("b_t",50),("b_l",10),("d_t",5),("d_l",5),
            ("del_t",0.01), ("del_l",0.02), ("beta_t", 0.02), ("beta_l", 0.01)]
    rxn1 = ([], ["m_t"], "hillnegative", {"s1":"p_l", "k":"K", "K":"b_t", "n":2})
    rxn2 = (["m_t"],[], "massaction", {"k":"d_t"})
    rxn3 = ([], ["m_l"], "hillnegative", {"s1":"p_t", "k":"K", "K":"b_l", "n":2})
    rxn4 = (["m_l"],[], "massaction", {"k":"d_l"})

    rxn5 = (["p_t"],[], "massaction", {"k":"del_t"})
    rxn6 = ([], ["p_t"], "general", {"rate":"beta_t*m_t"})

    rxn7 = (["p_l"],[], "massaction", {"k":"del_l"})
    rxn8 = ([], ["p_l"], "general", {"rate":"beta_l*m_l"})

    reactions = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6, rxn7, rxn8]

    x0 = {"m_t":1e-5, "m_l":1e-5, "p_t":1e-5, "p_l":1e-5}
    M = Model(species = species, parameters = params, reactions = reactions, initial_condition_dict = x0)
    timepoints = np.linspace(0,100,100)
    return M, timepoints

def test_jacobian(model_setup):
    """Jacobian test
    """
    M, timepoints = model_setup 
    states = np.array([2, 4, 5, 10])
    jacobian_true = np.array([[-5, 0, 0, -7.39644970e-01],
                             [0, -5, -6.4,  0],
                             [2e-02,  0, -1e-02,  0],
                             [0,  1e-02, 0, -2e-02]])
    jacobian_4thorder = py_get_jacobian(M, states, method = 'fourth_order_central_difference')
    assert np.allclose(jacobian_true, jacobian_4thorder, rtol = 1e-2)
    jacobian_central = py_get_jacobian(M, states, method = 'central_difference')
    assert np.allclose(jacobian_true, jacobian_central, rtol = 1e-2)
    jacobian_backward = py_get_jacobian(M, states, method = 'backward_difference')
    assert np.allclose(jacobian_true, jacobian_backward, rtol = 1e-2)
    jacobian_forward = py_get_jacobian(M, states, method = 'forward_difference')
    assert np.allclose(jacobian_true, jacobian_forward, rtol = 1e-2)

def test_sensitivity_to_parameter(model_setup):
    """Sensitivity to parameter test
    """
    M, _ = model_setup 
    states = np.array([20, 4, 5, 10])
    sensitivity_to_parameter_true = np.array([0,  0,  states[0], 0])
    sensitivity_to_parameter_4thorder = py_get_sensitivity_to_parameter(M, states, 'beta_t', method = 'fourth_order_central_difference')
    assert np.allclose(sensitivity_to_parameter_true, sensitivity_to_parameter_4thorder, rtol = 1e-2)
    sensitivity_to_parameter_central = py_get_sensitivity_to_parameter(M, states, 'beta_t', method = 'central_difference')
    assert np.allclose(sensitivity_to_parameter_true, sensitivity_to_parameter_central, rtol = 1e-2)
    sensitivity_to_parameter_backward = py_get_sensitivity_to_parameter(M, states, 'beta_t', method = 'backward_difference')
    assert np.allclose(sensitivity_to_parameter_true, sensitivity_to_parameter_backward, rtol = 1e-2)
    sensitivity_to_parameter_forward = py_get_sensitivity_to_parameter(M, states, 'beta_t', method = 'forward_difference')
    assert np.allclose(sensitivity_to_parameter_true, sensitivity_to_parameter_forward, rtol = 1e-2)

def test_SSM(model_setup):
    """Test sensitivty coefficient matrix for all time points for all parameters for all states.
    """
    M, timepoints = model_setup
    params_values = M.get_parameter_values()
    n_species = 4
    SSM = py_sensitivity_analysis(M, timepoints, normalize = True)
    assert np.shape(SSM) == (len(timepoints), len(params_values), n_species)
    # Check out the Sensitivity Analysis ipython notebook in bioscrape/examples for more.
