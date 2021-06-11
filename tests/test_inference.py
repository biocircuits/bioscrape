
import pytest
import test_utils
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

from bioscrape.types import Model
from bioscrape.simulator import py_simulate_model
from bioscrape.inference import py_inference
from bioscrape.inference_setup import InferenceSetup
from emcee import EnsembleSampler

np.random.seed(123)

def test_safe_modelsiminterface_deterministic():
    # Choose the "true" parameters.
    m_true = -0.9594
    b_true = 4.294
    # f_true = 0.534

    species = ['y']
    parameters = {'m':m_true, 'b': b_true}
    rule = ('assignment',{'equation':'y = _m*t + _b'})
    x0 = {'y':0}
    M = Model(species = species, parameters = parameters, rules = [rule], initial_condition_dict = x0)


    #Simulate the Model deterministically
    x0 = np.linspace(0, 10, 50)
    results_det = py_simulate_model(x0, Model = M) #Returns a Pandas DataFrame

    # Generate some synthetic data from the model.
    N = 50
    x = np.sort(10 * np.random.rand(N))
    yerr = 0.1 + 0.6 * np.random.rand(N)
    y = m_true * x + b_true
    # y += np.abs(f_true * y) * np.random.randn(N)
    y += yerr * np.random.randn(N)

    x0 = np.linspace(0, 10, 50)
    exp_data = pd.DataFrame()
    exp_data['x'] = x0
    exp_data['y'] = y

    prior = {'m' : ['gaussian', m_true, 500],'b' : ['gaussian', b_true, 1000]}
    sampler, pid = py_inference(Model = M, exp_data = exp_data, measurements = ['y'],
                        time_column = ['x'], params_to_estimate = ['m','b'],
                        nwalkers = 32, nsteps = 2000, init_seed = 1e-4, prior = prior,
                        sim_type = 'deterministic', plot_show = False)
    assert(isinstance(sampler, EnsembleSampler) == True)
    assert(isinstance(pid, InferenceSetup) == True)
    assert np.array(sampler.get_autocorr_time())[0] < 50
    assert np.array(sampler.get_autocorr_time())[1] < 50
    assert np.array(sampler.acceptance_fraction).all() < 2