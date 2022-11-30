
import pytest
import test_utils
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd

from bioscrape.types import Model
from bioscrape.simulator import py_simulate_model
from bioscrape.inference import py_inference
from bioscrape.inference_setup import InferenceSetup
from bioscrape.pid_interfaces import PIDInterface
from emcee import EnsembleSampler
from lmfit.minimizer import MinimizerResult

np.random.seed(123)

@pytest.fixture(scope = "module")
def model_setup():
    """
    Set up the model for inference tests
    """
    m_true = -0.9594
    b_true = 4.294
    # f_true = 0.534

    species = ['y']
    parameters = {'m':m_true, 'b': b_true}
    rule = ('assignment',{'equation':'y = m*t + b'})
    x0 = {'y':0}
    M = Model(species = species, parameters = parameters, rules = [rule], initial_condition_dict = x0)
    params_to_estimate = ['m','b']
    return M, params_to_estimate

def test_getstate(model_setup):
    M, params_to_estimate = model_setup

    IS = InferenceSetup(Model = M, params_to_estimate = params_to_estimate)

    state = IS.__getstate__()

    assert type(state[0]) == type(M)
    assert state[1] == params_to_estimate
    assert state[2] == IS.init_seed
    assert state[3] == IS.prior
    assert state[4] == IS.nwalkers
    assert state[5] == IS.nsteps
    assert state[6] == IS.dimension
    assert state[7] == IS.exp_data
    assert state[8] == IS.sim_type
    assert state[9] == IS.method
    assert state[10] == IS.timepoints
    assert state[11] == IS.time_column
    assert state[12] == IS.measurements
    assert state[13] == IS.initial_conditions
    assert state[14] == IS.parameter_conditions
    assert state[15] == IS.norm_order
    assert state[16] == IS.N_simulations
    assert state[17] == IS.LL_data
    assert state[18] == IS.debug
    assert state[19] == IS.cost_progress
    assert state[20] == IS.cost_params
    assert state[21] == IS.hmax

def test_setstate(model_setup):
    M, params_to_estimate = model_setup
    
    IS = InferenceSetup(Model = M, params_to_estimate = params_to_estimate)

    init_conds = {'y':1.0}
    param_conds = {'b':1.0}
    params_to_estimate = ['m']
    exp_data = pd.DataFrame()
    exp_data['y'] = np.linspace(0, 10, 50)
    exp_data['T'] = np.linspace(0, 10, 50)
    timepoints = None

    state = (
        M, params_to_estimate, .11, None, 11, 111, 2, exp_data, 'stochastic', 'lmfit', 
        timepoints, 'T', ['y'], init_conds, None, 1, 11, None, True, None, None, 1.0
        )

    IS.__setstate__(state)

    assert type(state[0]) == type(M)
    assert state[1] == IS.params_to_estimate
    assert state[2] == IS.init_seed
    assert state[3] == IS.prior
    assert state[4] == IS.nwalkers
    assert state[5] == IS.nsteps
    assert state[6] == IS.dimension
    assert 'y' in IS.exp_data and 'T' in IS.exp_data
    assert state[8] == IS.sim_type
    assert state[9] == IS.method
    #print("state[10]", state[10])
    #print("IS.timepoints", IS.timepoints)
    assert (np.linspace(0, 10, 50) == IS.timepoints).all() #this is extracted from the data
    assert state[11] == IS.time_column
    assert state[12] == IS.measurements
    assert state[13] in IS.initial_conditions #this gets multiplied out
    assert state[14] == IS.parameter_conditions
    assert state[15] == IS.norm_order
    assert state[16] == IS.N_simulations
    #assert state[17] == IS.LL_data #SKIP this one as LL.data is extracted
    assert state[18] == IS.debug
    assert state[19] == IS.cost_progress
    assert state[20] == IS.cost_params
    assert state[21] == IS.hmax

def test_basic_inference(model_setup):
    M, _ = model_setup
    m_true = -0.9594
    b_true = 4.294
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
                        nwalkers = 32, nsteps = 2000, discard = 100, init_seed = 1e-4, prior = prior,
                        sim_type = 'deterministic', plot_show = False)
    assert(isinstance(sampler, EnsembleSampler) == True)
    assert(isinstance(pid, InferenceSetup) == True)
    assert np.array(sampler.get_autocorr_time())[0] < 50
    assert np.array(sampler.get_autocorr_time())[1] < 50
    assert np.array(sampler.acceptance_fraction).all() < 2

def test_uniform_priors(model_setup):
    """ Uniform prior testing
    """
    M, params_to_estimate = model_setup
    prior = {'m' : ['uniform', -100, 100],'b' : ['uniform', 0, 1000]}
    test_pid_interface = PIDInterface(params_to_estimate, M, prior)
    params_dict = {'m':-10,'b':1}
    lp = 0
    log_prior = 0
    for param in prior:
        rv = scipy.stats.uniform(loc = prior[param][1], scale = prior[param][2]-prior[param][1])
        log_prior += np.log(rv.pdf(params_dict[param]))
    lp = test_pid_interface.check_prior(params_dict)
    np.testing.assert_allclose(lp, log_prior, rtol = 0.1)

def test_gaussian_priors(model_setup):
    """ Gaussian prior testing
    """
    M, params_to_estimate = model_setup
    prior = {'m' : ['gaussian', -0.9, 500],'b' : ['gaussian', 4.2, 1000]}
    test_pid_interface = PIDInterface(params_to_estimate, M, prior)
    params_dict = {'m':0,'b':1}
    lp = 0
    log_prior = 0
    for param in prior:
        rv = scipy.stats.norm(loc = prior[param][1], scale = prior[param][2])
        log_prior += np.log(rv.pdf(params_dict[param]))
    lp = test_pid_interface.check_prior(params_dict)
    np.testing.assert_allclose(lp, log_prior, rtol = 0.1)

def test_exponential_priors(model_setup):
    """ Exponential prior testing
    """
    M, params_to_estimate = model_setup
    prior = {'m' : ['exponential', 2],'b' : ['exponential', 5]}
    test_pid_interface = PIDInterface(params_to_estimate, M, prior)
    params_dict = {'m':1,'b':2}
    lp = 0
    log_prior = 0
    for param in prior:
        rv = scipy.stats.expon(scale = 1/prior[param][1])
        log_prior += np.log(rv.pdf(params_dict[param]))
    lp = test_pid_interface.check_prior(params_dict)
    np.testing.assert_allclose(lp, log_prior, rtol = 0.1)

def test_gamma_priors(model_setup):
    """ Gamma prior testing
    """
    M, params_to_estimate = model_setup
    prior = {'m' : ['gamma', 2, 4],'b' : ['gamma', 5, 4]}
    test_pid_interface = PIDInterface(params_to_estimate, M, prior)
    params_dict = {'m':1,'b':2}
    lp = 0
    log_prior = 0
    for param in prior:
        rv = scipy.stats.gamma(prior[param][1], scale = 1/prior[param][2])
        log_prior += np.log(rv.pdf(params_dict[param]))
    lp = test_pid_interface.check_prior(params_dict)
    np.testing.assert_allclose(lp, log_prior, rtol = 0.1)

def test_beta_priors(model_setup):
    """ Beta prior testing
    """
    M, params_to_estimate = model_setup
    prior = {'m' : ['beta', 2, 4],'b' : ['beta', 5, 4]}
    test_pid_interface = PIDInterface(params_to_estimate, M, prior)
    params_dict = {'m':0.1,'b':0.2}
    lp = 0
    log_prior = 0
    for param in prior:
        rv = scipy.stats.beta(prior[param][1], prior[param][2])
        log_prior += np.log(rv.pdf(params_dict[param]))
    lp = test_pid_interface.check_prior(params_dict)
    np.testing.assert_allclose(lp, log_prior, rtol = 0.1)
    
def test_log_uniform_priors(model_setup):
    """ Log-Uniform prior testing
    """
    M, params_to_estimate = model_setup
    prior = {'m' : ['log-uniform', 1, 400],'b' : ['log-uniform', 10, 100]}
    test_pid_interface = PIDInterface(params_to_estimate, M, prior)
    params_dict = {'m':20,'b':20}
    lp = 0
    log_prior = 0
    for param in prior:
        rv = scipy.stats.loguniform(prior[param][1], prior[param][2])
        log_prior += np.log(rv.pdf(params_dict[param]))
    lp = test_pid_interface.check_prior(params_dict)
    print(log_prior)
    print(lp)
    np.testing.assert_allclose(lp, log_prior, rtol = 0.1)

def test_log_gaussian_priors(model_setup):
    """ Log-Gaussian prior testing
    """
    M, params_to_estimate = model_setup
    prior = {'m' : ['log-gaussian', 2, 4],'b' : ['log-gaussian', 5, 4]}
    test_pid_interface = PIDInterface(params_to_estimate, M, prior)
    params_dict = {'m':1,'b':2}
    lp = 0
    log_prior = 0
    for param in prior:
        rv = scipy.stats.lognorm(scale = np.exp(prior[param][1]), s = prior[param][2])
        log_prior += np.log(rv.pdf(params_dict[param]))
    lp = test_pid_interface.check_prior(params_dict)
    np.testing.assert_allclose(lp, log_prior, rtol = 0.1)

@pytest.fixture
def birth_death_process():
    """Set up the Model and data to be used for inference tests
    with state initial conditions varying for each trajectory.
    """
    # Create a bioscrape model
    species = ['I','X', 'Y']
    reactions = [(['X'], [], 'massaction', {'k':'d1'}), 
                ([], ['X'], 'hillpositive', {'s1':'I', 'k':'k1', 'K':'KR', 'n':2}),
                (['X'],['Y'],'massaction', {'k':'k2'})]
    k1 = 50.0
    d1 = 0.5
    k2 = 10
    params = [('k1', k1), ('d1', d1), ('KR', 20), ('k2',k2)]
    initial_condition = {'X':0, 'I':0}
    M = Model(species = species, reactions = reactions, parameters = params, 
            initial_condition_dict = initial_condition)
    num_trajectories = 4 # each with different initial condition
    initial_condition_list = [{'I':5},{'I':10},{'I':15},{'I':20}] 
    timepoints = np.linspace(0,5,100)
    result_list = []
    for init_cond in initial_condition_list:
        M.set_species(init_cond)
        result = py_simulate_model(timepoints, Model = M)
        result_list.append(result)
    exp_data = pd.DataFrame()
    exp_data['timepoints'] = timepoints
    for i in range(num_trajectories):
        exp_data['X' + str(i)] = result_list[i]['X'] + np.random.normal(5, 1, size = np.shape(result_list[i]['X']))
        exp_data['Y' + str(i)] = result_list[i]['Y'] + np.random.normal(5, 1, size = np.shape(result_list[i]['Y']))
    exp_data.to_csv('tests/temp/birth_death_data_multiple_conditions.csv')
    return M, num_trajectories, initial_condition_list

@pytest.fixture
def birth_death_process_params():
    """Set up the Model and data to be used for inference tests
    with parameter conditions varying for each trajectory.
    """
    # Create a bioscrape model
    species = ['I','X']
    reactions = [(['X'], [], 'massaction', {'k':'d1'}), 
                ([], ['X'], 'hillpositive', {'s1':'I', 'k':'k1', 'K':'KR', 'n':2})]
    k1 = 50.0
    d1 = 0.5
    params = [('k1', k1), ('d1', d1), ('KR', 20)]
    initial_condition = {'X':0, 'I':0}
    M = Model(species = species, reactions = reactions, parameters = params, 
            initial_condition_dict = initial_condition)
    num_trajectories = 4 # each with different initial condition
    parameter_condition_list = [{'d1':0.01},{'d1':0.1},{'d1':1},{'d1':10}] 
    timepoints = np.linspace(0,5,100)
    result_list = []
    M.set_species({"I":20})
    for param_cond in parameter_condition_list:
        M.set_params(param_cond)
        result = py_simulate_model(timepoints, Model = M)['X']
        result_list.append(result)
    exp_data = pd.DataFrame()
    exp_data['timepoints'] = timepoints
    for i in range(num_trajectories):
        exp_data['X' + str(i)] = result_list[i] + np.random.normal(5, 1, size = np.shape(result))
    exp_data.to_csv('tests/temp/birth_death_data_multiple_param_conditions.csv')
    return M, num_trajectories, parameter_condition_list

def test_multiple_parameter_conditions_inference(birth_death_process_params):
    M, num_trajectories, parameter_condition_list = birth_death_process_params
    # Import data from CSV
    exp_data = []
    for i in range(num_trajectories):
        df = pd.read_csv('tests/temp/birth_death_data_multiple_param_conditions.csv', usecols = ['timepoints', 'X'+str(i)])
        df.columns = ['timepoints', 'X']
        exp_data.append(df)
        
    prior = {'k1' : ['uniform', 0, 100]}
    sampler, pid = py_inference(Model = M, exp_data = exp_data, measurements = ['X'], time_column = ['timepoints'],
                                parameter_conditions = parameter_condition_list,
                                nwalkers = 5, init_seed = 0.15, nsteps = 2000, sim_type = 'deterministic',
                                convergence_check = True,
                                params_to_estimate = ['k1'], prior = prior)
    assert(isinstance(sampler, EnsembleSampler) == True)
    assert(isinstance(pid, InferenceSetup) == True)
    assert np.array(sampler.get_autocorr_time())[0] < 50 
    assert np.array(sampler.acceptance_fraction).all() < 2

def test_multiple_conditions_inference(birth_death_process):
    """ Run MCMC inference on birth death model with data for multiple conditions

    Args:
        birth_death_process ([pointer]): Returned from the setup function
    """
    M, num_trajectories, initial_condition_list = birth_death_process
    # Import data from CSV
    exp_data = []
    for i in range(num_trajectories):
        df = pd.read_csv('tests/temp/birth_death_data_multiple_conditions.csv', usecols = ['timepoints', 'X'+str(i), 'Y'+str(i)])
        df.columns = ['timepoints', 'X','Y']
        exp_data.append(df)
    prior = {'d1' : ['uniform', 0.1, 10], 'k1' : ['uniform',0,100], 'KR' : ['uniform',0,100], 'k2':['uniform', 0, 100]}
    sampler, pid = py_inference(Model = M, exp_data = exp_data, measurements = ['X','Y'], time_column = ['timepoints'],
                                initial_conditions = initial_condition_list,
                                nwalkers = 8, init_seed = 0.15, nsteps = 100, sim_type = 'deterministic', discard = 10,
                                params_to_estimate = ['d1','k1','KR', 'k2'], prior = prior, convergence_check = False)
    assert(isinstance(sampler, EnsembleSampler) == True)
    assert(isinstance(pid, InferenceSetup) == True)
    # Takes a long time to run the test if convergence is also checked
    #### Recommend to uncomment when inference ####
    #### changes are made to inference module directly. ####
    #### Change nsteps = 9k, when uncommented. ####
    # assert np.array(sampler.get_autocorr_time())[0] < 200 
    # assert np.array(sampler.get_autocorr_time())[1] < 200
    # assert np.array(sampler.get_autocorr_time())[2] < 200
    # assert np.array(sampler.get_autocorr_time())[3] < 200
    # assert np.array(sampler.acceptance_fraction).all() < 2

    
def test_stochastic_inference(birth_death_process):
    """ Run stochastic inference using emcee for the birth death process

    Args:
        birth_death_process ([pointer]): Returned from the setup function
    """
    M, num_trajectories, initial_condition_list = birth_death_process
    # Import data from CSV
    exp_data = []
    for i in range(num_trajectories):
        df = pd.read_csv('tests/temp/birth_death_data_multiple_conditions.csv', usecols = ['timepoints', 'X'+str(i), 'Y'+str(i)])
        df.columns = ['timepoints', 'X', 'Y']
        exp_data.append(df)
    prior = {'d1' : ['uniform', 0.1, 10], 'k1' : ['uniform',0,100], 'KR' : ['uniform',0,100], 'k2':['uniform', 0, 100]}
    sampler, pid = py_inference(Model = M, exp_data = exp_data, measurements = ['X','Y'], time_column = ['timepoints'],
                                initial_conditions = initial_condition_list,
                                nwalkers = 10, init_seed = np.array([0.1,50,20,0.5]), nsteps = 100, sim_type = 'stochastic', discard = 10,
                                params_to_estimate = ['d1','k1','KR', 'k2'], prior = prior, debug=False)
    assert(isinstance(sampler, EnsembleSampler) == True)
    assert(isinstance(pid, InferenceSetup) == True)

def test_lmfit_inference(birth_death_process):
    """Run LMFit package for deterministic inference. 

    Args:
        birth_death_process ([pointer]): As returned by the setup function.
    """
    M, num_trajectories, initial_condition_list = birth_death_process
    # Import data from CSV
    exp_data = []
    for i in range(num_trajectories):
        df = pd.read_csv('tests/temp/birth_death_data_multiple_conditions.csv', usecols = ['timepoints', 'X'+str(i)])
        df.columns = ['timepoints', 'X']
        exp_data.append(df)
    prior = {'k1' : ['uniform', 0, 100]}
    minimizer_result = py_inference(Model = M, exp_data = exp_data, measurements = ['X'], time_column = ['timepoints'],
                                initial_conditions = initial_condition_list, nwalkers = 5, 
                                init_seed = 0.15, nsteps = 40, sim_type = 'deterministic',
                                params_to_estimate = ['k1'], prior = prior, 
                                inference_type = 'lmfit', method = 'leastsq', plot_show = False)
    assert(isinstance(minimizer_result[0], MinimizerResult) == True)
    