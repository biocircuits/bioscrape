
from bioscrape.inference import DeterministicLikelihood as DLL
from bioscrape.inference import BulkData
from bioscrape.inference import StochasticTrajectoriesLikelihood as STLL
from bioscrape.inference import StochasticTrajectories
import numpy as np


# Add a new class similar to this to create new interfaces.
class StochasticInference(object):
    def __init__(self, params_to_estimate, M, priors):
        self.params = params_to_estimate
        self.M = M
        self.priors = priors
        self.InitialConditions = None 
        self.MultipleTimepoints = None 
        self.MultipleInitialConditions = None
        self.N_simulations = 3
        self.norm_order = 2
        return

    def get_likelihood_function(self, log_params, data, timepoints, measurements, initial_conditions):
        self.InitialConditions = initial_conditions
        M = self.M
        N = np.shape(data)[0]
        #Ceate Likelihood objects:
        N_simulations = self.N_simulations # Number of simulations per sample to compare to
        norm_order = self.norm_order # (integer) Which norm to use: 1-Norm, 2-norm, etc.
        if type(initial_conditions) is list and initial_conditions:
            self.MultipleInitialConditions = True 
        elif type(initial_conditions) is np.ndarray and np.shape(initial_conditions)[0]:
            self.MultipleInitialConditions = False
        else:
            raise ValueError('MultipleInitialConditions attribute is not set.')
        if type(timepoints) is list and timepoints:
            self.MultipleTimepoints = True 
        elif type(timepoints) is np.ndarray and np.shape(timepoints)[0]:
            self.MultipleTimepoints = False
        else:
            raise ValueError('MultipleTimepoints attribute is not set.')
        # Create a data Objects
        if self.MultipleTimepoints:
            # In this case the timepoints should be a list of timepoints vectors for each iteration
            dataStoch = StochasticTrajectories(np.array(timepoints), data, measurements, N)
        else:
            dataStoch = StochasticTrajectories(timepoints, data, measurements, N)
        #If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
        #Note len(X0_list) must be equal to the number of trajectories N
        X0_list = self.InitialConditions
        initial_condition = self.InitialConditions
        if self.MultipleInitialConditions:
            LL_stoch = STLL(model = M, init_state = X0_list,
            data = dataStoch, N_simulations = N_simulations, norm_order = norm_order)
        #Multiple samples with a single initial only require a single initial condition.
        else:
            LL_stoch = STLL(model = M, init_state = initial_condition,
            data = dataStoch, norm_order = norm_order, N_simulations = N_simulations)
        params_dict = {}
        params_exp = np.exp(log_params)
        for key, p in zip(M.get_params2index().keys(),params_exp):
            params_dict[key] = p
        # Priors (uniform priors only implemented)
        priors = self.priors
        # Check prior
        if self.check_priors(params_dict, priors) == False:
            return -np.inf
        # Set params here and return the likelihood object.
        if LL_stoch:
            LL_stoch.set_init_params(params_dict)
            return -LL_stoch.py_log_likelihood()

    def check_priors(self, param_dict, prior):
        for key,value in param_dict.items():
            range = prior[key]
            if value > max(range) or value < min(range):
                return False
        return True
        

# Add a new class similar to this to create new interfaces.
class DeterministicInference(object):
    def __init__(self, params_to_estimate, M, priors):
        self.params = params_to_estimate
        self.M = M
        self.priors = priors
        self.InitialConditions = None 
        self.MultipleTimepoints = None 
        self.MultipleInitialConditions = None 
        self.norm_order = 2
        return

    def get_likelihood_function(self, log_params, data, timepoints, measurements, initial_conditions):
        # print('even here')
        self.InitialConditions = initial_conditions
        M = self.M
        N = np.shape(data)[0]
        if type(initial_conditions) is list and initial_conditions:
            self.MultipleInitialConditions = True 
        elif type(initial_conditions) is np.ndarray and np.shape(initial_conditions)[0]:
            self.MultipleInitialConditions = False
        else:
            raise ValueError('MultipleInitialConditions attribute is not set.')
        if type(timepoints) is list and timepoints:
            self.MultipleTimepoints = True 
        elif type(timepoints) is np.ndarray and np.shape(timepoints)[0]:
            self.MultipleTimepoints = False
        else:
            raise ValueError('MultipleTimepoints attribute is not set.')
 
        #Ceate Likelihood objects:
        norm_order = self.norm_order # (integer) Which norm to use: 1-Norm, 2-norm, etc.
        # Create a data Objects
        if self.MultipleTimepoints:
            # In this case the timepoints should be a list of timepoints vectors for each iteration
            dataDet = BulkData(np.array(timepoints), data, measurements, N)
        else:
            dataDet = BulkData(timepoints, data, measurements, N)


        #If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
        #Note len(X0_list) must be equal to the number of trajectories N
        X0_list = self.InitialConditions
        initial_condition = self.InitialConditions
        if self.MultipleInitialConditions:
            LL_det = DLL(model = M, init_state = X0_list,
            data = dataDet, norm_order = norm_order)
        #Multiple samples with a single initial only require a single initial condition.
        else:
            LL_det = DLL(model = M, init_state = initial_condition,
            data = dataDet, norm_order = norm_order)
        params_dict = {}
        params_exp = np.exp(log_params)
        for key, p in zip(M.get_params2index().keys(),params_exp):
            params_dict[key] = p
        # Priors (uniform priors only implemented)
        priors = self.priors
        # Check prior
        if self.check_priors(params_dict, priors) == False:
            return -np.inf
        # Set params here and return the likelihood object.
        if LL_det:
            LL_det.set_init_params(params_dict)
            return -LL_det.py_log_likelihood()

    def check_priors(self, param_dict, prior):
        for key,value in param_dict.items():
            range = prior[key]
            if value > max(range) or value < min(range):
                return False
        return True
        