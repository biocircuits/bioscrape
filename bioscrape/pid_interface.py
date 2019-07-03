
from bioscrape.inference import DeterministicLikelihood as DLL
from bioscrape.inference import BulkData
from bioscrape.inference import StochasticTrajectoriesLikelihood as STLL
from bioscrape.inference import StochasticTrajectories
import numpy as np
class PIDInterface(object):
    def __init__(self, params_to_estimate, M, priors):
        self.params = params_to_estimate
        self.M = M
        self.priors = priors
        self.MultipleTimepoints = False
        self.InitialConditions = M.get_initial_conditions() # TODO (implement this)
        self.MultipleInitialConditions = False
        self.type = 'stochastic'
        self.N_simulations = 3
        self.norm_order = 2
        return 
    def get_likelihood_function(self, log_params, Data, timepoints, measurements):
        M = self.M
        N = np.shape(Data)[0]
        #Ceate Likelihood objects:
        N_simulations = self.N_simulations # Number of simulations per sample to compare to
        norm_order = self.norm_order # (integer) Which norm to use: 1-Norm, 2-norm, etc.
        # Create a Data Objects
        if self.MultipleTimepoints:
            # In this case the timepoints should be a list of timepoints vectors for each iteration
            if self.type == 'stochastic':
                DataStoch = StochasticTrajectories(np.array(timepoints), Data,
                measurements, N)
            elif self.type == 'deterministic':
                DataDet = BulkData(np.array(timepoints), Data, measurements, N)
            else:
                raise ValueError('Invalid type argument in get_likelihood_function call')
        else:
            if self.type == 'stochastic':
                DataStoch = StochasticTrajectories(timepoints, Data, measurements, N)
            elif self.type == 'deterministic':
                DataDet = BulkData(timepoints, Data, measurements, N)
            else:
                raise ValueError('Invalid type argument in get_likelihood_function call')


        #If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
        #Note len(X0_list) must be equal to the number of trajectories N
        X0_list = self.InitialConditions
        initial_condition = self.InitialConditions
        if self.MultipleInitialConditions:
            if self.type == 'stochastic':    
                LL_stoch = STLL(model = M, init_state = X0_list,
                data = DataStoch, N_simulations = N_simulations, norm_order = norm_order)
            elif self.type == 'deterministic':
                LL_det = DLL(model = M, init_state = X0_list,
                data = DataDet, norm_order = norm_order)
            else:
                raise ValueError('Invalid type argument in get_likelihood_function call')

        #Multiple samples with a single initial only require a single initial condition.
        else:
            if self.type == 'stochastic':
                LL_stoch = STLL(model = M, init_state = initial_condition,
                data = DataStoch, norm_order = norm_order, N_simulations = N_simulations)
            elif self.type == 'deterministic':
                LL_det = DLL(model = M, init_state = initial_condition,
                data = DataDet, norm_order = norm_order)
            else:
                raise ValueError('Invalid type argument in get_likelihood_function call')

        # Priors? (TODO)
        priors = self.priors
        # Check prior
        if self.log_prior(param_dict, priors) == False:
            return -np.inf

        # Set params here 
        params_dict = {}
        params_exp = np.exp(log_params)
        for key, p in zip(M.params_dict.keys(),params_exp):
            params_dict[key] = p
        if LL_stoch:
            LL_stoch.set_init_params(params_dict)
            return -LL_stoch.py_log_likelihood()
        if LL_det:
            LL_det.set_init_params(params_dict)
            return -LL_det.py_log_likelihood()

    def log_prior(self, param_dict, prior):
        for key,value in param_dict.items():
            range = prior[key]
            if value > max(range) or value < min(range):
                return False
        return True
        