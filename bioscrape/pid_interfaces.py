
from bioscrape.inference import DeterministicLikelihood as DLL
from bioscrape.inference import StochasticTrajectoriesLikelihood as STLL
from bioscrape.inference import StochasticTrajectories
from bioscrape.inference import BulkData

import numpy as np

def check_priors(param_dict, prior):
    import scipy.stats
    for key,value in param_dict.items():
        type = prior[key][0]
        if type == 'uniform':
            if len(prior[key]) != 3:
                raise ValueError('For uniform distribution, the dictionary entry must be : [type, lower_bound, upper_bound]')
            lb = prior[key][1]
            ub = prior[key][2]
            if value > ub or value < lb:
                return False
        elif type == 'gaussian':
            if len(prior[key]) != 4:
                raise ValueError('For Gaussian distribution, the dictionary entry must be : [type, mean, std_dev, threshold]')
            mu = prior[key][1]
            sig = prior[key][2]
            prob_threshold = prior[key][3]
            distrib = scipy.stats.norm(mu, sig)
            # Check if value lies is a valid sample of (mu, sigma) Gaussian distribution
            if distrib.pdf(value) < prob_threshold:
                return False

    return True



# Add a new class similar to this to create new interfaces.
class StochasticInference(object):
    def __init__(self, params_to_estimate, M, priors):
        self.params = params_to_estimate
        self.M = M
        self.priors = priors
        return

    def get_likelihood_function(self, log_params, data, timepoints, measurements, initial_conditions, norm_order = 2, N_simulations = 3, debug = False):
        M = self.M
        params_dict = {}
        params_exp = np.exp(log_params)
        for key, p in zip(M.get_params2index().keys(),params_exp):
            params_dict[key] = p
        # Priors (uniform priors only implemented)
        priors = self.priors
        # Check prior
        if check_priors(params_dict, priors) == False:
            return -np.inf

        N = np.shape(data)[0]

        dataStoch = StochasticTrajectories(np.array(timepoints), data, measurements, N)
        #If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
        #Note len(initial_conditions) must be equal to the number of trajectories N
        LL_stoch = STLL(model = M, init_state = initial_conditions,
        data = dataStoch, N_simulations = N_simulations, norm_order = norm_order)
        # Set params here and return the likelihood object.
        if LL_stoch:
            # if debug:
            #     print('setting {0} to LL_stoch object'.format(params_dict))
            LL_stoch.set_init_params(params_dict)
            return -LL_stoch.py_log_likelihood()

       

# Add a new class similar to this to create new interfaces.
class DeterministicInference(object):
    def __init__(self, params_to_estimate, M, priors):
        self.params = params_to_estimate
        self.M = M
        self.priors = priors
        return

    def get_likelihood_function(self, log_params, data, timepoints, measurements, initial_conditions, norm_order = 2, debug = False ):
        M = self.M
        params_dict = {}
        params_exp = np.exp(log_params)
        for key, p in zip(M.get_params2index().keys(),params_exp):
            params_dict[key] = p
        priors = self.priors
        # Check prior
        if check_priors(params_dict, priors) == False:
            return -np.inf

        N = np.shape(data)[0]
        #Ceate Likelihood objects:
        # Create a data Objects
        # In this case the timepoints should be a list of timepoints vectors for each iteration
        dataDet = BulkData(np.array(timepoints), data, measurements, N)
        #If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
        #Note len(initial_conditions) must be equal to the number of trajectories N
        LL_det = DLL(model = M, init_state = initial_conditions,
        data = dataDet, norm_order = norm_order)
        #Multiple samples with a single initial only require a single initial condition.
        # Set params here and return the likelihood object.
        if LL_det:
            # if debug:
            #     print('setting {0} to LL_det object'.format(params_dict))
            LL_det.set_init_params(params_dict)
            return -LL_det.py_log_likelihood()
        




