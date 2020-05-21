from bioscrape.inference import DeterministicLikelihood as DLL
from bioscrape.inference import StochasticTrajectoriesLikelihood as STLL
from bioscrape.inference import StochasticTrajectories
from bioscrape.inference import BulkData

import numpy as np

def check_priors(param_dict, prior):
    import scipy.stats
    for key,value in param_dict.items():
        prior_type = prior[key][0]
        if prior_type == 'uniform':
            if len(prior[key]) != 3:
                raise ValueError('For uniform distribution, the dictionary entry must be : [prior_type, lower_bound, upper_bound]')
            lb = prior[key][1]
            ub = prior[key][2]
            if value > ub or value < lb:
                return False
        elif prior_type == 'gaussian':
            if len(prior[key]) != 4:
                raise ValueError('For Gaussian distribution, the dictionary entry must be : [prior_type, mean, std_dev, threshold]')
            mu = prior[key][1]
            sig = prior[key][2]
            prob_threshold = prior[key][3]
            # Check if value lies is a valid sample of (mu, sigma) normal distribution
            # Using probability density function for normal distribution
            prob = 1/(np.sqrt(2*np.pi) * sig) * np.exp(-0.5*(value - mu)**2/sig**2)
            if prob < prob_threshold:
                return False
    return True

class PIDInterface(object):
    def __init__(self, params_to_estimate, M, priors):
        self.params_to_estimate = params_to_estimate
        self.M = M
        self.priors = priors
        return
# Add a new class similar to this to create new interfaces.
class StochasticInference(object):
    def __init__(self, params_to_estimate, M, priors):
        self.params_to_estimate = params_to_estimate
        self.M = M
        self.priors = priors
        return

    def get_likelihood_function(self, params_values, data, timepoints, measurements, initial_conditions, norm_order = 2, N_simulations = 3, debug = False):
        M = self.M
        params_dict = {}
        # params_exp = np.exp(log_params)
        # if debug:
        #     print(params)
        for key, p in zip(self.params_to_estimate, params_values):
            params_dict[key] = p
        # Priors
        priors = self.priors
        # Check prior
        if check_priors(params_dict, priors) is False:
            return -np.inf

        N = np.shape(data)[0]

        if debug:
            print('The timepoints shape is {0}'.format(np.shape(timepoints)))
            print('The data shape is {0}'.format(np.shape(data)))
            print('The measurmenets is {0}'.format(measurements))
            print('The N is {0}'.format(N))

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
            return LL_stoch.py_log_likelihood()

       

# Add a new class similar to this to create new interfaces.
class DeterministicInference(object):
    def __init__(self, params_to_estimate, M, priors):
        self.params_to_estimate = params_to_estimate
        self.M = M
        self.priors = priors
        return

    def get_likelihood_function(self, params_values, data, timepoints, measurements, initial_conditions, norm_order = 2, debug = False):
        M = self.M
        params_dict = {}
        # params_exp = np.exp(log_params)
        for key, p in zip(self.params_to_estimate, params_values):
            params_dict[key] = p
        priors = self.priors
        # Check prior
        if check_priors(params_dict, priors) is False:
            return -np.inf
        else:
            lp = 0.0

        N = np.shape(data)[0]
        #Ceate Likelihood objects:
        # Create a data Objects
        # In this case the timepoints should be a list of timepoints vectors for each iteration
        dataDet = BulkData(np.array(timepoints), data, measurements, N)
        #If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
        #Note len(initial_conditions) must be equal to the number of trajectories N
        if debug:
            print('The timepoints shape is {0}'.format(np.shape(timepoints)))
            print('The data shape is {0}'.format(np.shape(data)))
            print('The measurmenets is {0}'.format(measurements))
            print('The N is {0}'.format(N))
        # TODO: Initial conditions not going through correctly?
        # TODO: priors needs all parameters of the models when it should only need those that are being identified. Attempted a fix.
        # TODO: Need to fix how multiple initial conditions will be handled because in pid_interfaces only one at a time can go through.
        LL_det = DLL(model = M, init_state = initial_conditions,
        data = dataDet, norm_order = norm_order)
        #Multiple samples with a single initial only require a single initial condition.
        # Set params here and return the likelihood object.
        if LL_det:
            # if debug:
            #     print('setting {0} to LL_det object'.format(params_dict))
            LL_det.set_init_params(params_dict)
            LL_det_cost = LL_det.py_log_likelihood()
            return lp + LL_det_cost
        




