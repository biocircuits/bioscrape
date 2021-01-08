from bioscrape.inference import DeterministicLikelihood as DLL
from bioscrape.inference import StochasticTrajectoriesLikelihood as STLL
from bioscrape.inference import StochasticTrajectories
from bioscrape.inference import BulkData
import warnings
import numpy as np

class PIDInterface():
    '''
    PID Interface : Parameter identification interface.
    Super class to create parameter identification (PID) interfaces. Two PID interfaces currently implemented: 
    Deterministic and Stochastic inference using time-series data.
    To add a new PIDInterface - simply add a new subclass of this parent class with your desired 
    log-likelihood functions. You can even have your own check_prior function in that class if you do not 
    prefer to use the built in priors with this package.
    '''
    def __init__(self, params_to_estimate, M, prior):
        '''
        Parent class for all PID interfaces.
        Arguments:
        * `params_to_estimate` : List of parameter names to be estimated 
        * `M` : The bioscrape Model object to use for inference
        * `prior` : A dictionary specifying prior distribution. 
        Two built-in prior functions are `uniform_prior` and `gaussian_prior`.
        Each prior has its own syntax for accepting the distribution parameters in the dictionary. 
        New priors may be added. The suggested format for prior dictionaries:
        prior_dict = {'parameter_name': ['prior_name', prior_distribution_parameters]}
        For built-in uniform prior, use {'parameter_name':['uniform', lower_bound, upper_bound]}
        For built-in gaussian prior, use {'parameter_name':['gaussian', mean, standard_deviation, probability threshold]}

        New PID interfaces can be added by creating child classes of PIDInterface class as shown for 
        Built-in PID interfaces : `StochasticInference` and `DeterministicInference`
        '''
        self.params_to_estimate = params_to_estimate
        self.M = M
        self.prior = prior
        return
    
    def check_prior(self, params_dict):
        '''
        To add new prior functions: simply add a new function similar to ones that exist and then 
        call it here.
        '''
        lp = 0.0
        for key,value in params_dict.items():
            if 'positive' in self.prior[key] and value  < 0:
                return np.inf
            prior_type = self.prior[key][0]
            if prior_type == 'uniform':
                lp += self.uniform_prior(key, value)
            elif prior_type == 'gaussian':
                lp += self.gaussian_prior(key, value)
            elif prior_type == 'exponential':
                lp += self.exponential_prior(key, value)
            elif prior_type == 'gamma':
                lp += self.gamma_prior(key, value)
            elif prior_type == 'log-uniform':
                lp += self.log_uniform_prior(key, value)
            elif prior_type == 'log-gaussian':
                lp += self.log_gaussian_prior(key, value)
            elif prior_type == 'beta':
                lp += self.beta_prior(key, value)
            elif prior_type == 'custom':
                # The last element in the prior dictionary must be a callable function
                # The callable function shoud have the following signature :
                # Arguments: param_name (str), param_value(float) 
                # Returns: log prior probability (float or numpy inf)
                custom_fuction = self.prior[key][-1]
                lp += custom_fuction(key, value)
            else:
                raise ValueError('Prior type undefined.')
        return lp

    def uniform_prior(self, param_name, param_value):
        '''
        Check if given param_value is valid according to the prior distribution.
        Returns np.Inf if the param_value is outside the prior range and 0.0 if it is inside. 
        param_name is used to look for the parameter in the prior dictionary.
        '''
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        lower_bound = prior_dict[param_name][1]
        upper_bound = prior_dict[param_name][2]
        if param_value > upper_bound or param_value < lower_bound:
            return np.inf
        else:
            return np.log( 1/(upper_bound - lower_bound) )

    def gaussian_prior(self, param_name, param_value):
        '''
        Check if given param_value is valid according to the prior distribution.
        Returns the log prior probability or np.Inf if the param_value is invalid. 
        '''
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        mu = prior_dict[param_name][1]
        sigma = prior_dict[param_name][2]
        if sigma < 0:
            raise ValueError('The standard deviation must be positive.')
        # Using probability density function for normal distribution
        # Using scipy.stats.norm has overhead that affects speed up to 2x
        prob = 1/(np.sqrt(2*np.pi) * sigma) * np.exp(-0.5*(param_value - mu)**2/sigma**2)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking Gaussian prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

    def exponential_prior(self, param_name, param_value):
        '''
        Check if given param_value is valid according to the prior distribution.
        Returns the log prior probability or np.inf if the param_value is invalid. 
        '''
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        lambda_p = prior_dict[param_name][1]

        prob = lambda_p * np.exp(-lambda_p * param_value)

        if prob < 0:
            warnings.warn('Probability less than 0 while checking Exponential prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)
    
    def gamma_prior(self, param_name, param_value):
        '''
        Check if given param_value is valid according to the prior distribution.
        Returns the log prior probability or np.inf if the param_value is invalid. 
        '''
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        alpha = prior_dict[param_name][1]
        beta = prior_dict[param_name][2]
        from scipy.special import gamma
        prob = (beta**alpha)/gamma(alpha) * param_value**(alpha - 1) * np.exp(-1 * beta*param_value)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking Exponential prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

    def beta_prior(self, param_name, param_value):
        '''
        Check if given param_value is valid according to the prior distribution.
        Returns the log prior probability or np.inf if the param_value is invalid. 
        '''
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        alpha = prior_dict[param_name][1]
        beta = prior_dict[param_name][2]
        import scipy.special.beta as beta_func
        prob = (param_value**(alpha-1) * (1 - param_value)**(beta - 1) )/beta_func(alpha, beta)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking Exponential prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))
            return np.inf
        else:
            return np.log(prob)

    def log_uniform_prior(self, param_name, param_value):
        '''
        Check if given param_value is valid according to the prior distribution.
        Returns the log prior probability or np.inf if the param_value is invalid. 
        '''
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        lower_bound = prior_dict[param_name][1]
        upper_bound = prior_dict[param_name][2]

        if lower_bound < 0 or upper_bound < 0:
            raise ValueError('Upper and lower bounds for log-uniform prior must be positive.')

        if param_value > upper_bound or param_value < lower_bound:
            return np.inf

        prob = 1/(param_value* (np.log(upper_bound) - np.log(lower_bound)))

        if prob < 0:
            warnings.warn('Probability less than 0 while checking Log-Uniform prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))

            return np.inf
        else:
            return np.log(prob)

    def log_gaussian_prior(self, param_name, param_value):
        '''
        Check if given param_value is valid according to the prior distribution.
        Returns the log prior probability or np.inf if the param_value is invalid. 
        '''
        prior_dict = self.prior
        if prior_dict is None:
            raise ValueError('No prior found')
        mu = prior_dict[param_name][1]
        sigma = prior_dict[param_name][2]
        if sigma < 0:
            raise ValueError('The standard deviation must be positive.')
        # Using probability density function for log-normal distribution
        prob = 1/(param_value * np.sqrt(2*np.pi) * sigma) * np.exp((-0.5 * (np.log(param_value) - mu)**2)/sigma**2)
        if prob < 0:
            warnings.warn('Probability less than 0 while checking log-normal prior! Current parameter name and value: {0}:{1}.'.format(param_name, param_value))

            return np.inf
        else:
            return np.log(prob)

# Add a new class similar to this to create new interfaces.
class StochasticInference(PIDInterface):
    def __init__(self, params_to_estimate, M, prior):
        self.LL_stoch = None
        self.dataStoch = None
        super().__init__(params_to_estimate, M, prior)
        return

    def setup_likelihood_function(self, data, timepoints, measurements, initial_conditions, norm_order = 2, N_simulations = 3, debug = False, **kwargs):
        N = np.shape(data)[0]
        if debug:
            print('Stochastic inference attributes:')
            print('The timepoints shape is {0}'.format(np.shape(timepoints)))
            print('The data shape is {0}'.format(np.shape(data)))
            print('The measurmenets is {0}'.format(measurements))
            print('The N is {0}'.format(N))
            print('Using the initial conditions: {0}'.format(initial_conditions))
        self.dataStoch = StochasticTrajectories(np.array(timepoints), data, measurements, N)
        #If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
        #Note len(initial_conditions) must be equal to the number of trajectories N
        self.LL_stoch = STLL(model = self.M, init_state = initial_conditions,
        data = self.dataStoch, N_simulations = N_simulations, norm_order = norm_order)

    def get_likelihood_function(self, params):
        # Set params here and return the likelihood object.
        if self.LL_stoch is None:
            raise RuntimeError("Must call StochasticInference.setup_likelihood_function before using StochasticInference.get_likelihood_function.")

        #Set params
        params_dict = {}
        for key, p in zip(self.params_to_estimate, params):
            params_dict[key] = p
        self.LL_stoch.set_init_params(params_dict)

        #Prior
        lp = self.check_prior(params_dict)
        if not np.isfinite(lp):
            return -np.inf

        LL_stoch_cost = self.LL_stoch.py_log_likelihood()
        ln_prob = lp + LL_stoch_cost
        return ln_prob
       
# Add a new class similar to this to create new interfaces.
class DeterministicInference(PIDInterface):
    def __init__(self, params_to_estimate, M, prior):
        self.LL_det = None
        self.dataDet = None
        super().__init__(params_to_estimate, M, prior)
        return

    def setup_likelihood_function(self, data, timepoints, measurements, initial_conditions, norm_order = 2, debug = False, **kwargs):
        N = np.shape(data)[0]
        #Create a data Objects
        # In this case the timepoints should be a list of timepoints vectors for each iteration
        self.dataDet = BulkData(np.array(timepoints), data, measurements, N)
        #If there are multiple initial conditions in a data-set, should correspond to multiple initial conditions for inference.
        #Note len(initial_conditions) must be equal to the number of trajectories N
        if debug:
            print('The deterministic inference attributes:')
            print('The timepoints shape is {0}'.format(np.shape(timepoints)))
            print('The data shape is {0}'.format(np.shape(data)))
            print('The measurmenets is {0}'.format(measurements))
            print('The N is {0}'.format(N))
            print('Using the initial conditions: {0}'.format(initial_conditions))


        #Create Likelihood object
        self.LL_det = DLL(model = self.M, init_state = initial_conditions, data = self.dataDet, norm_order = norm_order)

    def get_likelihood_function(self, params):
        if self.LL_det is None:
            raise RuntimeError("Must call DeterministicInference.setup_likelihood_function before using DeterministicInference.get_likelihood_function.")
        #this part is the only part that is called repeatedly
        params_dict = {}
        for key, p in zip(self.params_to_estimate, params):
            params_dict[key] = p
        self.LL_det.set_init_params(params_dict)

        # Check prior
        lp = 0
        lp = self.check_prior(params_dict)
        if not np.isfinite(lp):
            return -np.inf

        #apply cost function
        LL_det_cost = self.LL_det.py_log_likelihood()
        ln_prob = lp + LL_det_cost
        return ln_prob
        




