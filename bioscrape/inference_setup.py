import numpy as np
import warnings
try:
    import pandas as pd
except:
    warnings.warn('Pandas package not found.')
import sys
import matplotlib.pyplot as plt
from bioscrape.types import Model
from bioscrape.sbmlutil import import_sbml as sbmlutil_import_sbml
from bioscrape.simulator import ModelCSimInterface, DeterministicSimulator, SSASimulator
from bioscrape.pid_interfaces import *

def initialize_inference(**kwargs):
    return InferenceSetup(**kwargs)

class InferenceSetup(object):
    def __init__(self, **kwargs):
        self.M = None
        self.M = kwargs.get('Model', None)
        self.params_to_estimate = kwargs.get('params_to_estimate', [])
        self.init_seed = kwargs.get('init_seed', 0.01)
        self.prior = kwargs.get('prior', None)
        self.nwalkers = kwargs.get('nwalkers', 100)
        self.nsteps = kwargs.get('nsteps', 1000)
        self.dimension = kwargs.get('dimension', 1)
        self.exp_data = kwargs.get('exp_data', None)
        self.sim_type = kwargs.get('sim_type', 'deterministic')
        self.method = kwargs.get('method', 'emcee')
        self.timepoints = kwargs.get('timepoints', None)
        self.time_column = kwargs.get('time_column', 'time')
        self.measurements = kwargs.get('measurements', [''])
        self.initial_conditions = kwargs.get('initial_conditions', None)
        self.parameter_conditions = kwargs.get('parameter_conditions', None)
        self.norm_order = kwargs.get('norm_order', 2)
        self.N_simulations = kwargs.get('N_simulations', 3)
        self.LL_data = None
        self.debug = kwargs.get('debug', False)
        self.cost_progress = []
        self.cost_params = []
        self.hmax = kwargs.get('hmax', None)
        if self.exp_data is not None:
            self.prepare_inference()
            self.setup_cost_function()
        return 

    #Whenever new settings are updated in the constructor, please add them to getstate and setstate
    def __getstate__(self):
        return (
            self.M,
            self.params_to_estimate,
            self.init_seed,
            self.prior,
            self.nwalkers,
            self.nsteps,
            self.dimension,
            self.exp_data,
            self.sim_type,
            self.method,
            self.timepoints,
            self.time_column,
            self.measurements,
            self.initial_conditions,
            self.parameter_conditions,
            self.norm_order,
            self.N_simulations,
            self.LL_data,
            self.debug,
            self.cost_progress,
            self.cost_params,
            self.hmax
            )

    def __setstate__(self, state):
       
        self.M = state[0]
        self.params_to_estimate = state[1]
        self.init_seed = state[2]
        self.prior = state[3]
        self.nwalkers = state[4]
        self.nsteps = state[5]
        self.dimension = state[6]
        self.exp_data = state[7]
        self.sim_type = state[8]
        self.method = state[9]
        self.timepoints = state[10]
        self.time_column = state[11]
        self.measurements = state[12]
        self.initial_conditions = state[13]
        self.parameter_conditions = state[14]
        self.norm_order = state[15]
        self.N_simulations = state[16]
        self.LL_data = state[17]
        self.debug = state[18]
        self.cost_progress = state[19]
        self.cost_params = state[20]
        self.hmax = state[21]
        if self.exp_data is not None:
            self.prepare_inference()
            self.setup_cost_function()
            

    def set_model(self, M):
        '''
        Set the bioscrape Model to the inference object
        '''
        self.M = M
        return True

    def get_model(self):
        '''
        Get the bioscrape Model for the inference object
        '''
        return self.M 

    def set_prior(self, prior):
        '''
        Set the prior distribution for the parameter inference
        '''
        self.prior = prior
        if len(list(self.prior.keys())) != len(self.params_to_estimate):
            raise ValueError('Prior keys length must be equal to the length of params_to_estimate.')
        return True 

    def set_nsteps(self, nsteps: int):
        '''
        Set the nsteps for the parameter inference using emcee
        '''
        self.nsteps = nsteps 
        return True 
        
    def set_nwalkers(self, nwalkers: int):
        '''
        Set the number of walkers for the parameter inference using emcee
        '''
        self.nwalkers = nwalkers 
        return True 

    def set_dimension(self, dimension: int):
        '''
        Set the dimension of parameter set to identify 
        '''
        self.dimension = dimension 
        return True 

    def set_init_seed(self, init_seed: float):
        '''
        Set the init_seed of parameter set to initialize the parameter identification routine
        '''
        self.init_seed = init_seed 
        return True 

    def set_timepoints(self, timepoints):
        '''
        Set (force) the timepoints to use for parameter inference. By default, the time points are loaded in from exp_data
        '''
        self.timepoints = timepoints 
        return True 

    def set_params_to_estimate(self, params_to_estimate: list):
        '''
        Set the list of parameters to estimate
        '''
        self.params_to_estimate = params_to_estimate 
        return True 

    def set_sim_type(self, sim_type: str):
        '''
        Set the sim_type of simulations to run (deterministic or stochastic) when doing parameter inference
        '''
        self.sim_type = sim_type 
        return True 

    def set_method(self, method: str):
        '''
        Set the parameter identification method to use. 
        Supported method keywords:
        * 'emcee' : Uses Python emcee: https://emcee.readthedocs.io/en/stable/
        * 'lmfit' : Uses Python non-linear least squares package: https://lmfit.github.io/lmfit-py/
        '''
        self.method = method 
        return True 
        
    def set_N_simulations(self, N_simulations: int):
        '''
        Set the number of simulations to run for each condition when doing stochastic inference.
        '''
        if self.sim_type == 'stochastic':
            self.N_simulations = N_simulations 
            return True 
        else:
            warnings.warn('N_simulations needs to be set only for stochastic inference which is not currently the case.')

    def set_initial_conditions(self, initial_conditions):
        '''
        Set the initial conditions for parameter inference. 
        Must be a dictionary object with keys pointing 
        to the species measurements (in order of self.measurements).
        The dictionary value is the initial condition. 
        A list of such dictionary may be passed
        in corresponding to each measurement. 
        '''
        # Get initial_conditions from the model if not given explicitly
        if initial_conditions is None: 
            initial_conditions = self.M.get_species_dictionary()
        if type(initial_conditions) is list and len(initial_conditions):
            for curr_ic_index in range(len(initial_conditions)):
                ic = initial_conditions[curr_ic_index]
                if type(ic) is not dict:
                    raise ValueError('All entries in the initial condition list must be dictionaries.')
        self.initial_conditions = initial_conditions
        return True 

    def set_parameter_conditions(self, parameter_conditions):
        '''
        Set the parameter conditions for parameter inference. 
        Must be a dictionary object with keys pointing to the 
        parameters that change with each measurement
        Value of dictionary is the parameter value. 
        A list of such dictionaries may be passed in 
        corresponding to each measurement.
        '''
        self.parameter_conditions = parameter_conditions
        return True 


    def set_measurements(self, measurements: list):
        '''
        Set the list of measurements (outputs) to look for in exp_data
        '''
        self.measurements = measurements
        return True 

    def set_time_column(self, time_column: str):
        '''
        Set the time_column attribute to specify the name of the time column in exp_data
        '''
        self.time_column = time_column 
        return True 

    def set_exp_data(self, exp_data):
        '''
        Set the experimental data to the inference object. Must be a list of Pandas data frame objects.
        '''
        if isinstance(exp_data, (pd.DataFrame, list)):
            self.exp_data = exp_data
        else:
            raise ValueError('exp_data must be either a Pandas dataframe or a list of dataframes.')
        return True 

    def set_norm_order(self, norm_order: int):
        '''
        Set the norm_order used to compute log likelihood
        '''
        self.norm_order = norm_order 
        return True 

    def get_parameters(self):
        '''
        Returns the list of parameters to estimate that are set for the inference object
        '''
        return self.params_to_estimate

    def run_mcmc(self, **kwargs):
        self.prepare_inference(**kwargs)
        sampler = self.run_emcee(**kwargs)
        return sampler

    def prepare_inference(self, **kwargs):
        timepoints = kwargs.get('timepoints')
        norm_order = kwargs.get('norm_order')
        N_simulations = kwargs.get('N_simulations')
        debug = kwargs.get('debug')
        if N_simulations:
            self.set_N_simulations(N_simulations)
        if norm_order:
            # (integer) Which norm to use: 1-Norm, 2-norm, etc.
            self.set_norm_order(norm_order)
        if debug:
            self.debug = debug
        if timepoints is not None:
            if isinstance(timepoints, (list, np.ndarray)):
                self.set_timepoints(timepoints)
            else:
                raise ValueError('Expected type list or np.ndarray for timepoints.')
        self.prepare_initial_conditions()
        self.prepare_parameter_conditions()
        self.LL_data = self.extract_data()

    def prepare_initial_conditions(self, ):
        # Create initial conditions as required
        N = 1 if type(self.exp_data) is dict else len(self.exp_data)
        if type(self.initial_conditions) is dict:
            all_initial_conditions = [self.initial_conditions]*N
        elif type(self.initial_conditions) is list:
            if len(self.initial_conditions) != N:
                raise ValueError('For a list of initial conditions,'
                                    'each item must be a dictionary and'
                                    'the length of the list must be the'
                                    'same as the number of trajectories.')
            all_initial_conditions = self.initial_conditions
        self.initial_conditions = all_initial_conditions
        return

    def prepare_parameter_conditions(self):
        # Create parameter conditions as required
        N = 1 if type(self.exp_data) is dict else len(self.exp_data)
        if type(self.parameter_conditions) is dict:
            all_parameter_conditions = [self.parameter_conditions]*N
        elif type(self.parameter_conditions) is list:
            if len(self.parameter_conditions) != N:
                raise ValueError('For a list of parameter conditions,'
                                    'each item must be a dictionary and'
                                    'the length of the list must be the'
                                    'same as the number of trajectories.')
            all_parameter_conditions = self.parameter_conditions
        else:
            all_parameter_conditions = None
        self.parameter_conditions = all_parameter_conditions
        # Make sure that parameters to estimate do not intersect with parameters
        if self.parameter_conditions is not None:
            # that are changing through parameter conditions
            for param_condition in self.parameter_conditions:
                for param in self.params_to_estimate:
                    assert param not in param_condition.keys()
        return

    def extract_data(self):
        exp_data = self.exp_data
        # Get timepoints from given experimental data
        if isinstance(self.timepoints, (list, np.ndarray)):
            warnings.warn('Timepoints given by user, not using the data to extract the timepoints automatically.')
        M = len(self.measurements)# Number of measurements
        if type(exp_data) is list:
            if len(exp_data) == 1:
                exp_data = exp_data[0]
        # Multiple trajectories case 
        if type(exp_data) is list:
            N = len(exp_data)# Number of trajectories
            data_list_final = []
            timepoints_list = []
            for df in exp_data:
                data_list = []
                if type(df) is not pd.DataFrame:
                    raise TypeError('All elements of exp_data attribute of an InferenceSetup object must be Pandas DataFrame objects.')
                # Extract timepoints
                if self.time_column:
                    timepoint_i = np.array(df.get(self.time_column), dtype='double').flatten()
                    timepoints_list.append(timepoint_i)
                else:
                    raise TypeError('time_column attribute of InferenceSetup object must be a string.')
                # Extract measurements    
                if type(self.measurements) is list and len(self.measurements) == 1:
                    data_list.append(np.array(df.get(self.measurements[0]), dtype='double'))
                elif type(self.measurements) is list and len(self.measurements) > 1:
                    for m in self.measurements:
                        # Error in multiple measurements
                        data_list.append(np.array(df.get(m), dtype='double'))
                # Number of timepoints
                T = len(timepoints_list[0])
                if T != len(timepoint_i):
                    warnings.warn('The length of timepoints for all experimental trajectories must be the same, they can have different timepoints but not length of timepoints.')
                data_i = np.array(data_list)
                data_i = np.reshape(data_i, (T, M))
                data_list_final.append(data_i)
            data = np.array(data_list_final)
            self.timepoints = timepoints_list
            T = len(timepoints_list[0])
            data = np.reshape(data, (N,T,M))
            if self.debug:
                print('N (Number of trajectories) = {0}'.format(N))
                print('T (Length of timepoints) = {0}'.format(T))
                print('M (Number of measured species) = {0}'.format(M))
                print('The shape of data is {0}'.format(np.shape(data)))
            assert np.shape(data) == (N,T,M)
        elif type(exp_data) is pd.DataFrame:
            # Extract time
            if self.time_column:
                self.timepoints = np.array(exp_data.get(self.time_column), dtype='double').flatten()
            else:
                raise TypeError('time_column attribute of InferenceSetup object must be a string.')
            
            # Extract measurements
            if type(self.measurements) is list and len(self.measurements) == 1:
                data = np.array(exp_data.get(self.measurements[0]), dtype='double')
            elif type(self.measurements) is list and len(self.measurements) > 1:
                data_list = []
                for m in self.measurements:
                    data_list.append(np.array(exp_data.get(m), dtype='double'))
                data = np.array(data_list)
            else:
                raise ValueError('Something wrong with experimental data input to inference.')
            N = 1 # Number of trajectories
            T = len(self.timepoints) # Number of timepoints
            M = len(self.measurements)# Number of measurements
            data = np.reshape(data, (N,T,M))
        else:
            raise TypeError('exp_data attribute of InferenceSetup object must'
                            'be a list of Pandas DataFrames or a single Pandas DataFrame. ')
        return data

    def setup_cost_function(self, **kwargs):
        if self.sim_type == 'stochastic':
            self.pid_interface = StochasticInference(self.params_to_estimate, self.M, self.prior, **kwargs)
            self.pid_interface.setup_likelihood_function(self.LL_data, self.timepoints, self.measurements, 
                                                         initial_conditions=self.initial_conditions, 
                                                         parameter_conditions=self.parameter_conditions,
                                                         norm_order = self.norm_order,
                                                         N_simulations = self.N_simulations, **kwargs)
        elif self.sim_type == 'deterministic':
            self.pid_interface = DeterministicInference(self.params_to_estimate, self.M, self.prior, **kwargs)
            self.pid_interface.setup_likelihood_function(self.LL_data, self.timepoints, self.measurements, 
                                                         initial_conditions=self.initial_conditions,
                                                         parameter_conditions=self.parameter_conditions,
                                                         norm_order=self.norm_order, **kwargs)

    def cost_function(self, params):
        if self.pid_interface is None:
            raise RuntimeError("Must call InferenceSetup.setup_cost_function() before InferenceSetup.cost_function(params) can be used.")

        cost_value = self.pid_interface.get_likelihood_function(params)
        self.cost_progress.append(cost_value)
        self.cost_params.append(params)
        return cost_value

    def seed_parameter_values(self, **kwargs):
        if 'init_seed' in kwargs:
            self.set_init_seed(kwargs['init_seed'])
        ndim = len(self.params_to_estimate)
        params_values = []
        for p in self.params_to_estimate:
            value = self.M.get_parameter_dictionary()[p]
            params_values.append(value)
        #Sample a one percent ball around a given initial value
        if (isinstance(self.init_seed, np.ndarray) \
            or isinstance(self.init_seed, list)) \
            and len(self.init_seed) == ndim:
            p0 = np.array(self.init_seed) + 0.01*np.array(self.init_seed) \
                * np.random.randn(self.nwalkers, ndim)
        #Use this exact start value
        elif isinstance(self.init_seed, np.ndarray) \
            and self.init_seed.shape == (self.nwalkers, ndim):
            p0 =  np.array(self.init_seed)
        #Sample the Prior Distributions to determine initial values
        elif self.init_seed == "prior":
            p0 = np.zeros((self.nwalkers, ndim))
            for i, p in enumerate(self.params_to_estimate):
                prior = self.prior[p]
                if prior[0] == "uniform":
                    p0[:, i] = np.random.rand(self.nwalkers)*(prior[2]-prior[1])+prior[1]
                elif prior[0] == "gaussian":
                    p0[:, i] = prior[2]*np.random.randn(self.nwalkers)+prior[1]
                elif prior[0] == "log-uniform":
                    a = np.log(prior[1])
                    b = np.log(prior[2])

                    u = np.random.randn(self.nwalkers)*(b - a)+a
                    p0[:, i] = np.exp(u)
                else:
                    raise ValueError("Can only sample uniform and gaussian priors"
                                     "when 'init_seed' is set to prior. "
                                     "Try setting intial seed to a number [0, 1]"
                                     "to sample a gaussian ball around the model"
                                     "parameters instead.")
        #sample a gaussian ball around the initial model parameters
        elif isinstance(self.init_seed, float):
            p0 = np.array(params_values) + self.init_seed * np.array(params_values) * np.random.randn(self.nwalkers, ndim)
        else:
            raise ValueError("init_seed must be a float (will sample a gaussian ball"
                             " of this percent around the model initial condition), "
                             "array (of size parameters or walkers x parameters), "
                             "or the string 'prior' (will sample from uniform "
                             "and guassian priors)")
        
        #Ensure parameters are positive, if their priors are declared to be positive
        #When working in log space, a small non-zero value is used
        if hasattr(self.pid_interface, "log_space_parameters") and self.pid_interface.log_space_parameters:
            epsilon = 10**-8
        else:
            epsilon = 0

        for i, p in enumerate(self.params_to_estimate):
            if "positive" in self.prior[p]:
                p0[:, i] = p0[:, i]*(p0[:, i] > 0) + (p0[:, i] <= 0)*epsilon


        #convert to log space, if pid_interface.log_space_parameters
        if hasattr(self.pid_interface, "log_space_parameters") and self.pid_interface.log_space_parameters:
            p0 = np.log(p0)
        return p0

    def run_emcee(self, **kwargs):
        if kwargs.get("reuse_likelihood", False) is False: 
            self.setup_cost_function(**kwargs)

        progress = kwargs.get('progress')
        convergence_check = kwargs.get('convergence_check', False)
        convergence_diagnostics = kwargs.get('convergence_diagnostics', convergence_check)
        skip_initial_state_check = kwargs.get('skip_initial_state_check', False)
        progress = kwargs.get('progess', True)
        # threads = kwargs.get('threads', 1)
        fname_csv = kwargs.get('filename_csv', 'mcmc_results.csv')
        if 'results_filename' in kwargs:
            warnings.warn('The keyword results_filename is deprecated and'
                          'is being replaced by filename_csv and filename_txt where CSV is for'
                          'the MCMC samples and cost function progress respectively.', DeprecationWarning)
            fname_csv = kwargs.get('results_filename', 'mcmc_results.csv')
        fname_txt = kwargs.get('filename_txt', 'mcmc_results.txt')
        printout = kwargs.get('printout', True)

        try:
            import emcee
        except:
            raise ImportError('emcee package not installed.')
        ndim = len(self.params_to_estimate)

        p0 = self.seed_parameter_values(**kwargs)

        assert p0.shape == (self.nwalkers, ndim)
        
        pool = kwargs.get('pool', None)
        if printout: print("creating an ensemble sampler with multiprocessing pool=", pool)

        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.cost_function, pool = pool)
        sampler.run_mcmc(p0, self.nsteps, progress=progress,
                         skip_initial_state_check=skip_initial_state_check)
        if convergence_check:
            self.autocorrelation_time = sampler.get_autocorr_time()
        if convergence_diagnostics:
            if not convergence_check:
                warnings.warn('MCMC diagnostics cannot be printed when convergence check is False.')
                self.convergence_diagnostics = {}
            else:
                self.convergence_diagnostics = {'Autocorrelation time for each parameter':self.autocorrelation_time,
                                    'Acceptance fraction (fraction of steps that were accepted)':sampler.acceptance_fraction}
        # Write results
        import csv
        with open(fname_csv,'w', newline = "") as f:
            writer = csv.writer(f)
            writer.writerows(sampler.get_chain(flat = True))
            f.close()
        with open(fname_txt, 'w', newline = "") as f:
            f.write('\nCost function progress\n')
            f.write(str(self.cost_progress))
            if convergence_diagnostics:
                f.write('\nMCMC convergence diagnostics\n')
                f.write(str(self.convergence_diagnostics))
            f.close()
        if printout: print("Results written to" + fname_csv + " and " + fname_txt)
        if printout: print('Successfully completed MCMC parameter identification procedure.'
                           'Check the MCMC diagnostics to evaluate convergence.')
        return sampler
    
    def plot_mcmc_results(self, sampler, **kwargs):
        """Plots MCMC results collected in emcee.EnsembleSampler object as a corner plot.
        Args:
            sampler (emcee.EnsembleSampler): Go to emcee.EnsembleSampler documentation for more.
    
        Returns:
            param_report: A dictionary consisting of true values and uncertainties associated with them
                          for each parameter.
                truth list: The truth values for each parameter are computed as 
                            16th, 50th, and 84th percentiles from sampled data. 
                            The percentiles may be computed at different values by passing
                            in the `percentiles` keyword to the function with a list of values.
                            Default `percentiles = [16,50,84]`.
                uncertainty list: The uncertainties aaround the true values 
                                  computed from the samples.
            figure_objects: A list of three figure related objects:
                            `[fig, axes, corner_fig]` where `fig` and `axes` correspond 
                            to the MCMC chain plots and `corner_fig` is the corner figure
                            object returned by the call to `corner.corner`. If corner plot is not 
                            shown, then only `[fig, axes]` is returned.
        """
        print('Parameter posterior distribution convergence plots:')
        figsize = kwargs.get('figsize', (10,7))
        sharex = kwargs.get('sharex', True)
        ndim = sampler.ndim
        fig, axes = plt.subplots(ndim, figsize=figsize, sharex=sharex)
        figure_objects = []
        samples = sampler.get_chain()
        labels = kwargs.get('labels', list(self.params_to_estimate))
        alpha = kwargs.get('alpha', 0.3)
        for i in range(ndim):
            if type(axes) is np.ndarray:
                ax = axes[i]
            else:
                ax = axes
            ax.plot(samples[:, :, i], "k", alpha=alpha)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
        if type(axes) is np.ndarray:
            axes[-1].set_xlabel("step number")
        else:
            axes.set_xlabel("step number")
        figure_objects.append(fig)
        figure_objects.append(axes)
        # arbitrarily discard 2*nwalkers steps 
        discard = kwargs.get('discard', 2*self.nwalkers)
        thin = int(kwargs.get('thin', 1))
        flat = kwargs.get('flat', True)
        
        flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=flat)
        param_report = {}
        # Percentiles to compute for numpy.percentile
        percentiles = kwargs.get('percentiles', [16, 50, 84]) # Set percentiles to compute by default to q
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], q = percentiles)
            q = np.diff(mcmc)
            param_report[self.params_to_estimate[i] + '_true'] = mcmc[1]
            param_report[self.params_to_estimate[i] + '_uncertainties'] = q
            # uncertainty_list.append([-1.0*q[0], q[1]])
            # uncertainty_list.append(q)
        try:
            import corner
            corner_fig = corner.corner(
                flat_samples, labels=labels,
            )
            figure_objects.append(corner_fig)
        except:
            warnings.warn('corner package not found - cannot plot parameter distributions.')
        return param_report, figure_objects 


    def run_lmfit(self, method = 'leastsq', plot_show = True, **kwargs):
        """
        Run the Python LMFit package for a bioscrape model.
        """
        self.prepare_inference(**kwargs)
        N = np.shape(self.LL_data)[0]
        if self.sim_type == 'stochastic':
            stochastic = True
        else:
            stochastic = False
        self.pid_interface = LMFitInference(self.params_to_estimate, self.M, self.prior)
        minimizer_result = [None]*N
        if N == 1:
            minimizer_result[0] = self.pid_interface.\
                get_minimizer_results(self.LL_data[0,:,:], 
                                      self.timepoints, self.measurements,
                                      initial_conditions = self.initial_conditions,
                                      parameter_conditions = self.parameter_conditions,
                                      stochastic = stochastic, debug = self.debug,
                                      method = method, plot_show = plot_show, **kwargs)
        else:   
            if self.parameter_conditions is None:
                self.parameter_conditions = [None]*N
            for i in range(N):
                minimizer_result[i] = self.pid_interface.\
                    get_minimizer_results(self.LL_data[i,:,:], 
                                          self.timepoints[i], self.measurements,
                                          initial_conditions = self.initial_conditions[i],
                                          parameter_conditions = self.parameter_conditions[i],
                                          stochastic = stochastic, debug = self.debug,
                                          method = method, plot_show = plot_show, **kwargs)
        print('Successfully completed parameter identification'
              ' procedure using LMFit. Parameter values and fitness' 
              ' reports written to lmfit_results.csv file. '
              'Check the minimizer_results object returned to '
              'further statistically evaluate the goodness of fit.')
        return minimizer_result

    def write_lmfit_results(self, minimizer_result, **kwargs):
        """
        Process minimizer_result list obtained from LMFit package
        and plot parameter values on a scatter plot using the corner package.
        Arguments:
        * minimizer_result: List of MinimizerResult object (from LMFit package)
        * kwargs: Passed to the corner package. 
        Output:
        * A CSV file: "lmfit_results.csv" is written with each row consisting of the optimized parameter
        value corresponding to each trajectory in the data.
        * Corner histogram plot showing the parameter values obtained.
        """
        import csv
        try:
            from lmfit import fit_report
        except:
            raise ImportError('Package lmfit not found.')
        values_dict = {}
        for param_name in self.params_to_estimate:
            values_dict[param_name] = []
        for result in minimizer_result:
            for param_name in self.params_to_estimate:
                values_dict[param_name].append(dict(result.params.valuesdict())[param_name])
        df = pd.DataFrame.from_dict(values_dict)
        # df.to_csv('lmfit_results.csv', columns = self.params_to_estimate, index = False)
        with open('lmfit_results.csv','w') as f:
            f.write(str(df))

        convergence_diagnostics = kwargs.get('convergence_diagnostics', True) 
        if convergence_diagnostics:
            count = 0
            for result in minimizer_result:
                with open('lmfit_results.csv','a') as f:
                    # writer = csv.writer(f)
                    f.write('\nFor trajectory: {0}\n'.format(count))
                    f.write(fit_report(result))
                    count += 1
                    f.close()
