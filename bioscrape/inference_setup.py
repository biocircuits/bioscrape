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
    obj = InferenceSetup(**kwargs)
    return obj

class InferenceSetup(object):
    def __init__(self, **kwargs):
        self.M = None
        self.pid_interface = None
        if 'Model' in kwargs:
            self.set_model(kwargs.get('Model'))
        self.params_to_estimate = []
        if 'params_to_estimate' in kwargs:
            self.set_params_to_estimate(kwargs.get('params_to_estimate'))
        self.init_seed = 0.01 
        if 'init_seed' in kwargs:
            self.set_init_seed(kwargs.get('init_seed'))
        self.prior = None
        if 'prior' in kwargs:
            self.set_prior(kwargs.get('prior'))
        self.nwalkers = 100
        if 'nwalkers' in kwargs:
            self.set_nwalkers(kwargs.get('nwalkers'))
        self.nsteps = 2000
        if 'nsteps' in kwargs:
            self.set_nsteps(kwargs.get('nsteps'))
        self.dimension = 0
        if 'dimension' in kwargs:
            self.set_dimension(kwargs.get('dimension'))
        self.exp_data = None # Pandas DataFrame object
        if 'exp_data' in kwargs:
            self.set_exp_data(kwargs.get('exp_data'))
        self.sim_type = 'deterministic'
        if 'sim_type' in kwargs:
            self.set_sim_type(kwargs.get('sim_type'))
        self.method = 'emcee'
        if 'method' in kwargs:
            self.set_method(kwargs.get('method'))
        self.timepoints = None
        if 'timepoints' in kwargs:
            self.set_timepoints(kwargs.get('timepoints'))
        self.time_column = ''
        if 'time_column' in kwargs:
            self.set_time_column(kwargs.get('time_column'))
        self.measurements = ['']
        if 'measurements' in kwargs:
            self.set_measurements(kwargs.get('measurements'))
        self.initial_conditions = None
        if 'initial_conditions' in kwargs:
            self.set_initial_conditions(kwargs.get('initial_conditions'))
        self.norm_order = 2
        if 'norm_order' in kwargs:
            self.set_norm_order(kwargs.get('norm_order'))
        self.N_simulations = 3
        if 'N_simulations' in kwargs:
            self.set_N_simulations(kwargs.get('N_simulations'))
        self.LL_data = None
        self.debug = False
        if 'debug' in kwargs:
            self.debug = kwargs.get('init_seed')
        self.cost_progress = []
        self.hmax = kwargs.get('hmax', None)
        return 

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
        Set the initial_conditions for parameter inference. 
        Must be a dictionary object with keys pointing to the species measurements (in order of self.measurements).
        Values of dictionary may be list of different initial conditions for which the data was collected. 
        '''
        if type(initial_conditions) is dict:
            if len(list(initial_conditions.keys())) != len(self.M.get_species_array()):
                new_ic = dict(self.M.get_species_dictionary())
                for key, value in initial_conditions.items():
                    new_ic[key] = value
                self.initial_conditions = new_ic
            else:
                self.initial_conditions = initial_conditions
        elif type(initial_conditions) is list and len(initial_conditions):
            for curr_ic_index in range(len(initial_conditions)):
                ic = initial_conditions[curr_ic_index]
                if type(ic) is not dict:
                    raise ValueError('All entries in the initial condition list must be dictionaries.')
                elif len(list(ic.keys())) != len(self.M.get_species_array()):
                    new_ic = dict(self.M.get_species_dictionary())
                    for key, value in ic.items():
                        new_ic[key] = value
                    initial_conditions[curr_ic_index] = new_ic
            self.initial_conditions = initial_conditions
        # Get initial_conditions from the model if not given explicitly
        initial_conditions = self.initial_conditions
        if initial_conditions is None: 
            self.initial_conditions = self.M.get_species_dictionary()
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
        self.LL_data = self.extract_data(self.exp_data)

    def extract_data(self, exp_data):
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
                    timepoint_i = np.array(df.get(self.time_column)).flatten()
                    timepoints_list.append(timepoint_i)
                else:
                    raise TypeError('time_column attribute of InferenceSetup object must be a string.')

                # Extract measurements    
                if type(self.measurements) is list and len(self.measurements) == 1:
                    data_list.append(np.array(df.get(self.measurements[0])))
                elif type(self.measurements) is list and len(self.measurements) > 1:
                    for m in self.measurements:
                        # Error in multiple measurements
                        data_list.append(np.array(df.get(m)))
                # Number of timepoints
                T = len(timepoints_list[0])
                if T != len(timepoint_i):
                    warnings.warn('The length of timepoints for all experimental trajectories must be the same, they can have different timepoints but not length of timepoints.')
                data_i = np.array(data_list)
                data_i = np.reshape(data_i, (T, M))
                data_list_final.append(data_i)

            # Create initial conditions as required
            if type(self.initial_conditions) is dict:
                all_initial_conditions = [self.initial_conditions]*N
            elif type(self.initial_conditions) is list:
                if len(self.initial_conditions) != N:
                    raise ValueError('For a list of initial conditions, each item must be a dictionary and the length of the list must be the same as the number of trajectories.')
                all_initial_conditions = self.initial_conditions
            self.initial_conditions = all_initial_conditions
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
                self.timepoints = np.array(exp_data.get(self.time_column)).flatten()
            else:
                raise TypeError('time_column attribute of InferenceSetup object must be a string.')
            
            # Extract measurements
            if type(self.measurements) is list and len(self.measurements) == 1:
                data = np.array(exp_data.get(self.measurements[0]))
            elif type(self.measurements) is list and len(self.measurements) > 1:
                data_list = []
                for m in self.measurements:
                    data_list.append(np.array(exp_data.get(m)))
                data = np.array(data_list)
            else:
                raise ValueError('Something wrong with experimental data input to inference.')
            N = 1 # Number of trajectories
            T = len(self.timepoints) # Number of timepoints
            M = len(self.measurements)# Number of measurements
            data = np.reshape(data, (N,T,M))
        else:
            raise TypeError('exp_data attribute of InferenceSetup object must be a list of Pandas DataFrames or a single Pandas DataFrame. ')
        return data

    def setup_cost_function(self, **kwargs):
        if self.sim_type == 'stochastic':
            self.pid_interface = StochasticInference(self.params_to_estimate, self.M, self.prior, **kwargs)
            self.pid_interface.setup_likelihood_function(self.LL_data, self.timepoints, self.measurements, 
                                                            self.initial_conditions, norm_order = self.norm_order, 
                                                            N_simulations = self.N_simulations, debug = self.debug, **kwargs)
        elif self.sim_type == 'deterministic':
            self.pid_interface = DeterministicInference(self.params_to_estimate, self.M, self.prior, **kwargs)
            self.pid_interface.setup_likelihood_function(self.LL_data, self.timepoints, self.measurements, 
                                                            self.initial_conditions, norm_order = self.norm_order, debug = self.debug, **kwargs)

    def cost_function(self, params):
        if self.pid_interface is None:
            raise RuntimeError("Must call InferenceSetup.setup_cost_function() before InferenceSetup.cost_function(params) can be used.")

        cost_value = self.pid_interface.get_likelihood_function(params)
        self.cost_progress.append(cost_value)
        return cost_value

    def run_emcee(self, **kwargs):
        self.setup_cost_function(**kwargs)
        progress = kwargs.get('progress')
        convergence_check = kwargs.get('convergence_check', True)
        convergence_diagnostics = kwargs.get('convergence_diagnostics', convergence_check)
        skip_initial_state_check = kwargs.get('skip_initial_state_check', False)
        progress = kwargs.get('progess', True)
        threads = kwargs.get('threads', 1)
        #if not 'convergence_check' in kwargs:
        #    convergence_check = True
        #if not 'convergence_diagnostics' in kwargs:
        #    convergence_diagnostics = True
        #    if not convergence_check:
        #        convergence_diagnostics = False
        #if not 'progress' in kwargs:
        #    progress = True

        try:
            import emcee
        except:
            raise ImportError('emcee package not installed.')
        ndim = len(self.params_to_estimate)
        params_values = []
        for p in self.params_to_estimate:
            value = self.M.get_parameter_dictionary()[p]
            params_values.append(value)
        if isinstance(self.init_seed, np.ndarray) or isinstance(self.init_seed, list):
            p0 = np.array(self.init_seed) + 0.01*np.random.randn(self.nwalkers, ndim)
        else:
            p0 = np.array(params_values) + self.init_seed * np.random.randn(self.nwalkers, ndim)
        assert p0.shape == (self.nwalkers, ndim)
        print("creating an ensemble sampler with threads=", threads)
        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.cost_function, threads = threads)
        sampler.run_mcmc(p0, self.nsteps, progress = progress, skip_initial_state_check = skip_initial_state_check)
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
        with open('mcmc_results.csv','w', newline = "") as f:
            writer = csv.writer(f)
            writer.writerows(sampler.get_chain(flat = True))
            if convergence_diagnostics:
                writer.writerow('\nMCMC convergence diagnostics\n')
                writer.writerow(self.convergence_diagnostics)
            writer.writerow('\nCost function progress\n')
            writer.writerow(self.cost_progress)
            f.close()
        print('Successfully completed MCMC parameter identification procedure. Parameter distribution data written to mcmc_results.csv file. Check the MCMC diagnostics to evaluate convergence.')
        return sampler
    
    def plot_mcmc_results(self, sampler, plot_show = True, **kwargs):
        print('Parameter posterior distribution convergence plots:')
        ndim = sampler.ndim
        convergence_check = kwargs.get('convergence_check')
        if not 'convergence_check' in kwargs:
            convergence_check = True
        if 'figsize' in kwargs:
            figsize = kwargs.get('figsize')
        else:
            figsize = (10, 7)
        if 'sharex' in kwargs:
            sharex = kwargs.get('sharex')
        else:
            sharex = True
        fig, axes = plt.subplots(ndim, figsize=figsize, sharex=sharex)

        samples = sampler.get_chain()
        if 'labels' in kwargs.keys():
            labels = kwargs.get('labels')
        else:
            labels = list(self.params_to_estimate)
        if 'alpha' in kwargs:
            alpha = kwargs.get('alpha')
        else:
            alpha = 0.3 
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

        if 'discard' in kwargs.keys():
            discard = kwargs.get('discard')
        else:
            discard = 100 #arbitrarily discard the first 100 steps
        if 'thin' in kwargs.keys():
            thin = kwargs.get('thin')
        else:
            if convergence_check:
                thin = np.mean(np.array(self.autocorrelation_time)) / 2 #thin by half the autocorrelation time
            else:
                thin = 1
            if not np.isfinite(thin):
                thin = 1
            else:
                thin = int(thin)
        if 'flat' in kwargs.keys():
            flat = kwargs.get('flat')
        else:
            flat = True
        
        flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=flat)
        truth_list = [] 
        uncertainty_list = [] 
        # Percentiles to compute for numpy.percentile
        if 'percentiles' in kwargs.keys():
            percentiles = kwargs.get('percentiles')
        else:
            percentiles = [16, 50, 84] # Set percentiles to compute by default to q
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], q = percentiles)
            q = np.diff(mcmc)
            truth_list.append(mcmc[1])
            # uncertainty_list.append([-1.0*q[0], q[1]])
            uncertainty_list.append(q)
        try:
            import corner
            fig = corner.corner(
                flat_samples, labels=labels,
            )
        except:
            warnings.warn('corner package not found - cannot plot parameter distributions.')
        return truth_list, uncertainty_list


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
            minimizer_result[0] = self.pid_interface.get_minimizer_results(self.LL_data[0,:,:], self.timepoints, self.measurements, 
                                                                        self.initial_conditions, 
                                                                        stochastic = stochastic, debug = self.debug, 
                                                                        method = method, plot_show = plot_show, **kwargs)
        else:   
            for i in range(N):
                minimizer_result[i] = self.pid_interface.get_minimizer_results(self.LL_data[i,:,:], self.timepoints[i], 
                                                                            self.measurements, self.initial_conditions[i], 
                                                                            stochastic = stochastic, debug = self.debug, 
                                                                            method = method, plot_show = plot_show, **kwargs)
        print('Successfully completed parameter identification procedure using LMFit. Parameter values and fitness reports written to lmfit_results.csv file. Check the minimizer_results object returned to further statistically evaluate the goodness of fit.')
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

        convergence_diagnostics = kwargs.get('convergence_diagnostics') 
        if convergence_diagnostics is None:
            convergence_diagnostics = True
        if convergence_diagnostics:
            count = 0
            for result in minimizer_result:
                with open('lmfit_results.csv','a') as f:
                    # writer = csv.writer(f)
                    f.write('\nFor trajectory: {0}\n'.format(count))
                    f.write(fit_report(result))
                    count += 1
                    f.close()
