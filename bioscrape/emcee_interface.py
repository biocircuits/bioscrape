import numpy as np
try:
    import pandas as pd
except:
    print('Pandas package not found.')
import sys
import warnings
import emcee
import matplotlib.pyplot as plt
from bioscrape.types import Model
from bioscrape.sbmlutil import import_sbml as sbmlutil_import_sbml
from bioscrape.simulator import ModelCSimInterface, DeterministicSimulator, SSASimulator
from bioscrape.pid_interfaces import StochasticInference, DeterministicInference

def initialize_mcmc(**kwargs):
    obj = MCMC(**kwargs)
    return obj

class MCMC(object):
    def __init__(self, **kwargs):
        self.M = None
        if 'Model' in kwargs:
            self.set_model(kwargs.get('Model'))
        self.params_to_estimate = []
        if 'params_to_estimate' in kwargs:
            self.set_params_to_estimate(kwargs.get('params_to_estimate'))
        self.init_seed = 0.1 
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
        self.sim_type = 'stochastic'
        if 'sim_type' in kwargs:
            self.set_sim_type(kwargs.get('sim_type'))
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
        Set the dimension of parameter set to identify when using Python emcee
        '''
        self.dimension = dimension 
        return True 

    def set_init_seed(self, init_seed: float):
        '''
        Set the init_seed of parameter set to initialize Python emcee
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
        self.initial_conditions = initial_conditions 
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
        self.exp_data = exp_data 
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
        # Get initial_conditions from the model if not given explicitly
        initial_conditions = self.initial_conditions
        if initial_conditions == None: 
            self.initial_conditions = self.M.get_species_dictionary()
        self.prepare_mcmc(**kwargs)
        sampler = self.run_emcee(**kwargs)
        return sampler

    def prepare_mcmc(self, **kwargs):
        timepoints = kwargs.get('timepoints')
        exp_data = kwargs.get('exp_data')
        params = kwargs.get('params')
        prior = kwargs.get('prior')
        nwalkers = kwargs.get('nwalkers')
        init_seed = kwargs.get('init_seed')
        nsteps = kwargs.get('nsteps')
        penalty = kwargs.get('penalty')
        cost = kwargs.get('cost')
        measurements = kwargs.get('measurements')
        initial_conditions = kwargs.get('initial_conditions')

        norm_order = kwargs.get('norm_order')
        N_simulations = kwargs.get('N_simulations')
        debug = kwargs.get('debug')
        if N_simulations:
            self.N_simulations = N_simulations # Number of simulations per sample to compare to
        if norm_order:
            self.norm_order = norm_order # (integer) Which norm to use: 1-Norm, 2-norm, etc.
        if debug:
            self.debug = debug
        if type(timepoints) is list:
            if len(timepoints):
                self.timepoints = timepoints
        elif type(timepoints) is np.ndarray:
            if list(timepoints):
                self.timepoints = timepoints
        if isinstance(exp_data, (pd.DataFrame, list)):
            self.exp_data = exp_data
        if isinstance(params, list):
            self.params_to_estimate = params
        if isinstance(prior, dict):
            self.prior = prior
        if nwalkers:
            self.nwalkers = nwalkers
        if init_seed:
            self.init_seed = init_seed 
        if nsteps:
            self.nsteps = nsteps
        if penalty:
            self.penalty = penalty
        if cost:
            self.cost = cost
        if isinstance(measurements, list):
            self.measurements = measurements
        if type(initial_conditions) is dict and len(list(initial_conditions.keys())):
            self.initial_conditions = initial_conditions
        elif type(initial_conditions) is list and len(initial_conditions):
            self.initial_conditions = initial_conditions
        self.LL_data = self.extract_data(self.exp_data)

    def extract_data(self, exp_data):
        # Get timepoints from given experimental data
        if isinstance(self.timepoints, (list, np.ndarray)):
            warnings.warn('Timepoints given by user, not using the data to extract the timepoints automatically.')
        M = len(self.measurements)# Number of measurements
        # Multiple trajectories case 
        if type(self.exp_data) is list:
            data_list_final = []
            timepoints_list = []
            for df in exp_data:
                data_list = []
                if type(df) is not pd.DataFrame:
                    raise TypeError('All elements of exp_data attribute of an MCMC object must be Pandas DataFrame objects.')
                # Extract timepoints
                if self.time_column:
                    timepoint_i = np.array(df.get(self.time_column)).flatten()
                    timepoints_list.append(timepoint_i)
                else:
                    raise TypeError('time_column attribute of MCMC object must be a string.')

                # Extract measurements    
                if type(self.measurements) is list and len(self.measurements) == 1:
                    data_list.append(np.array(df.get(self.measurements[0])))
                elif type(self.measurements) is list and len(self.measurements) > 1:
                    for m in self.measurements:
                        data_list.append(np.array(df.get(m)))
                # Number of timepoints
                T = len(timepoints_list[0])
                if T != len(timepoint_i):
                    warnings.warn('The length of timepoints for all experimental trajectories must be the same, they can have different timepoints but not length of timepoints.')
                data_i = np.array(data_list)
                data_i = np.reshape(data_i, (T, M))
                data_list_final.append(data_i)
            data = np.array(data_list_final)
            self.timepoints = timepoints_list
            N = len(exp_data)# Number of trajectories
            T = len(timepoints_list[0])
            data = np.reshape(data, (N,T,M))
            if self.debug:
                print('N = {0}'.format(N))
                print('T = {0}'.format(T))
                print('M = {0}'.format(M))
                # print('The shape of data is {0}'.format(np.shape(data)))
            assert np.shape(data) == (N,T,M)
        elif type(exp_data) is pd.DataFrame:
            # Extract time
            if self.time_column:
                self.timepoints = np.array(exp_data.get(self.time_column)).flatten()
            else:
                raise TypeError('time_column attribute of MCMC object must be a string.')
            
            # Extract measurements
            if type(self.measurements) is list and len(self.measurements) == 1:
                data = np.array(exp_data.get(self.measurements[0]))
            elif type(self.measurements) is list and len(self.measurements) > 1:
                data_list = []
                for m in self.measurements:
                    data_list.append(np.array(df.get(m)))
                data = np.array(data_list)
            N = 1 # Number of trajectories
            T = len(self.timepoints) # Number of timepoints
            M = len(self.measurements)# Number of measurements
            # if self.debug:
                # print('timepoints in extract_data : {0}'.format(self.timepoints))
                # print('N = {0}'.format(N))
                # print('T= {0}'.format(T))
                # print('M = {0}'.format(M))
            data = np.reshape(data, (N,T,M))
        else:
            raise TypeError('exp_data attribute of MCMC object must be a list of Pandas DataFrames or a single Pandas DataFrame. ')
        return data
 
    def cost_function(self, params):
        if self.sim_type == 'stochastic':
            pid_interface = StochasticInference(self.params_to_estimate, self.M, self.prior)
            cost_value = pid_interface.get_likelihood_function(params, self.LL_data, self.timepoints, self.measurements, 
                                                            self.initial_conditions, norm_order = self.norm_order, 
                                                            N_simulations = self.N_simulations, debug = self.debug)
            self.cost_progress.append(cost_value)
            return cost_value
        elif self.sim_type == 'deterministic':
            pid_interface = DeterministicInference(self.params_to_estimate, self.M, self.prior)
            cost_value = pid_interface.get_likelihood_function(params, self.LL_data, self.timepoints, self.measurements,
                                                            self.initial_conditions, norm_order = self.norm_order, debug = self.debug)
            self.cost_progress.append(cost_value)
            return cost_value

    def run_emcee(self, **kwargs):
        plot_show = kwargs.get('plot_show')
        progress = kwargs.get('progress')
        if not 'progress' in kwargs:
            progress = True
        if not plot_show:
            plot_show = False
        try:
            import emcee
        except:
            print('emcee package not installed')
        ndim = len(self.params_to_estimate)
        # p0 = []
        # for walker in range(self.nwalkers):
        #     plist = []
        #     for p in self.params_to_estimate:
        #         value = self.M.get_parameter_dictionary()[p]
        #         # pinit = value + self.init_seed * np.random.randn(self.nwalkers, ndim)
        #         plist.append(pinit)
        #     p0.append(np.array(plist))   
        #     if kwargs.get('debug'):
        #         print('Sample log-like: {0}'.format(self.cost_function(np.array(plist))))
        params_values = []
        for p in self.params_to_estimate:
            value = self.M.get_parameter_dictionary()[p]
            params_values.append(value)
        p0 = np.array(params_values) + self.init_seed * np.random.randn(self.nwalkers, ndim)
        # if p0 is None:
        #     p0 = np.random.randn(ndim*self.nwalkers).reshape((self.nwalkers,self.dimension)) / 20.0
        # p0 = np.array(p0)
        assert p0.shape == (self.nwalkers, ndim)
        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.cost_function)
        sampler.run_mcmc(p0, self.nsteps, progress = progress)
        self.autocorrelation_time = sampler.get_autocorr_time()
        # for iteration, (pos,lnp,state) in enumerate(sampler.sample(p0,iterations=self.nsteps)):
            # if progress:
                # l = self.nsteps
                # printProgressBar(float(iteration), l, prefix = 'Progress:', suffix = 'Complete', length = 50)
        # Write results
        import csv
        with open('mcmc_results.csv','w', newline = "") as f:
            writer = csv.writer(f)
            writer.writerows(sampler.flatchain)
            writer.writerow('\nAutocorrelation times\n')
            writer.writerow(self.autocorrelation_time)
            writer.writerow('\nCost function progress\n')
            writer.writerow(self.cost_progress)
            f.close()
        print('Successfully completed MCMC parameter identification procedure. Parameter distribution data written to mcmc_results.csv file')
        # fitted_model, params = self.plot_mcmc_results(sampler, plot_show, **kwargs)
        return sampler
    
    def plot_mcmc_results(self, sampler, plot_show = True, **kwargs):
        print('Parameter posterior distribution convergence plots:')
        ndim = sampler.ndim
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
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=alpha)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
        axes[-1].set_xlabel("step number")

        if 'discard' in kwargs.keys():
            discard = kwargs.get('discard')
        else:
            discard = 100 #arbitrarily discard the first 100 steps
        if 'thin' in kwargs.keys():
            thin = kwargs.get('thin')
        else:
            thin = np.mean(np.array(self.autocorrelation_time)) / 2 #thin by half the autocorrelation time
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
        if 'q' in kwargs.keys():
            q = kwargs.get('q')
        else:
            q = [16, 50, 84] # Set percentiles to compute by default to q
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], q = q)
            q = np.diff(mcmc)
            # print(q)
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
        # best_p = []
        # for i in range(len(self.params_to_estimate)):
        #     my_list = [tup[i] for tup in sampler.flatchain]
        #     new_list = []
        #     for x in my_list:
        #         if x > 0:
        #             new_list.append(x)
        #     if plot_show:
        #         n, bins, patches = plt.hist(new_list, density = True, histtype = "bar")
        #         plt.title('Parameter inference distribution for parameter #{0}'.format(i))
        #     else:
        #         fig = plt.figure()
        #         n, bins, patches = plt.hist(new_list, density = True, histtype = "bar")
        #         plt.close(fig)
        #     # n, bins, patches = plt.hist(new_list, density = True, bins = 80, histtype = "bar")
        #     # Find best p
        #     best_p_ind = np.where(n == np.max(n))
        #     if np.shape(best_p_ind)[-1] != 1:
        #         warnings.warn('Multiple parameter values for {0} found with max distribution, choosing the first one. The results might be misleadig.'.format(params_to_estimate[i]))
        #         best_p_ind = np.array(best_p_ind[0].tolist()[0])
        #     assert len(best_p_ind) == 1
        #     best_p.append(bins[best_p_ind])
        #     # Plot
        #     if plot_show:
        #         plt.savefig('parameter - ' + str(self.params_to_estimate[i]) +' .svg')
        #         plt.show()

        # # Write fitted model
        # best_p = list(best_p)
        # import copy
        # fitted_model = copy.copy(self)
        # params_names = fitted_model.params_to_estimate
        # params = {}
        # for i in range(len(params_names)):
        #     p_name = params_names[i]
        #     p_sampled_value = best_p[i]
        #     params[p_name] = p_sampled_value
        # fitted_model.M.set_params(params)
        # if len(self.cost_progress) and plot_show:
        #     plt.plot(self.cost_progress, linewidth = 3)
        #     plt.title('Cost function progress')
        #     plt.xlabel('MCMC samples')
        #     plt.ylabel('Cost function')
        #     plt.savefig('Cost function progress')
        #     plt.show()
        # return fitted_model, params
    
    def simulate(self, timepoints, **kwargs):
        ''' 
        To simulate using bioscrape.
        '''
        sim_type = kwargs.get('sim_type')
        species_to_plot = kwargs.get('species_to_plot')
        plot_show = kwargs.get('plot_show')
        if not plot_show:
            plot_show = False
        if self.M:
            # If bioscrape model
            M = self.M
            s = ModelCSimInterface(M)
            if sim_type == 'deterministic':
                s.py_prep_deterministic_simulation()
                s.py_set_initial_time(timepoints[0])
                sim = DeterministicSimulator()
                result = sim.py_simulate(s, timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = M.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, M
            elif sim_type == 'stochastic':
                warnings.warn('For stochastic simulation of SBML models using bioscrape, it is highly recommended to NOT use reversible reactions as the SSA algorithm might not work for such cases.')
                sim = SSASimulator()
                s.py_set_initial_time(timepoints[0])
                result = sim.py_simulate(s,timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = M.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, M
            else:
                raise ValueError('Optional argument "sim_type" must be either deterministic or stochastic')
        else:
            raise ValueError('Model not found')

    def export_sbml(self, filename):
        raise NotImplementedError

    def import_sbml(self, filename):
        M = sbmlutil_import_sbml(filename)
        self.M = M
        return self.M

# # Print iterations progress
# def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
#     """
#     Call in a loop to create terminal progress bar
#     @params:
#         iteration   - Required  : current iteration (Int)
#         total       - Required  : total iterations (Int)
#         prefix      - Optional  : prefix string (Str)
#         suffix      - Optional  : suffix string (Str)
#         decimals    - Optional  : positive number of decimals in percent complete (Int)
#         length      - Optional  : character length of bar (Int)
#         fill        - Optional  : bar fill character (Str)
#     """
#     percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
#     filledLength = int(length * iteration // total)
#     bar = fill * filledLength + '-' * (length - filledLength)
#     print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
#     # Print New Line on Complete
#     if iteration == total: 
#         print()


class ExpData:
    def __init__(self, name, data_type, data):
        '''
        name : string representing the name of the data set 
        data_type : data_type of data set - whether time series (Use string 'timeseries') or distributional (use 'distrib')
        data : Pandas data frame
        '''
        self.name = name
        self.data_type = data_type
        self.data = data 

    def get_df(self):
        ''' 
        Returns the Pandas data frame object 
        '''
        return self.data

    def get_keys(self):
        '''
        Returns the key list of the Pandas data frame data
        '''
        return list(self.data.keys())

    def get_values(self, key):
        '''
        Returns the values as a list of the Pandas data frame object data with given key 
        '''
        return list(self.data.get(key))

 
def import_timeseries(filename, time_column, value_column, properties = {}, plot_show = False, **kwargs):
    '''
    filename : csv file with columns for data values 
    (The column numbers start at 1)
    time_column : the column number in the file that has all the time series indexes that you want to import
    value_column : the column number in the file that has all the corresponding values that you want to import 
    properties : Optional dictionary to specify other properties that the imported data must satisfy. For example, 
    properties = {3 : 'abc'}, would only impor those rows that have column 3 value equal to 'abc'
    '''
    try:
        import csv
        from operator import itemgetter
        from itertools import groupby
        import math
    except:
        print('Packages not found. Make sure csv, operator, itertool, and math are installed.')

    delimiter = kwargs.get('delimiter')
    if not delimiter:
        delimiter = ','
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter= delimiter)
        data_dict = {}
        data_dict_list = []
        for row in csv_reader:
            if row and row[time_column - 1] and row[value_column - 1]:
                if properties:
                    for col, value in properties.items():
                        if row[col - 1] == str(value):
                            data_dict[float(row[time_column - 1])] = float(row[value_column - 1])
                        else:
                            break 
                else:
                    cell_t = row[time_column - 1]
                    cell_v = row[value_column - 1]
                    temp_str_t = cell_t.replace('.','',1).replace('e','',1).replace('-','',1)
                    temp_str_v = cell_v.replace('.','',1).replace('e','',1).replace('-','',1)
                    if temp_str_t.isdigit() and temp_str_v.isdigit():
                        data_dict[float(cell_t)] = float(cell_v)
            data_dict_list.append(data_dict)
        # Create Pandas dataframe out of dictionary
        data_pd = pd.DataFrame(data_dict_list)
        data_obj = ExpData(filename, 'timeseries', data_pd)
        if plot_show:
            time = list(data_obj.get_keys())
            values = list(data_obj.get_values(data_pd.keys()))
            try:
                import matplotlib.pyplot as plt
            except:
                raise Exception('matplotlib not installed.')
            max_time = math.floor(max(time))
            index = []
            for i in range(len(time)):
                if int(math.floor(float(time[i]))) == max_time:
                    index.append(i)
            final_index = []
            for k, g in groupby(enumerate(index), lambda x:x[1]-x[0]):
                map_t = map(itemgetter(1), g)
                final_index.append(max(map_t)+1)
            init_time_index = 0
            for i in final_index:
                plt.plot(time[init_time_index:i],values[init_time_index:i])
                plt.show()
                init_time_index = i
    return data_obj
    
