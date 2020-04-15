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

def initialize_mcmc():
    obj = MCMC()
    return obj
class MCMC(object):
    def __init__(self):
        self.M = None
        self.params_to_estimate = []
        self.init_seed = 0.1 
        self.prior = None
        self.nwalkers = 100
        self.nsteps = 200
        self.dimension = 0
        self.exp_data = None # Pandas DataFrame object
        self.type = 'stochastic'
        self.timepoints = None
        self.time_column = ''
        self.measurements = ['']
        self.initial_conditions = None
        self.norm_order = 2
        self.N_simulations = 3
        self.LL_data = None
        self.debug = False
        return 

    def get_parameters(self):
        return self.params_to_estimate

    def run_mcmc(self, **kwargs):
        # Get initial_conditions from the model if not given explicitly
        initial_conditions = self.initial_conditions
        if initial_conditions == None: 
            self.initial_conditions = self.M.get_species_dictionary()
        self.prepare_mcmc(**kwargs)
        fitted_model, params = self.run_emcee(**kwargs)
        return fitted_model.M, params

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
                print('M = {0}'.format(M))
                print('T = {0}'.format(T))
                print('The shape of data is {0}'.format(np.shape(data)))
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
            if self.debug:
                # print('timepoints in extract_data : {0}'.format(self.timepoints))
                print('N = {0}'.format(N))
                print('M = {0}'.format(M))
                print('T = {0}'.format(T))
            data = np.reshape(data, (N,T,M))
        else:
            raise TypeError('exp_data attribute of MCMC object must be a list of Pandas DataFrames or a single Pandas DataFrame. ')
        return data
 
    def cost_function(self, log_params):
        if self.type == 'stochastic':
            pid_interface = StochasticInference(self.params_to_estimate, self.M, self.prior)
            return pid_interface.get_likelihood_function(log_params, self.LL_data, self.timepoints, self.measurements, self.initial_conditions, norm_order = self.norm_order, N_simulations = self.N_simulations, debug = self.debug)
        elif self.type == 'deterministic':
            pid_interface = DeterministicInference(self.params_to_estimate, self.M, self.prior)
            return pid_interface.get_likelihood_function(log_params, self.LL_data, self.timepoints, self.measurements, self.initial_conditions, norm_order = self.norm_order, debug = self.debug)

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
        p0 = []
        for walker in range(self.nwalkers):
            plist = []
            ploglist = []
            for p in self.params_to_estimate:
                value = self.M.get_parameter_dictionary()[p]
                pinit = np.random.normal(value, self.init_seed*value)
                plist.append(pinit)
                ploglist.append(np.log(pinit))
            p0.append(np.array(plist))   
            if progress and kwargs.get('debug'):
                print('Sample log-like: {0}'.format(self.cost_function(np.array(ploglist))))

        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.cost_function)
        if p0 is None:
            p0 = np.random.randn(ndim*self.nwalkers).reshape((self.nwalkers,self.dimension)) / 20.0

        for iteration, (pos,lnp,state) in enumerate(sampler.sample(p0,iterations=self.nsteps)):
            if progress:
                l = self.nsteps
                printProgressBar(float(iteration), l, prefix = 'Progress:', suffix = 'Complete', length = 50)
        # Write results
        import csv
        with open('mcmc_results.csv','w', newline = "") as f:
            writer = csv.writer(f)
            writer.writerows(sampler.flatchain)
            f.close()
        print('Successfully completed MCMC parameter identification procedure. Parameter distribution data written to mcmc_results.csv file')
        fitted_model, params = self.plot_mcmc_results(sampler, plot_show)
        return fitted_model, params
    
    def plot_mcmc_results(self, sampler, plot_show = True):
        best_p = []
        for i in range(len(self.params_to_estimate)):
            my_list = [tup[i] for tup in sampler.flatchain]
            new_list = []
            for x in my_list:
                if x > 0:
                    new_list.append(x)
            if plot_show:
                n, bins, patches = plt.hist(new_list, density = True, histtype = "bar")
                plt.title('Parameter inference distribution for parameter #{0}'.format(i))
            else:
                fig = plt.figure()
                n, bins, patches = plt.hist(new_list, density = True, histtype = "bar")
                plt.close(fig)
            # n, bins, patches = plt.hist(new_list, density = True, bins = 80, histtype = "bar")
            # Find best p
            best_p_ind = np.where(n == np.max(n))
            if np.shape(best_p_ind)[-1] != 1:
                warnings.warn('Multiple parameter values for {0} found with max distribution, choosing the first one. The results might be misleadig.'.format(params_names[i]))
                best_p_ind = np.array(best_p_ind[0].tolist()[0])
            assert len(best_p_ind) == 1
            best_p.append(bins[best_p_ind])
            # Plot
            if plot_show:
                plt.savefig('parameter - ' + str(self.params_to_estimate[i]) +' .svg')
                plt.show()

        # Write fitted model
        best_p = list(best_p)
        import copy
        fitted_model = copy.copy(self)
        params_names = fitted_model.params_to_estimate
        params = {}
        for i in range(len(params_names)):
            p_name = params_names[i]
            p_sampled_value = best_p[i]
            params[p_name] = p_sampled_value
        fitted_model.M.set_params(params)
        return fitted_model, params
    
    def simulate(self, timepoints, **kwargs):
        ''' 
        To simulate using bioscrape.
        '''
        type = kwargs.get('type')
        species_to_plot = kwargs.get('species_to_plot')
        plot_show = kwargs.get('plot_show')
        if not plot_show:
            plot_show = False
        if self.M:
            # If bioscrape model
            M = self.M
            s = ModelCSimInterface(M)
            if type == 'deterministic':
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
            elif type == 'stochastic':
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
                raise ValueError('Optional argument "type" must be either deterministic or stochastic')
        else:
            raise ValueError('Model not found')

    def export_sbml(self, filename):
        raise NotImplementedError

    def import_sbml(self, filename):
        M = sbmlutil_import_sbml(filename)
        self.M = M
        return self.M

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()


class ExpData:
    def __init__(self, name, type, data):
        '''
        name : string representing the name of the data set 
        type : type of data set - whether time series (Use string 'timeseries') or distributional (use 'distrib')
        data : Pandas data frame
        '''
        self.name = name
        self.type = type
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
    
