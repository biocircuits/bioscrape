import numpy as np
import sys
import warnings
import emcee
import matplotlib.pyplot as plt
from bioscrape.types import Model, read_model_from_sbml
from bioscrape.simulator import ModelCSimInterface, DeterministicSimulator, SSASimulator
from bioscrape.pid_interfaces import StochasticInference, DeterministicInference

def initialize_mcmc():
    obj = MCMC()
    return obj
class MCMC(object):
    def __init__(self):
        self.M = None
        self.params_to_estimate = []
        self.prior = None
        self.nwalkers = 100
        self.nsteps = 1000
        self.nsamples = 500
        self.dimension = 0
        self.exp_data = None
        self.type = 'stochastic'
        self.timepoints = []
        self.measurements = ['']
        self.initial_conditions = None
        return 

    def get_parameters(self):
        return self.params_to_estimate

    def run_mcmc(self, **kwargs):
        self.prepare_mcmc(params = self.params_to_estimate, prior = self.prior, 
                        timepoints = self.timepoints, exp_data = self.exp_data, nwalkers = self.nwalkers, 
                        nsteps = self.nsteps, nsamples = self.nsamples, measurements = self.measurements, 
                        initial_conditions = self.initial_conditions, **kwargs)
        self.run_emcee(**kwargs)

    def prepare_mcmc(self, **kwargs):
        
        timepoints = kwargs.get('timepoints')
        exp_data = kwargs.get('exp_data')
        params = kwargs.get('params')
        prior = kwargs.get('prior')
        nwalkers = kwargs.get('nwalkers')
        nsamples = kwargs.get('nsamples')
        nsteps = kwargs.get('nsteps')
        penalty = kwargs.get('penalty')
        cost = kwargs.get('cost')
        measurements = kwargs.get('measurements')
        initial_conditions = kwargs.get('initial_conditions')

        if timepoints:
            self.timepoints = timepoints
        if exp_data.size:
            self.exp_data = exp_data
        if len(params):
            self.params_to_estimate = params
        if len(prior):
            self.prior = prior
        if nwalkers:
            self.nwalkers = nwalkers
        if nsamples:
            self.nsamples = nsamples    
        if nsteps:
            self.nsteps = nsteps
        if penalty:
            self.penalty = penalty
        if cost:
            self.cost = cost
        if len(measurements):
            self.measurements = measurements
        if type(initial_conditions) is dict and len(list(initial_conditions.keys())):
            self.MultipleInitialConditions = False 
            self.initial_conditions = initial_conditions
        elif type(initial_conditions) is list and len(initial_conditions):
            self.MultipleInitialConditions = True
            self.initial_conditions = initial_conditions
        # Create a wrapper for this to make this available to the user.
        # print(self.measurements)

    def cost_function(self, log_params):
        if self.type == 'stochastic':
            pid_interface = StochasticInference(self.params_to_estimate, self.M, self.prior)
        elif self.type == 'deterministic':
            pid_interface = DeterministicInference(self.params_to_estimate, self.M, self.prior)
        # TODO : self.Multiple .. ? needed?
        # print(pid_interface.MultipleTimepoints)
        # print(pid_interface.MultipleInitialConditions)
        # print(pid_interface.InitialConditions)
        # print(pid_interface.type)
        exp_data = self.exp_data
        timepoints = self.timepoints
        measurements = self.measurements
        initial_conditions = self.initial_conditions
        return pid_interface.get_likelihood_function(log_params, exp_data, timepoints, measurements, initial_conditions)

    def run_emcee(self, **kwargs):
        plot_show = kwargs.get('plot_show')
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
            for key, value in self.params_to_estimate.items():
                pinit = np.random.normal(value, 0.25*value)
                plist.append(pinit)
                ploglist.append(np.log(pinit))
            p0.append(np.array(plist))   
            print('Sample log-like: {0}'.format(self.cost_function(np.array(ploglist))))

        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.cost_function)
        # for i, junk in enumerate(sampler.sample(p0, iterations=nsteps)):
        #     print('Step %d' % i)
        # TODO: Add progress percentage update display code here 
        sampler.run_mcmc(p0, self.nsteps)    
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
            else:
                fig = plt.figure()
                n, bins, patches = plt.hist(new_list, density = True, histtype = "bar")
                plt.close(fig)
            # n, bins, patches = plt.hist(new_list, density = True, bins = 80, histtype = "bar")
            # Find best p
            best_p_ind = np.where(n == np.max(n))
            best_p.append(bins[best_p_ind])
            # Plot
            if plot_show:
                plt.savefig('parameter - ' + str(list(self.params_to_estimate.keys())[i]) +' .svg')
                plt.show()

        # Write fitted model
        best_p = list(best_p)
        fitted_model = self
        params_names = list(fitted_model.params_to_estimate.keys())
        params = {}
        for i in range(len(params_names)):
            p_name = params_names[i]
            p_sampled_value = best_p[i]
            params[p_name] = p_sampled_value

        fitted_model.parameters = params
        # Simulate again 
        fitted_model.simulate(self.timepoints, type = 'stochastic', species_to_plot = self.measurements, plot_show = plot_show)
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
        return 

    def import_sbml(self, filename):
        M = bioscrape.types.read_model_from_sbml(filename)
        self.M = M
        return self.M

class Data:
    def __init__(self, name, type, data = {}):
        '''
        name : string representing the name of the data set 
        type : type of data set - whether time series (Use string 'timeseries') or distributional (use 'distrib')
        data : dictionary with key and value for the data
        '''
        self.name = name
        self.type = type
        self.data = data 

    def get_keys(self):
        '''
        Returns the key list of the data dictionary
        '''
        return list(self.data.keys())

    def get_values(self):
        '''
        Returns the values list of the data dictionary
        '''
        return list(self.data.values())

 
def import_timeseries(filename, time_column, value_column, properties = {}, plot_show = False):
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

    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        data_dict = {}
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
        data_obj = Data(filename, 'timeseries', data_dict)
        time = list(data_obj.get_keys())
        values = list(data_obj.get_values())
        if plot_show:
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
    
def import_distribution(filename, index_column, value_column, properties = {}, plot_show = False):
    '''
    TODO: To be implemented
    ''' 
    return

