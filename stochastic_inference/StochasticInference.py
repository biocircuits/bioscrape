import numpy as np
import sys
import warnings
import emcee
import matplotlib.pyplot as plt
from bioscrape.types import Model, read_model_from_sbml
from bioscrape.simulator import ModelCSimInterface, DeterministicSimulator, SSASimulator

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

    
class StochasticInference(object):
    def __init__(self):
        self.m = None
        self.sbml = None
        self.params_to_estimate = {}
        self.prior = None
        self.num_walkers = 500
        self.num_iterations = 100
        self.num_samples = 100
        self.likelihoods = []
        self.dimension = 0
        self.cost = 'L2norm'
        self.penalty = 1
        self.exp_data = None
        self.timepoints = []
        self.nsamples = 50
        self.measurements = ['']
        self.species = {}
        self.parameters = {}
        return 

    def set_prior(self, prior):
        self.prior = prior
        return self.prior

    def get_parameters(self):
        return self.params_to_estimate

    def get_prior(self):
        return self.prior
   
    def log_prior(self, param_dict, prior):
        for key,value in param_dict.items():
            range = prior[key]
            if value > max(range) or value < min(range):
                return False
        return True
        
    
    def log_likelihood(self, log_params):
        measurements = self.measurements
        timepoints = self.timepoints
        nsamples = self.nsamples
        exp_data = self.exp_data
        cost = self.cost
        penalty = self.penalty
        param_dict = {}
        params_exp = np.exp(log_params)
        for key, p in zip(self.parameters.keys(),params_exp):
            param_dict[key] = p
        prior = self.get_prior()
        # Check prior
        if self.log_prior(param_dict, prior) == False:
            return -np.inf

        if self.sbml:
            try:
                import libsbml
            except:
                print('python-libsbml must be installed to use SBML models')
            filename = 'temp_simulate.xml'
            libsbml.writeSBML(self.sbml,filename)
            self.m = read_model_from_sbml(filename)
            m = self.m
            outputs = []
            for species in measurements:
                outputs.append(m.get_species_index(species))
            results = np.zeros((len(timepoints), len(outputs), nsamples))
            for sample in range(nsamples):
                sim,m  = self.simulate(timepoints, type = 'stochastic')
                for i in range(len(outputs)):
                    out = outputs[i]
                    results[:,i,sample] = sim[:,out]
        else:
            raise NotImplementedError('SBML models only (for now).')

        total_error = 0
        for i in range(len(outputs)):
            for j in range(nsamples):
                d1 = results[:,i,j]
                diff = np.abs(d1 - exp_data.get_values()) 
                '''
                TODO : Experimental data having multiple output case is important to implement here. 
                We can already handle multiple measurements. Right now what we are doing is WRONG (for multiple outputs case), because we are 
                subtracting the same experimental data from different output responses. 
                '''
                if cost == 'inf':
                    infinity_error = np.max(diff)
                    total_error += infinity_error**2
                elif cost == 'L2norm':
                    L2_norm_error = diff**2
                    L2_norm_error = np.linalg.norm(diff)
                    total_error += L2_norm_error
        # print(total_error)
        return -total_error*penalty
    

    def run_mcmc(self, params, prior, timepoints, expdata, **kwargs):

        nwalkers = kwargs.get('nwalkers')
        nsteps = kwargs.get('nsteps')
        penalty = kwargs.get('penalty')
        cost = kwargs.get('cost')
        measurements = kwargs.get('measurements')
        plot_show = kwargs.get('plot_show')

        if not params:
            params = self.params_to_estimate
        if not nwalkers:
            nwalkers = self.num_walkers
        if not nsteps:
            nsteps = self.num_iterations
        if not penalty:
            penalty = self.penalty
        if not cost:
            cost = self.cost
        if not measurements:
            measurements = self.measurements
        if not plot_show:
            plot_show = False

        # TODO: Refactor this into a new function: setup_inference that does all of this globally 
        self.params = params
        self.prior = prior
        self.timepoints = timepoints
        self.exp_data = expdata


        # Run emcee
        try:
            import emcee
        except:
            print('emcee package not installed')

        ndim = len(params)
        p0 = []
        for walker in range(nwalkers):
            plist = []
            ploglist = []
            for key, value in params.items():
                pinit = np.random.normal(value, 0.5*value)
                plist.append(pinit)
                ploglist.append(np.log(pinit))
            p0.append(np.array(plist))
            # print(p0)
        # print('going to run emcee now')    
            print('Sample log-like: {0}'.format(self.log_likelihood(np.array(ploglist))))

        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_likelihood)
        # for i, junk in enumerate(sampler.sample(p0, iterations=nsteps)):
        #     print('Step %d' % i)
        # TODO: Add progress bar here 
        sampler.run_mcmc(p0, nsteps)    
        # Write results
        import csv
        with open('mcmc_results.csv','w', newline = "") as f:
            writer = csv.writer(f)
            writer.writerows(sampler.flatchain)
            f.close()
            
        print('Successfully completed MCMC parameter identification procedure.')

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
        fitted_model.simulate(timepoints, type = 'stochastic', species_to_plot = ['c1'], plot_show = plot_show)
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
        if self.sbml:
            try:
                import libsbml
            except:
                print('libsbml-python must be installed to use SBML models')
            # If SBML model,
            filename = 'temp_simulate.xml'
            libsbml.writeSBML(self.sbml, filename) 
            m = read_model_from_sbml(filename)
            s = ModelCSimInterface(m)
            if type == 'deterministic':
                s.py_prep_deterministic_simulation()
                s.py_set_initial_time(timepoints[0])
                sim = DeterministicSimulator()
                result = sim.py_simulate(s, timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = m.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, m
            elif type == 'stochastic':
                warnings.warn('For stochastic simulation of SBML models using bioscrape, it is highly recommended to NOT use reversible reactions as the SSA algorithm might not work for such cases.')
                sim = SSASimulator()
                s.py_set_initial_time(timepoints[0])
                result = sim.py_simulate(s,timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = m.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, m
            else:
                raise ValueError('Optional argument "type" must be either deterministic or stochastic')

        elif self.m:
            # If bioscrape model
            m = self.m
            s = ModelCSimInterface(m)
            if type == 'deterministic':
                s.py_prep_deterministic_simulation()
                s.py_set_initial_time(timepoints[0])
                sim = DeterministicSimulator()
                result = sim.py_simulate(s, timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = m.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, m
            elif type == 'stochastic':
                warnings.warn('For stochastic simulation of SBML models using bioscrape, it is highly recommended to NOT use reversible reactions as the SSA algorithm might not work for such cases.')
                sim = SSASimulator()
                s.py_set_initial_time(timepoints[0])
                result = sim.py_simulate(s,timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = m.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, m
            else:
                raise ValueError('Optional argument "type" must be either deterministic or stochastic')

    def export_sbml(self, filename):
        try:
            import libsbml
        except:
            print('libsbml-python must be installed to use SBML models')
        if self.sbml:
            model = self.sbml.getModel()
            params = self.params_to_estimate
            for pid,pval in params.items():
                if isinstance(pval, (list, np.ndarray)):
                    pval = pval[0]
                model.getElementBySId(pid).setValue(float(pval))
            libsbml.writeSBML(self.sbml,filename)
        elif self.m:
            model_doc = self.m.convert_to_sbml()

            libsbml.writeSBML(model_doc, filename)
        else:
            raise ValueError('Model must be SBML or bioscrape XML. Other models not supported.')

    def import_sbml(self, filename):
        try:
            import libsbml
        except:
            print('Unable to import libsbml. Make sure that python-libsbml package is installed and working')
        reader = libsbml.SBMLReader()
        doc = reader.readSBML(filename)
        self.sbml = doc
        params_dict = {}
        species_dict = {}
        # Get global parameters
        for param in doc.getModel().getListOfParameters():
            params_dict[param.getId()] = param.getValue()
        count = 0
        # Get local parameters
        for reaction in doc.getModel().getListOfReactions():
                for param in reaction.getKineticLaw().getListOfLocalParameters():
                    count = count + 1
                    params_dict[param.getId() + '_local_' + reaction.getId() + str(count)] = param.getValue()
        self.parameters = params_dict

        for species in doc.getModel().getListOfSpecies():
            species_dict[species.getId()] = species.getInitialAmount()

        self.species = species_dict
        return self.sbml

