cimport numpy as np
import numpy as np
from libc.math cimport log

from types import Model
from types cimport Model
import sys

import emcee

##################################################                ####################################################
######################################              DISTRIBUTION                      ################################
#################################################                     ################################################


cdef class Distribution:
    cdef double unprob(self,np.ndarray a):
        raise NotImplementedError("Calling unprob() on Distribution")
    cdef double prob(self,np.ndarray a):
        raise NotImplementedError("Calling prob() on Distribution")
    cdef unsigned dim(self):
        raise NotImplementedError("Calling dim() on Distribution")

    def py_unprob(self, np.ndarray a):
        return self.unprob(a)

    def py_prob(self, np.ndarray a):
        return self.prob(a)

    def py_dim(self):
        return self.dim()


cdef class UniformDistribution(Distribution):
    def __init__(self, np.ndarray lb, np.ndarray ub):
        self.dimension = lb.shape[0]
        self.lower_bounds = lb.copy()
        self.upper_bounds = ub.copy()
        self.volume = 1.0

        if lb.shape[0] != ub.shape[0]:
            raise SyntaxError("Wrong bound array sizes")

        cdef unsigned i
        for i in range(self.dimension):
            if self.lower_bounds[i] > self.upper_bounds[i]:
                raise SyntaxError("Lower bound > Upper bound")
            self.volume *= (self.upper_bounds[i]-self.lower_bounds[i])

    cdef double unprob(self,np.ndarray[np.double_t, ndim=1] a):
        cdef unsigned i
        for i in range(self.dimension):
            if a[i] > self.upper_bounds[i] or a[i] < self.lower_bounds[i]:
                return 0

        return 1

    cdef double prob(self,np.ndarray a):
        cdef unsigned i
        for i in range(self.dimension):
            if a[i] > self.upper_bounds[i] or a[i] < self.lower_bounds[i]:
                return 0

        return 1/self.volume


    cdef unsigned dim(self):
        return self.dimension


##################################################                ####################################################
######################################              DATA                              ################################
#################################################                     ################################################

cdef class BulkData:
    def __init__(self):
        self.measured_species = []
        self.timepoints = None
        self.measurements = None

    def set_data(self,np.ndarray timepoints, np.ndarray measurements, list measured_species):
        self.timepoints = timepoints
        self.measurements = measurements
        self.measured_species = measured_species


cdef class Data:
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

    def import_timeseries(self, filename, time_column, value_column, properties = {}, plot_show = False):
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
    
    def import_distribution(self, filename, index_column, value_column, properties = {}, plot_show = False):
        '''
        TODO: To be implemented
        ''' 
        return
##################################################                ####################################################
######################################              LIKELIHOOD                        ################################
#################################################                     ################################################

cdef class Likelihood:
    cdef double get_log_likelihood(self):
        raise NotImplementedError("Calling get_log_likelihood() for Likelihood")

cdef class DeterministicLikelihood(Likelihood):
    def __init__(self):
        self.m = None
        self.init_state_vals = np.zeros(0)
        self.init_param_vals = np.zeros(0)
        self.init_param_indices = np.zeros(0,dtype=int)
        self.init_state_indices = np.zeros(0,dtype=int)
        self.csim = None
        self.propagator = None


    def set_model(self, Model m, CSimInterface csim = None, RegularSimulator prop = None):
        self.m = m
        self.csim = csim
        self.propagator = prop

        if csim is None:
            csim = ModelCSimInterface(m)
        if prop is None:
            prop = DeterministicSimulator()


    def set_data(self, BulkData bd):
        self.bd = bd
        the_list = bd.get_measured_species()
        self.meas_indices = np.zeros(len(the_list), dtype=int)
        for i in range(len(the_list)):
            self.meas_indices[i] = self.m.get_species_index( the_list[i] )


    def set_init_species(self, dict sd):
        pairs = sd.items()
        self.init_state_indices = np.zeros(len(pairs), dtype=int)
        self.init_state_vals = np.zeros(len(pairs),)

        index = 0
        for (key,val) in  pairs:
            self.init_state_indices[index] = self.m.get_species_index( key  )
            self.init_state_vals[index] = val
            index +=  1


    def set_init_params(self, dict pd):
        pairs = pd.items()
        self.init_param_indices = np.zeros(len(pairs), dtype=int)
        self.init_param_vals = np.zeros(len(pairs),)

        index = 0
        for (key,val) in pairs:
            self.init_param_indices[index] = self.m.get_param_index( key  )
            self.init_param_vals[index] = val
            index +=  1


    cdef double get_log_likelihood(self):
        # Write in the specific parameters and species values.
        cdef np.ndarray species_vals = self.m.get_species_values()
        cdef np.ndarray param_vals = self.m.get_params_values()

        cdef unsigned i
        for i in range(self.init_state_indices.shape[0]):
            species_vals[ self.init_state_indices[i] ] = self.init_state_vals[i]
        for i in range(self.init_param_indices.shape[0]):
            param_vals[ self.init_param_indices[i] ] = self.init_param_vals[i]

        # Do a simulation of the model with time points specified by the data.
        cdef np.ndarray ans = self.propagator.simulate(self.csim, self.bd.get_timepoints()).get_result()

        # Compare the data using l_1 norm and return the likelihood.
        cdef double error = 0.0
        for i in range(self.meas_indices.shape[0]):
            error += np.linalg.norm(  self.bd.get_measurements()[:,i] - ans[:,self.meas_indices[i]] )

        return -error

##################################################                ####################################################
######################################              INFERENCE                         ################################
#################################################                     ################################################


cdef class DeterministicInference:
    def __init__(self):
        self.m = None
        self.params_to_estimate = np.zeros(0,dtype=int)
        self.global_init_params = np.zeros(0,)
        self.global_init_state  = np.zeros(0,)
        self.prior = None
        self.num_walkers = 500
        self.num_iterations = 100
        self.likelihoods = []
        self.dimension = 0
        self.sigma = 10

    def set_model(self, Model m):
        self.m = m
        self.global_init_params = m.get_params_values().copy()
        self.global_init_state = m.get_species_values().copy()

    def set_sigma(self, double s):
        self.sigma = s

    def set_mcmc_params(self, unsigned walkers, unsigned iterations):
        self.num_walkers = walkers
        self.num_iterations = iterations

    def set_prior(self, Distribution prior):
        self.prior = prior

    def set_params_to_estimate(self, list params):
        self.params_to_estimate = np.zeros(len(params),dtype=int)
        cdef unsigned i
        for i in range(len(params)):
            self.params_to_estimate[i] =  self.m.get_param_index(params[i])

        self.dimension = len(params)

    def add_to_likelihood(self,list likelihoods):
        self.likelihoods.extend(likelihoods)

    def likelihood_function(self, np.ndarray params):
        cdef unsigned i
        cdef unsigned j
        cdef double answer = self.prior.unprob(params)

        if answer == 0.0:
            return -np.inf

        answer = log(answer)

        cdef np.ndarray actual_species = self.m.get_species_values()
        cdef np.ndarray actual_params = self.m.get_params_values()

        for i in range(len(self.likelihoods)):
            # Copy in the global initial conditions.
            np.copyto(actual_species, self.global_init_state)
            np.copyto(actual_params, self.global_init_params)

            # Then copy in the specified parameters
            for j in range(self.params_to_estimate.shape[0]):
                actual_params[ self.params_to_estimate[j] ]  = params[j]

            # The run the likelihood
            answer += (<Likelihood> self.likelihoods[i]).get_log_likelihood()

        return answer / self.sigma**2

    def py_likelihood_function(self, np.ndarray params):
        return self.likelihood_function(params)

    def run_mcmc(self, p0 = None):
        def lnprob(np.ndarray theta):
            loglhood = self.likelihood_function(np.exp(theta))
            if np.isnan(loglhood):
                sys.stderr.write('Parameters returned NaN likelihood: ' + str(np.exp(theta)) + '\n')
                sys.stderr.flush()
                return -np.Inf
            return loglhood

        sampler = emcee.EnsembleSampler(self.num_walkers, self.dimension, lnprob)
        if p0 is None:
            p0 = np.random.randn(self.dimension*self.num_walkers).reshape((self.num_walkers,self.dimension)) / 20.0

        for iteration, (pos,lnp,state) in enumerate(sampler.sample(p0,iterations=self.num_iterations)):
            print('%.1f percent complete' % (100*float(iteration)/self.num_iterations))

        return sampler


cdef class StochasticInference:
    def __init__(self):
        self.m = None
        self.params_to_estimate = np.zeros(0, dtype=int)
        self.global_init_params = np.zeros(0,)
        self.global_init_state  = np.zeros(0,)
        self.prior = None
        self.num_walkers = 500
        self.num_iterations = 100
        self.num_samples = 100
        self.likelihoods = []
        self.dimension = 0
        self.sigma = 10

    def set_model(self, Model m):
        self.m = m
        self.global_init_params = m.get_params_values().copy()
        self.global_init_state = m.get_species_values().copy()

    def set_sigma(self, double s):
        self.sigma = s

    def set_mcmc_params(self, unsigned walkers, unsigned iterations):
        self.num_walkers = walkers
        self.num_iterations = iterations

    def set_prior(self, Distribution prior):
        self.prior = prior

    def set_params_to_estimate(self, list params):
        self.params_to_estimate = np.zeros(len(params),dtype=int)
        cdef unsigned i
        for i in range(len(params)):
            self.params_to_estimate[i] =  self.m.get_param_index(params[i])
        self.dimension = len(params)

    def add_to_likelihood(self,list likelihoods):
        self.likelihoods.extend(likelihoods)

    def likelihood_function(self, np.ndarray params):
        cdef unsigned i
        cdef unsigned j
        cdef double answer = self.prior.unprob(params)

        if answer == 0.0:
            return -np.inf

        answer = log(answer)

        cdef np.ndarray actual_species = self.m.get_species_values()
        cdef np.ndarray actual_params = self.m.get_params_values()

        for i in range(len(self.likelihoods)):
            # Copy in the global initial conditions.
            np.copyto(actual_species, self.global_init_state)
            np.copyto(actual_params, self.global_init_params)

            # Then copy in the specified parameters
            for j in range(self.params_to_estimate.shape[0]):
                actual_params[ self.params_to_estimate[j] ]  = params[j]

            # The run the likelihood
            answer += (<Likelihood> self.likelihoods[i]).get_log_likelihood()

        return answer / self.sigma**2

    def py_likelihood_function(self, np.ndarray params):
        return self.likelihood_function(params)
    
    def log_prior(param_dict, priors):
        for key,value in param_dict.items():
            range = priors[key]
            if value > max(range) or value < min(range):
                return False
        return True
        
    def log_likelihood(log_params):
        param_dict = {}
        params_exp = np.exp(log_params)
        for key, p in zip(self.parameters.keys(),params_exp):
            param_dict[key] = p
        # Check prior
        if log_prior(param_dict, priors) == False:
            return -np.inf

        if simulator == 'bioscrape':
            try:
                import bioscrape
            except:
                print('bioscrape package must be installed to use bioscrape simulator')
            if self.sbml:
                try:
                    import libsbml
                except:
                    print('python-libsbml must be installed to use SBML models')
                filename = 'temp_simulate.xml'
                libsbml.writeSBML(self.sbml,filename)
                self.bioscrape = bioscrape.types.read_model_from_sbml(filename)
            outputs = []
            for species in measurements:
                outputs.append(self.bioscrape.get_species_index(species))
            results = np.zeros((len(timepoints), len(outputs), nsamples))
            for sample in range(nsamples):
                sim, m = self.simulate_bioscrape(timepoints, sim_type)
                for i in range(len(outputs)):
                    out = outputs[i]
                    results[:,i,sample] = sim[:,out]
        else:
            raise ValueError('Other simulators not implemented.')

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
        print(total_error)
        return -total_error*penalty
    

    def run_mcmc(self, params, priors, timepoints, measzurements, Data exp_data,
                cost = 'L2norm', penalty = 1, plot_show = False):

        params = self.params_to_estimate
        nwalkers = self.num_walkers
        nsteps = self.num_iterations
        nsamples = self.num_samples

        '''
        cost: L2norm for L2 norm or inf for infinity norm cost
        '''
        # Run emcee
        try:
            import emcee
        except:
            print('emcee package not installed')
        ndim = len(self.parameters)
        p0 = []
        for walker in range(nwalkers):
            plist = []
            ploglist = []
            for key, value in self.parameters.items():
                pinit = np.random.normal(value, 0.5*value)
                plist.append(pinit)
                ploglist.append(np.log(pinit))
            p0.append(np.array(plist))
            # print(p0)
        # print('going to run emcee now')    
            print('Sample log-like: {0}'.format(log_likelihood(np.array(ploglist))))

        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood)
        # for i, junk in enumerate(sampler.sample(p0, iterations=nsteps)):
        #     print('Step %d' % i)
        sampler.run_mcmc(p0, nsteps)    
        # Write results
        import csv
        with open('mcmc_results.csv','w', newline = "") as f:
            writer = csv.writer(f)
            writer.writerows(sampler.flatchain)
            # f.write(str(sampler.flatchain))
            f.close()
            
        print('Successfully completed MCMC parameter identification procedure.')

        best_p = []
        for i in range(len(self.parameters)):
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
                plt.savefig('parameter - ' + str(list(self.parameters.keys())[i]) +' .svg')
                plt.show()

        # Write fitted model
        best_p = list(best_p)
        fitted_model = self
        params_names = list(fitted_model.parameters.keys())
        params = {}
        for i in range(len(params_names)):
            p_name = params_names[i]
            p_sampled_value = best_p[i]
            params[p_name] = p_sampled_value

        fitted_model.parameters = params

        # Simulate again
        fitted_model.simulate_bioscrape(timepoints, 'stochastic', species_to_plot = ['c1'], plot_show = plot_show)
        return fitted_model, params


    def run_mcmc(self, p0 = None):
        def lnprob(np.ndarray theta):
            loglhood = self.likelihood_function(np.exp(theta))
            if np.isnan(loglhood):
                sys.stderr.write('Parameters returned NaN likelihood: ' + str(np.exp(theta)) + '\n')
                sys.stderr.flush()
                return -np.Inf
            return loglhood

        sampler = emcee.EnsembleSampler(self.num_walkers, self.dimension, lnprob)
        if p0 is None:
            p0 = np.random.randn(self.dimension*self.num_walkers).reshape((self.num_walkers,self.dimension)) / 20.0

        for iteration, (pos,lnp,state) in enumerate(sampler.sample(p0,iterations=self.num_iterations)):
            print('%.1f percent complete' % (100*float(iteration)/self.num_iterations))

        return sampler
