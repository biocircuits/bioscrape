from bioscrape.inference import StochasticInference

spid = StochasticInference()

# # Import a model that you want to use for inference
# spid.import_sbml('sbml_example.xml')
# # Supports SBML rate rules (TODO: Implement this), initial assignments, 
# # and other SBML level 3 core features (packages not supported)
# # (Optional) 
# model = spid.import_bioscrape('bs.xml')
# # (Optional) model = import_ode(ode_function) TODO: Not implemented yet
M = Model('test.xml') 
#########################################  Part - I  ######################################### 

# Handling the data

# properties_dict = {'Property1': 'Value1','Property2': 'Value2'}
# data = import_timeseries('data.csv',time_column = 1, value_column = 2, properties = properties_dict)
data = import_timeseries('test_data.csv', time_column = 2, value_column = 4, properties = {3 : 51})

# Plot experimental data, give a True value to the optional argument plot_show
# import_timeseries('test_data.csv', time_column = 2, value_column = 4, plot_show = True)

# Initial simulation data
import numpy as np
import matplotlib.pyplot as plt  

t_exp = data.get_keys()
timepoints = np.linspace(t_exp[0],t_exp[ -1],len(t_exp))
bs_data, m = spid.simulate(timepoints, type = 'deterministic')
simtime = timepoints
simdata = bs_data[:,m.get_species_index('c1')]
# simtime, simdata = model.simulate_roadrunner(timepoints)

# Create artificial data for model
# SKIP if you already have your data, store it inside exp_data

lines = []
exp_value = []
for t,d in zip(simtime, simdata):
    v = d + np.random.normal(0,0.5)
    exp_value.append(v)
    lines.append([t,v])
import csv 
with open('test_data_noisy.csv','w') as f:
    writer = csv.writer(f)
    writer.writerows(lines)
    f.close()

plt.plot(simtime, exp_value, label = 'Data')
plt.plot(simtime, simdata, label = 'Model simulation')
plt.legend()
plt.title('Simulated model and artificial (noisy) data')
plt.show()


#########################################  Part - II  ######################################### 

# Identifiability and sensitivity analysis (NotImplemented)
data = import_timeseries('test_data_noisy.csv', time_column = 1, value_column = 2)
#Analysis
list_of_measurements = ['c1'] #outputs being measured
list_of_timepoints = [t_exp] # time range for which the output is measured
list_of_intervals = [0.1] # time interval between each measurements above
# initial_params = {'kc' : 0.6, 'k1' : 1} #optionally set initial guesses for parameters here
# initial_species = {'A' : 500} # optionally set initial conditions for species
'''
TODO: Implement analysis_sens for sensitivity analysis
TODO: Implement analysis_ident for parameter identifiability analysis
'''
#### model.analysis_sens(list_of_measurements, list_of_timepoints, list_of_intervals, initial_params, initial_species)
# Returns parameters that are identifiable
#### params = model.analysis_ident(list_of_measurements, list_of_timepoints, list_of_intervals, initial_params, initial_species)

#########################################  Part - III  ######################################### 

# Parameter identification using MCMC

spid.prior = {'kc' : [1e-3, 1e3],'k1' : [1e-2, 1e5]}
spid.params_to_estimate = {'kc':6, 'k1':1}

# Prepare for MCMC and run
spid.prepare_mcmc(params = spid.params_to_estimate, prior = spid.prior, timepoints = timepoints, 
                  exp_data = data, nwalkers = 40, nsteps = 100, nsamples = 500, measurements = ['c1'])
fit_model, id_params = spid.run_mcmc(plot_show = True)
# fit_model is with identified parameters' "best" value substituted


#########################################  Part - IV  ######################################### 

# Post identification analysis

res_orig, m_orig = spid.simulate(timepoints, type = 'stochastic')
res, m = spid.simulate(timepoints, type = 'stochastic')

plt.figure()
simtime = timepoints
simdata_orig = res_orig[:,m_orig.get_species_index('c1')]
simdata = res[:,m.get_species_index('c1')]
plt.plot(simtime, exp_value, label ='exp')
plt.plot(simtime, simdata_orig, label = 'original model')
plt.plot(simtime, simdata, label = 'identified model')
plt.legend()
plt.title('Identified model with data')
plt.show()


# Write the fitted model to SBML/bioscrape
fit_model.export_sbml('sbml_fit.xml')
# Optionally, write a bioscrape file.
# fit_model.export_bioscrape('bs_fit.xml')

# TODO list:
# 1. Give option to create custom log-likelihood function and use with run_mcmc. - Done.
# 2. Fix the multiple output case in the built-in log-likelihood function. - Done. Works for single output. Do we need multiple outputs?
# 3. Include setup_inference for easier installation. 
# 4. Setup StochasticInference module along with bioscrape installation or provide a separate installtion command like we have with lineages
# 5. Fix all other TODO at other places. - Finished most.
# 6. Add bioscrape Model compatibility as the central model maybe and then just have one place where you import sbml (currently it's a bit messed up)
# 7. Finish integration with bioscrape.inference somehow.
# 8. Other kinds of bioscrape compatibilities and functionalities (that is, do things in the bioscrape-way)