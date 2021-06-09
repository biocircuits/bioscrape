
import warnings
from bioscrape.types import Model
from bioscrape.simulator import ModelCSimInterface, DeterministicSimulator
from scipy.integrate import odeint
import numpy as np

def py_sensitivity_analysis(model, timepoints, normalize, **kwargs):
    sens_obj = SensitivityAnalysis(model)
    ans_df = sens_obj.propagator.py_simulate(sens_obj.sim_interface, timepoints).py_get_dataframe(sens_obj.M)
    solutions_array = np.array(ans_df.iloc[:,range(0,len(ans_df.T) - 1)])
    return sens_obj.compute_SSM(solutions_array, timepoints, normalize, **kwargs)

def py_get_jacobian(model, state, **kwargs):
    return SensitivityAnalysis(model).compute_J(state, **kwargs)
def py_get_sensitivity_to_parameter(model, state, param_name, **kwargs):
    return SensitivityAnalysis(model).compute_Zj(state, param_name, **kwargs)

class SensitivityAnalysis(Model):
    def __init__(self, M, dx = 0.01):
        """
        Local Sensitivity Analysis for Bioscrape models.
        Arguments:
        * M: The Bioscrape Model object.
        * dx: Small parameter used in approximate computation methods. 
        """
        self.M = M
        sim = ModelCSimInterface(self.M)
        sim.py_prep_deterministic_simulation()
        self.sim_interface = sim
        self.propagator = DeterministicSimulator()
        self.num_equations = sim.py_get_num_species()
        self.dx = 0.01
        self.original_parameters = dict(M.get_parameter_dictionary())
    
    def _evaluate_model(self, states, params = None, time = 0.0):
        """
        Numerically evaluates the model at a given value for the states and time (if time-varying model).
        """
        sim = self.sim_interface
        if params is not None:
            self.M.set_params(params)
        states = np.array(states, dtype = 'float64')
        derivative_array = np.zeros((self.num_equations), dtype = 'float64')
        # add apply_repeated_rules
        sim.py_calculate_deterministic_derivative(states, derivative_array, time)
        return derivative_array

    def compute_J(self, x, **kwargs):
        """
        Compute the Jacobian J = df/dx at a point x.
        Returns a matrix of size n x n.
        Uses fourth-order central difference method to compute Jacobian
        """
        method = kwargs.get('method')
        if method is None:
            method = 'fourth_order_central_difference'
        x = np.array(x, dtype = 'float64')
        n = len(x)
        # initialize J
        J = np.zeros( (n, n) )   
        # Future: Use numdifftools to compute Jacobian matrix
        if method == 'numdifftools':
            # jself = nd.Jacobian(lambda x: self_ode(0, x, P), **kwargs) # Using numdifftools
            # return jself(x)
            return NotImplementedError
        # store the variable with respect to which we approximate the differentiation (df/dvar)
        state_input = np.array(x)
        for i in range(n):
            for j in range(n):
                h = state_input[j]*self.dx
                # h = self.dx
                # Gets O(h^4) central difference on df_i/dvar_j
                if h != 0:
                    x = np.array(state_input)
                    f_0 = self._evaluate_model(x)[i]
                    x = np.array(state_input)
                    x[j] = x[j] + 2*h
                    f_2h = self._evaluate_model(x)[i]
                    x = np.array(state_input)
                    x[j] = x[j] + h
                    f_h = self._evaluate_model(x)[i]
                    x = np.array(state_input)
                    x[j] = x[j] - h
                    f_mh = self._evaluate_model(x)[i]
                    x = np.array(state_input)
                    x[j] = x[j] - 2*h
                    f_m2h = self._evaluate_model(x)[i]
                    x = np.array(state_input)
                    if method == 'fourth_order_central_difference':
                        # compute these only when needed
                        J[i,j]= (-f_2h + 8*f_h - 8*f_mh + f_m2h)/(12*h)
                    if method == 'central_difference':
                        J[i,j]= (f_h - f_mh)/(2*h) 
                    if method == 'backward_difference':
                        J[i,j]= (f_0 - f_mh)/h
                    if method == 'forward_difference':
                        J[i,j]= (f_h - f_0)/h
                    # Error check
                    if J[i, j] == np.Inf:
                        warnings.warn('Inf found while computing the Jacobian.')
                        J[i, j] = 1
                    elif J[i, j] == np.NaN:
                        warnings.warn('NaN found while conputing the Jacobian.')
                        J[i, j] = 0
        return J
        
    def compute_Zj(self, x, param_name, **kwargs):
        """
        Compute Z_j, i.e. df/dparam_name at a particular point x
        Returns a vector of size n x 1. 
        """
        method = kwargs.get('method')
        if method is None:
            method = 'fourth_order_central_difference'
        x = np.array(x, dtype = 'float64')
        n = len(x)
        Z = np.zeros(n)    
        params_dict = dict(self.M.get_parameter_dictionary())
        # For each state
        for i in range(n):
            h = params_dict[param_name]*0.01 # Small parameter for this parameter
            # Gets O(4) central difference on dfi/dpj
            if h != 0:
                params_dict = self.original_parameters
                self.M.set_params(params_dict)
                params_dict[param_name] = params_dict[param_name] + 2*h
                f_2h = self._evaluate_model(x, params_dict)[i]

                params_dict = self.original_parameters
                self.M.set_params(params_dict)
                params_dict[param_name] = params_dict[param_name] + h
                f_h = self._evaluate_model(x, params_dict)[i]

                params_dict = self.original_parameters
                self.M.set_params(params_dict)
                params_dict[param_name] = params_dict[param_name] - h
                f_mh = self._evaluate_model(x, params_dict)[i]

                params_dict = self.original_parameters
                self.M.set_params(params_dict)
                params_dict[param_name] = params_dict[param_name] - 2*h
                f_m2h = self._evaluate_model(x, params_dict)[i]
                #Store approx. dfi/dp[param_name] into Z
                if method == 'fourth_order_central_difference':
                    Z[i]= (-f_2h + 8*f_h - 8*f_mh + f_m2h)/(12*h)
                if method == 'central_difference':
                    Z[i]= (f_h - f_mh)/(2*h) 
                if method == 'backward_difference':
                    Z[i]= (f_0 - f_mh)/h
                if method == 'forward_difference':
                    Z[i]= (f_h - f_0)/h
                # Error check
                if Z[i] == np.Inf:
                    warnings.warn('Inf found')
                    Z[i] = 1
                elif Z[i] == np.NaN:
                    warnings.warn('NaN found')
                    Z[i] = 0
        return Z

    def compute_SSM(self, solutions, timepoints, normalize = False, **kwargs):
        """
        Returns the sensitivity coefficients S_j for each parameter p_j. 
        Solutions is the ODE solution to self for timepoints.
        Solutions is of shape (len(timepoints), n), where n is the len(x).
        The sensitivity coefficients are written in a sensitivity matrix SSM of size len(timepoints) x len(params) x n
        If normalize argument is true, the coefficients are normalized by the nominal value of each paramneter.
        Arguments:
        * solutions: Pandas dataframe object returned by py_simulate_model that contains solutions for all model variables.
        * timepoints: The time points at which sensitivity coefficients are needed (this is the same as timepoints used for solutions).
        * normalize: (bool, default is False): When set to True, the returned sensitivity coefficients are normalized with state and parameter values.
        * kwargs: Other kwargs passed to `compute_J` and `compute_Z` functions.
        """
        def sensitivity_ode(t, x, J, Z):
            # ODE to solve for sensitivity coefficient S
            dsdt = J@x + Z
            return dsdt
        all_params = self.M.get_param_list()
        number_of_params = len(all_params)
        n = self.num_equations
        S0 = np.zeros(n) # Initial value for S_i  
        SSM = np.zeros( (len(timepoints), number_of_params, n) )
        xs = solutions
        xs = np.reshape(xs, (len(timepoints), n) )
        # Solve for SSM at each time point 
        for k in range(len(timepoints)): 
            timepoints_ssm = timepoints[0:k+1]
            if len(timepoints_ssm) == 1:
                continue
            # get the jacobian matrix
            J = self.compute_J(xs[k,:], **kwargs)
            #Solve for S = dx/dp for all x and all P (or theta, the parameters) at time point k
            for j in range(len(all_params)):
                param_name = all_params[j]
                # get the pmatrix
                Zj = self.compute_Zj(xs[k,:], param_name, **kwargs)
                # solve for S
                f_sensitivity_ode = lambda t, x : sensitivity_ode(t, x, J, Zj)
                sol = odeint(f_sensitivity_ode, S0, timepoints_ssm, tfirst = True)
                S = sol
                S = np.reshape(S, (len(timepoints_ssm), n))
                SSM[k,j,:] = S[k,:]
        if normalize:
            SSM = self.normalize_SSM(SSM, xs, self.M.get_parameter_values()) #Identifiablity was estimated using an normalized SSM
        return SSM

    def normalize_SSM(self, SSM, solutions, params_values):
        '''
        Returns normalized sensitivity coefficients. 
        Multiplies each sensitivity coefficient with the corresponding parameter p_j
        Divides the result by the corresponding state to obtain the normalized coefficient that is returned.
        '''
        n = np.shape(solutions)[1]
        SSM_normalized = np.zeros(np.shape(SSM))
        for j in range(len(params_values)):
            for i in range(n):
                SSM_normalized[:,j,i] = np.divide(SSM[:,j,i]*params_values[j], solutions[:,i]) 
        return SSM_normalized
