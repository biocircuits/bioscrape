import warnings
from bioscrape.types import Model
from bioscrape.simulator import ModelCSimInterface, DeterministicSimulator
from scipy.integrate import odeint
import numpy as np
from typing import List, Union

def py_sensitivity_analysis(model: Model, timepoints: np.ndarray, 
                            normalize: bool, **kwargs) -> np.ndarray:
    """User interface function to perform sensitivity analysis 
    on a bioscrape model. The sensitivity coefficients are computed 
    where each coefficient s_ij = rate of change of x_i with parameter p_j
    for each time point in timepoints.

    Args:
        model (bioscrape.types.Model): A bioscrape Model object
        timepoints (numpy.ndarray): Array of time points.
        normalize (bool): when `True` the sensitivity coefficients returned are 
                          normalized by state values at each time 
                          (divides each coefficient by x_i/p_j).
                          when `False` the sensitivity coefficients are not normalized.

    Returns:
        numpy.ndarray: A numpy array of size:
                       len(timepoints) x len(parameters) x len(states)
    """
    dx = kwargs.get("dx", 0.01)
    precision = kwargs.get("precision", 10) 
    sens_obj = SensitivityAnalysis(model, dx=dx, precision=precision)
    ans_df = sens_obj.propagator.py_simulate(sens_obj.sim_interface, 
                                             timepoints).py_get_dataframe(sens_obj.M)
    solutions_array = np.array(ans_df.iloc[:,range(0,len(ans_df.T) - 1)])
    return sens_obj.compute_SSM(solutions_array, timepoints, normalize, **kwargs)

def py_get_jacobian(model: Model, state: Union[list, np.ndarray], **kwargs) -> np.ndarray:
    """User interfacce function to compute Jacobian (df/dx) of the model.

    Args:
        model (Model): Bioscrape Model
        state (Union[list, np.ndarray]): The state values (vector of length n) 
                                         at which to compute the Jacobian

    Returns:
        np.ndarray: A (n x n) Jacobian matrix, where n = len(state)
    """
    return SensitivityAnalysis(model).compute_J(state, **kwargs)

def py_get_sensitivity_to_parameter(model: Model, state: Union[list, np.ndarray], 
                                    param_name: str, **kwargs) -> np.ndarray:
    """User interface function to compute the sensitivity to parameter (df/dp)
    where p is the parameter and f is the model

    Args:
        model (Model): Bioscrape Model
        state (Union[list, np.ndarray]): The state values (vector of length n) 
                                         at which to compute df/dp
        param_name (str): The parameter name for which df/dp is computed

    Returns:
        np.ndarray: A np.ndarray of size (n x 1), where n is the length of state 
    """
    return SensitivityAnalysis(model).compute_Zj(state, param_name, **kwargs)

class SensitivityAnalysis(Model):
    def __init__(self, M, dx = 0.01, precision = 10):
        """
        Local Sensitivity Analysis for Bioscrape models.
        Arguments:
        * M: The Bioscrape Model object.
        * dx: Small parameter used in approximate computation methods. 
        * precision: the number of decimal places to round to
        """
        self.M = M
        sim = ModelCSimInterface(self.M)
        sim.py_prep_deterministic_simulation()
        self.sim_interface = sim
        self.propagator = DeterministicSimulator()
        self.num_equations = sim.py_get_num_species()
        self.dx = 0.01
        self.original_parameters = dict(M.get_parameter_dictionary())
        self.precision = precision
    
    def _evaluate_model(self, states, params = None, time = 0.0):
        """
        Numerically evaluates the model at a given value for the states and time (if time-varying model).
        """
        sim = self.sim_interface
        if params is not None:
            self.M.set_params(params)
        states = np.array(states, dtype = 'float64')
        derivative_array = np.zeros((self.num_equations), dtype = 'float64')
        sim.py_apply_repeated_rules(states, time, True)
        sim.py_calculate_deterministic_derivative(states, derivative_array, time)
        return derivative_array

    def compute_J(self, x, time = 0.0, **kwargs):
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
            f_0 = self._evaluate_model(state_input, time = time)[i]
            for j in range(n):
                h = self.dx
                if h == 0:
                    raise ValueError('Small parameter exactly equal to 0, cannot compute Jacobian')
                x = np.array(state_input)
                x[j] = x[j] + h
                f_h = self._evaluate_model(x, time = time)[i]
                x = np.array(state_input)
                x[j] = x[j] - h
                f_mh = self._evaluate_model(x, time = time)[i]
                if method == 'fourth_order_central_difference':
                    # Gets O(h^4) central difference on df_i/dvar_j
                    x = np.array(state_input)
                    x[j] = x[j] + 2*h
                    f_2h = self._evaluate_model(x, time = time)[i]
                    x = np.array(state_input)
                    x[j] = x[j] - 2*h
                    f_m2h = self._evaluate_model(x, time = time)[i]
                    J[i,j]= (-f_2h + 8*f_h - 8*f_mh + f_m2h)/(12*h)
                if method == 'central_difference':
                    J[i,j]= (f_h - f_mh)/(2*h) 
                if method == 'backward_difference':
                    J[i,j]= (f_0 - f_mh)/h
                if method == 'forward_difference':
                    J[i,j]= (f_h - f_0)/h
                # Error check
                if J[i, j] == np.inf:
                    warnings.warn('inf found while computing the Jacobian. Replacing by 1. Check model.')
                    J[i, j] = 1
                elif J[i, j] == np.nan:
                    warnings.warn('nan found while conputing the Jacobian. Replacing 0. Check model.')
                    J[i, j] = 0
        return np.round(J, decimals = self.precision)
        
    def compute_Zj(self, x, param_name, time = 0.0, **kwargs):
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
        params_dict = dict(self.original_parameters)
        array_f_0 = self._evaluate_model(x, params_dict, time = time)
        h = self.dx # Small parameter for this parameter
        # For each state
        for i in range(n):
            if h == 0:
                raise ValueError(f'Small parameter exactly equal to 0, cannot compute Zj for parameter {param_name}')
            f_0 = array_f_0[i]
            params_dict[param_name] = params_dict[param_name] + h
            self.M.set_params(params_dict)
            f_h = self._evaluate_model(x, params_dict, time = time)[i]
            # Reset
            params_dict = dict(self.original_parameters)
            self.M.set_params(params_dict)
            # Update
            params_dict[param_name] = params_dict[param_name] - h
            self.M.set_params(params_dict)
            f_mh = self._evaluate_model(x, params_dict, time = time)[i]
            # Reset
            params_dict = dict(self.original_parameters)
            self.M.set_params(params_dict)
            if method == 'fourth_order_central_difference':
                # Gets O(4) central difference on dfi/dpj
                params_dict[param_name] = params_dict[param_name] + 2*h
                self.M.set_params(params_dict)
                f_2h = self._evaluate_model(x, params_dict, time = time)[i]
                params_dict = dict(self.original_parameters)
                params_dict[param_name] = params_dict[param_name] - 2*h
                self.M.set_params(params_dict)
                f_m2h = self._evaluate_model(x, params_dict, time = time)[i]
                params_dict = dict(self.original_parameters)
                self.M.set_params(params_dict)
                #Store approx. dfi/dp[param_name] into Z
                Z[i]= (-f_2h + 8*f_h - 8*f_mh + f_m2h)/(12*h)
            if method == 'central_difference':
                Z[i]= (f_h - f_mh)/(2*h) 
            if method == 'backward_difference':
                Z[i]= (f_0 - f_mh)/h
            if method == 'forward_difference':
                Z[i]= (f_h - f_0)/h
            # Error check
            if Z[i] == np.inf:
                warnings.warn('inf found while compute Zj, replacing by 1. Check model.')
                Z[i] = 1
            elif Z[i] == np.nan:
                warnings.warn('nan found while compute Zj, replacing by 0. Check model.')
                Z[i] = 0
        return np.round(Z, decimals = self.precision)

    def compute_SSM(self, solutions, timepoints, normalize = False, params = None, **kwargs):
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
        * params: (list of parameters, default is None): The parameters to compute sensitivty to. When None defaults to all model parameters
        * kwargs: Other kwargs passed to `compute_J` and `compute_Z` functions.
        """
        def sensitivity_ode(t, x, J, Z):
            # ODE to solve for sensitivity coefficient S
            dsdt = J@x + Z
            return dsdt

        if params is None:
            all_params = list(self.original_parameters.keys())
        else:
            all_params = params

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
            J = self.compute_J(xs[k,:], time = timepoints[k], **kwargs)
            #Solve for S = dx/dp for all x and all P (or theta, the parameters) at time point k
            for j in range(len(all_params)):
                param_name = all_params[j]
                # get the pmatrix
                Zj = self.compute_Zj(xs[k,:], param_name, time = timepoints[k], **kwargs)
                # solve for S
                f_sensitivity_ode = lambda t, x : sensitivity_ode(t, x, J, Zj)
                sol = odeint(f_sensitivity_ode, S0, timepoints_ssm, tfirst = True)
                S = sol
                S = np.reshape(S, (len(timepoints_ssm), n))
                SSM[k,j,:] = S[k,:]
        if normalize:
            param_dict = self.M.get_parameter_dictionary()
            param_vals = np.array([param_dict[p] for p in all_params])
            SSM = self.normalize_SSM(SSM, xs, param_vals) #Identifiablity was estimated using an normalized SSM
        return np.round(SSM, decimals = self.precision)

    def normalize_SSM(self, SSM, solutions, params_values):
        """
        Returns normalized sensitivity coefficients. 
        Multiplies each sensitivity coefficient with the corresponding parameter p_j
        Divides the result by the corresponding state to obtain the normalized coefficient that is returned.
        """
        n = np.shape(solutions)[1]
        SSM_normalized = np.zeros(np.shape(SSM))
        for j in range(len(params_values)):
            for i in range(n):
                SSM_normalized[:,j,i] = np.divide(SSM[:,j,i]*params_values[j], solutions[:,i]) 
        return SSM_normalized
