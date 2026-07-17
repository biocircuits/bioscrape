"""
Local sensitivity analysis of bioscrape models.

Provides finite-difference sensitivity analysis of a deterministically
simulated `~bioscrape.types.Model` with respect to its parameters, via
`SensitivityAnalysis` and the module-level convenience functions
`py_sensitivity_analysis`, `py_get_jacobian`, and
`py_get_sensitivity_to_parameter`.
"""
import warnings
from bioscrape.types import Model
from bioscrape.simulator import ModelCSimInterface, DeterministicSimulator
from scipy.integrate import odeint
import numpy as np
from typing import List, Union

def py_sensitivity_analysis(model: Model, timepoints: np.ndarray,
                            normalize: bool, **kwargs) -> np.ndarray:
    r"""
    Compute the sensitivity of a model's trajectories to its parameters.

    Simulates `model` deterministically over `timepoints` and returns the
    sensitivity coefficients s_ij = d(x_i)/d(p_j) for each
    state x_i, parameter p_j, and time point.

    Parameters
    ----------
    model : Model
        The bioscrape model to analyze.
    timepoints : numpy.ndarray
        Array of time points at which to compute the sensitivity coefficients.
    normalize : bool
        If True, each sensitivity coefficient is normalized by x_i/p_j at
        that time point. If False, the coefficients are not normalized.
    dx : float, default 0.01
        Finite-difference step size used when approximating derivatives.
    precision : int, default 10
        Number of decimal places to round returned sensitivity values to.
    method : str, optional
        The finite-difference method to use; see
        `SensitivityAnalysis.compute_J`.

    Returns
    -------
    numpy.ndarray
        Array of shape `(len(timepoints), len(parameters), len(states))`
        containing the sensitivity coefficients.
    """
    dx = kwargs.get("dx", 0.01)
    precision = kwargs.get("precision", 10) 
    sens_obj = SensitivityAnalysis(model, dx=dx, precision=precision)
    ans_df = sens_obj.propagator.py_simulate(sens_obj.sim_interface, 
                                             timepoints).py_get_dataframe(sens_obj.M)
    solutions_array = np.array(ans_df.iloc[:,range(0,len(ans_df.T) - 1)])
    return sens_obj.compute_SSM(solutions_array, timepoints, normalize, **kwargs)

def py_get_jacobian(model: Model, state: Union[list, np.ndarray], **kwargs) -> np.ndarray:
    r"""
    Compute the Jacobian J = d(f)/d(x) of a model.

    Parameters
    ----------
    model : Model
        The bioscrape model to analyze.
    state : list or numpy.ndarray
        The state values (vector of length `n`) at which to compute the
        Jacobian.
    method : str, optional
        The finite-difference method to use; see
        `SensitivityAnalysis.compute_J`.
    time : float, optional
        The time at which to evaluate the (possibly time-varying)
        model (default 0.0).

    Returns
    -------
    numpy.ndarray
        The `n` x `n` Jacobian matrix, where `n = len(state)`.
    """
    return SensitivityAnalysis(model).compute_J(state, **kwargs)

def py_get_sensitivity_to_parameter(model: Model, state: Union[list, np.ndarray],
                                    param_name: str, **kwargs) -> np.ndarray:
    r"""
    Compute a model's sensitivity Z_j = d(f)/d(p_j) to p_j.

    Parameters
    ----------
    model : Model
        The bioscrape model to analyze.
    state : list or numpy.ndarray
        The state values (vector of length `n`) at which to compute
        d(f)/d(p_j).
    param_name : str
        The name of the parameter p_j to differentiate with respect to.
    method : str, optional
        The finite-difference method to use; see
        `SensitivityAnalysis.compute_Zj`.
    time : float, optional
        The time at which to evaluate the (possibly time-varying)
        model (default 0.0).

    Returns
    -------
    numpy.ndarray
        Vector of length `n` giving d(f_i)/d(p_j) for each state i.
    """
    return SensitivityAnalysis(model).compute_Zj(state, param_name, **kwargs)

class SensitivityAnalysis(Model):
    r"""
    Local sensitivity analysis for bioscrape models.

    Computes the sensitivity of a model's deterministic trajectories to its
    parameters, using finite-difference approximations of the Jacobian
    J = d(f)/d(x) and the parameter-sensitivity vector
    Z_j = d(f)/d(p_j), where f is the model's right-hand side.
    See `py_sensitivity_analysis` for the higher-level, user-facing interface
    built on this class.

    Parameters
    ----------
    M : Model
        The bioscrape model to analyze.
    dx : float, optional
        Finite-difference step size used when approximating derivatives
        (default 0.01).
    precision : int, optional
        Number of decimal places to round returned sensitivity values to
        (default 10).
    """
    def __init__(self, M, dx = 0.01, precision = 10):
        """See class docstring."""
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
        r"""
        Compute the Jacobian J = d(f)/d(x) at a point `x`.

        Uses a finite-difference approximation of the model's right-hand side
        f; see `method` below for the available approximation orders.

        Parameters
        ----------
        x : array_like
            The state (vector of length `n`) at which to evaluate the
            Jacobian.
        time : float, optional
            The time at which to evaluate the (possibly time-varying) model
            (default 0.0).
        method : str, optional
            The finite-difference method to use:
            'fourth_order_central_difference' (default), 'central_difference',
            'backward_difference', or 'forward_difference'.

        Returns
        -------
        numpy.ndarray
            The `n` x `n` Jacobian matrix.
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
        r"""
        Compute Z_j = d(f)/d(p_j) at a point `x`.

        Uses the same finite-difference methods as `compute_J`, where p_j
        is the parameter named `param_name`.

        Parameters
        ----------
        x : array_like
            The state (vector of length `n`) at which to evaluate the
            sensitivity.
        param_name : str
            The name of the parameter p_j to differentiate with respect to.
        time : float, optional
            The time at which to evaluate the (possibly time-varying) model
            (default 0.0).
        method : str, optional
            The finite-difference method to use; see `compute_J`.

        Returns
        -------
        numpy.ndarray
            Vector of length `n` giving d(f_i)/d(p_j) for each state i.
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
        r"""
        Compute the sensitivity coefficients for each parameter.

        Computes S_j = d(x)/d(p_j) for each parameter p_j, at
        each time point, by integrating the sensitivity ODE dS/dt = J*S + Z_j
        forward in time using the Jacobian and parameter-sensitivity vectors
        from `compute_J` and `compute_Zj`.

        Parameters
        ----------
        solutions : array_like
            The deterministic trajectory of the model (e.g. as returned by
            `~bioscrape.simulator.py_simulate_model`) evaluated at
            `timepoints`, of shape `(len(timepoints), n)`.
        timepoints : numpy.ndarray
            The time points at which `solutions` was computed and at which
            sensitivity coefficients are returned.
        normalize : bool, optional
            If True, normalize the returned coefficients by the nominal
            parameter and state values (default False).
        params : list of str, optional
            The parameters to compute sensitivities for. Defaults to all of
            the model's parameters.
        method : str, optional
            The finite-difference method to use; see `compute_J`.

        Returns
        -------
        numpy.ndarray
            Array of shape `(len(timepoints), len(params), n)` containing the
            sensitivity coefficients S_j for each parameter and time point.
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
        Normalize sensitivity coefficients by parameter and state values.

        Multiplies each sensitivity coefficient by its corresponding parameter
        value and divides by the corresponding state value.

        Parameters
        ----------
        SSM : numpy.ndarray
            Sensitivity coefficients as returned by `compute_SSM`, of shape
            `(len(timepoints), len(params), n)`.
        solutions : array_like
            The state trajectory used to compute `SSM`, of shape
            `(len(timepoints), n)`.
        params_values : array_like
            The nominal value of each parameter in `SSM`, in the same order.

        Returns
        -------
        numpy.ndarray
            The normalized sensitivity coefficients, same shape as `SSM`.
        """
        n = np.shape(solutions)[1]
        SSM_normalized = np.zeros(np.shape(SSM))
        for j in range(len(params_values)):
            for i in range(n):
                SSM_normalized[:,j,i] = np.divide(SSM[:,j,i]*params_values[j], solutions[:,i]) 
        return SSM_normalized
