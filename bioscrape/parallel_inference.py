"""Persistent likelihood workers for batched parameter inference."""

_worker_pid_interface = None


def _initialize_likelihood_worker(config):
    """
    Build one likelihood per worker from prepared inference settings.

    Parameters
    ----------
    config : dict
        Model, prepared data, PID interface type, parameter names, default
        values, conditions and likelihood options from `InferenceSetup`.
        Contains no prior definitions or prior callbacks.
    """
    global _worker_pid_interface
    model = config['Model']
    model.set_params(config['default_parameters'])
    model.set_species(config['default_species'])
    _worker_pid_interface = config['interface_type'](
        config['params_to_estimate'], model, None,
        log_space_parameters=config['log_space_parameters'],
        debug=config['debug'])
    _worker_pid_interface.setup_likelihood_function(
        config['data'], config['timepoints'], config['measurements'],
        initial_conditions=config['initial_conditions'],
        parameter_conditions=config['parameter_conditions'],
        norm_order=config['norm_order'], **config['likelihood_kwargs'])


def _evaluate_likelihood(params):
    """
    Evaluate a tuple using this process's prepared likelihood.

    Parameters
    ----------
    params : array_like
        Candidate values in sampler coordinates, in parameter-name order.

    Returns
    -------
    float
        The log likelihood without a prior contribution.
    """
    return _worker_pid_interface.get_log_likelihood(params)
