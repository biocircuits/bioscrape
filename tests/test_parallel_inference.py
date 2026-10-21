"""Numerical and multiprocessing checks for batched emcee inference."""

import multiprocessing
import os
import pickle
import time

import numpy as np
import pandas as pd
import pytest

from bioscrape.inference import py_inference
from bioscrape.inference_setup import InferenceSetup
from bioscrape.parallel_inference import (
    _initialize_likelihood_worker, _evaluate_likelihood)
from bioscrape.pid_interfaces import DeterministicInference
from bioscrape.simulator import py_simulate_model
from bioscrape.types import Model


def _worker_identity(params):
    import bioscrape.parallel_inference as worker
    time.sleep(.01)
    return (os.getpid(), id(worker._worker_pid_interface),
            _evaluate_likelihood(params))


def _scalar_coupled_prior(params):
    return -.25*params['kprod'] if params['kprod'] < 2.5 else -np.inf


def _batch_coupled_prior(params):
    return np.where(params['kprod'] < 2.5, -.25*params['kprod'], -np.inf)


class FailingLikelihood(DeterministicInference):
    def get_log_likelihood(self, params):
        raise ValueError('worker likelihood failure')


@pytest.fixture
def inference_setup():
    model = Model(
        species=['X', 'Y'],
        reactions=[([], ['X'], 'massaction', {'k': 'kprod'}),
                   (['X'], [], 'massaction', {'k': 'kdeg'}),
                   (['X'], ['X', 'Y'], 'massaction', {'k': 'ktl'}),
                   (['Y'], [], 'massaction', {'k': 'kloss'})],
        parameters={'kprod': 2., 'kdeg': .3, 'ktl': .6, 'kloss': .2},
        initial_condition_dict={'X': 0., 'Y': 0.})
    times = np.linspace(0, 2, 9)
    initial = [{'X': 0., 'Y': 0.}, {'X': 1., 'Y': 2.}]
    conditions = [{'ktl': .4}, {'ktl': .8}]
    data = []
    for state, parameters in zip(initial, conditions):
        model.set_species(state)
        model.set_params(parameters)
        frame = py_simulate_model(times, Model=model)
        data.append(pd.DataFrame({'time': times, 'X': frame.X, 'Y': frame.Y}))
    model.set_species({'X': 0., 'Y': 0.})
    model.set_params({'ktl': .6})
    return InferenceSetup(
        Model=model, params_to_estimate=['kprod', 'kdeg'],
        exp_data=data, measurements=['X', 'Y'], time_column='time',
        initial_conditions=initial, parameter_conditions=conditions,
        prior={'kprod': ['uniform', .1, 5.], 'kdeg': ['uniform', .01, 2.]},
        hmax=.1)


@pytest.fixture
def likelihood_pool(inference_setup):
    with multiprocessing.get_context('spawn').Pool(
            processes=2, initializer=_initialize_likelihood_worker,
            initargs=(inference_setup._get_likelihood_worker_config(),)) as pool:
        yield pool


def test_batched_log_posterior_matches_scalar(inference_setup, likelihood_pool):
    params = np.array([[2., .3], [.2, 1.], [8., .3], [1., .1], [2., .3]])
    scalar = np.array([inference_setup.cost_function(row) for row in params])
    batch = inference_setup._cost_function_batch(params, likelihood_pool)
    np.testing.assert_allclose(batch, scalar, rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(inference_setup.cost_params[-len(params):], params)
    np.testing.assert_array_equal(inference_setup.cost_progress[-len(params):], batch)


def test_all_prior_types_in_batch(inference_setup, likelihood_pool):
    params = np.array([[.3, .5], [.5, .7]])
    for spec in (
            ['uniform', .01, 2.], ['gaussian', .4, .1, 'positive'],
            ['exponential', 2.], ['gamma', 2., 4.], ['beta', 2., 4.],
            ['log-uniform', .01, 2.], ['log-gaussian', -1., 1.],
            ['custom', lambda name, value: -.5*value**2]):
        inference_setup.pid_interface.prior = {
            name: spec for name in inference_setup.params_to_estimate}
        scalar = [inference_setup.cost_function(row) for row in params]
        batch = inference_setup._cost_function_batch(params, likelihood_pool)
        np.testing.assert_allclose(batch, scalar, rtol=1e-12, atol=1e-12)


def test_prior_callbacks_run_in_parent(inference_setup):
    parent_pid = os.getpid()
    calls = []

    def parameter_prior(name, value):
        assert os.getpid() == parent_pid
        return -.25*value

    def joint_prior(params):
        assert os.getpid() == parent_pid
        calls.append(params.copy())
        return -.75 if params['kprod'] < 3. else -np.inf

    inference_setup.prior['kprod'] = ['custom', parameter_prior]
    inference_setup.custom_joint_prior = joint_prior
    inference_setup.setup_cost_function()
    config = inference_setup._get_likelihood_worker_config()
    pickle.dumps(config)  # Local callbacks must not be serialized to workers.
    params = np.array([[2., .3], [4., .3]])
    with multiprocessing.get_context('spawn').Pool(
            1, _initialize_likelihood_worker, (config,)) as pool:
        batch = inference_setup._cost_function_batch(params, pool)
    assert len(calls) == 2
    expected = (-.25*2 - .75 + np.log(1/1.99)
                + inference_setup.pid_interface.get_log_likelihood(params[0]))
    np.testing.assert_allclose(batch[0], expected)
    assert batch[1] == -np.inf


def test_rejected_parameters_never_reach_workers(inference_setup):
    class RecordingPool:
        def map(self, function, params, chunksize):
            self.params = params.copy()
            assert chunksize == 1
            return [10., 20.]

    pool = RecordingPool()
    params = np.array([[2., .3], [8., .3], [1., .4]])
    batch = inference_setup._cost_function_batch(params, pool)
    np.testing.assert_array_equal(pool.params, params[[0, 2]])
    prior = np.array([inference_setup.pid_interface.get_log_prior(p)
                      for p in params[[0, 2]]])
    np.testing.assert_allclose(batch[[0, 2]], prior + [10., 20.])
    assert batch[1] == -np.inf


def test_completely_rejected_batch(inference_setup):
    class UnusedPool:
        def map(self, *args, **kwargs):
            raise AssertionError('A fully rejected batch must not dispatch tasks.')

    batch = inference_setup._cost_function_batch([[9., .3], [2., -1.]], UnusedPool())
    np.testing.assert_array_equal(batch, [-np.inf, -np.inf])


def test_worker_likelihood_initialized_once(inference_setup, likelihood_pool):
    rows = likelihood_pool.map(_worker_identity, [[2., .3]]*16, chunksize=1)
    second = likelihood_pool.map(_worker_identity, [[2., .3]]*16, chunksize=1)
    identities = {}
    for pid, identity, value in rows + second:
        identities.setdefault(pid, set()).add(identity)
        assert np.isfinite(value)
    assert len(identities) == 2
    assert all(len(values) == 1 for values in identities.values())


def test_worker_state_reset(inference_setup, likelihood_pool):
    first = likelihood_pool.map(_evaluate_likelihood, [[2., .3]]*4, chunksize=1)
    likelihood_pool.map(_evaluate_likelihood, [[.2, 1.5], [4.5, .02]]*4, chunksize=1)
    repeated = likelihood_pool.map(_evaluate_likelihood, [[2., .3]]*4, chunksize=1)
    np.testing.assert_allclose(repeated, first, rtol=1e-12, atol=1e-12)


def test_worker_preserves_likelihood_options(inference_setup):
    inference_setup.setup_cost_function(hmax=.07, rtol=1e-8, atol=1e-9,
                                       log_space_parameters=True)
    config = inference_setup._get_likelihood_worker_config()
    assert config['likelihood_kwargs']['hmax'] == .07
    assert config['likelihood_kwargs']['rtol'] == 1e-8
    params = np.log([[2., .3], [1., .4]])
    scalar = [inference_setup.cost_function(row) for row in params]
    with multiprocessing.get_context('spawn').Pool(
            1, _initialize_likelihood_worker, (config,)) as pool:
        batch = inference_setup._cost_function_batch(params, pool)
    np.testing.assert_allclose(batch, scalar, rtol=1e-12, atol=1e-12)


def test_worker_exceptions_propagate(inference_setup):
    config = inference_setup._get_likelihood_worker_config()
    config['interface_type'] = FailingLikelihood
    with multiprocessing.get_context('spawn').Pool(
            1, _initialize_likelihood_worker, (config,)) as pool:
        with pytest.raises(ValueError, match='worker likelihood failure'):
            inference_setup._cost_function_batch([[2., .3]], pool)


@pytest.mark.parametrize('joint_mode', ['none', 'scalar_and_batch', 'batch_only'])
def test_batched_emcee_matches_serial(inference_setup, tmp_path, joint_mode):
    start = np.array([2., .3]) * (1 + .05*np.random.default_rng(2026).normal(size=(8, 2)))
    chains, probabilities = [], []
    for mode in ('serial', 'pool', 'batched'):
        np.random.seed(2026)
        callbacks = {}
        if joint_mode != 'none':
            callbacks['custom_joint_prior_batch'] = _batch_coupled_prior
        if joint_mode == 'scalar_and_batch':
            callbacks['custom_joint_prior'] = _scalar_coupled_prior
        sampler, setup = py_inference(
            Model=inference_setup.M, exp_data=inference_setup.exp_data,
            params_to_estimate=inference_setup.params_to_estimate,
            prior=inference_setup.prior, measurements=['X', 'Y'],
            initial_conditions=inference_setup.initial_conditions,
            parameter_conditions=inference_setup.parameter_conditions,
            time_column='time', nwalkers=8, nsteps=8, init_seed=start,
            parallel=mode != 'serial', n_processes=2,
            parallel_method='batched' if mode == 'batched' else 'pool',
            plot_show=False, printout=False, progess=False,
            filename_csv=str(tmp_path / (mode+'.csv')),
            filename_txt=str(tmp_path / (mode+'.txt')), **callbacks)
        chains.append(sampler.get_chain())
        probabilities.append(sampler.get_log_prob())
        if mode == 'batched':
            assert len(setup.cost_progress) == 8 + 8*8
    for chain, prob in zip(chains[1:], probabilities[1:]):
        np.testing.assert_allclose(chain, chains[0], rtol=0, atol=0)
        np.testing.assert_allclose(prob, probabilities[0], rtol=1e-12, atol=1e-12)


def test_batch_joint_callback_stays_in_parent(inference_setup):
    parent_pid = os.getpid()
    calls = []

    def batch_prior(params):
        assert os.getpid() == parent_pid
        calls.append(params)
        return np.where(params['kprod'] < 3., -.75, -np.inf)

    inference_setup.custom_joint_prior_batch = batch_prior
    inference_setup.setup_cost_function()
    config = inference_setup._get_likelihood_worker_config()
    assert 'custom_joint_prior_batch' not in config['likelihood_kwargs']
    pickle.dumps(config)
    params = np.array([[2., .3], [8., .3], [4., .3], [1., .4]])
    with multiprocessing.get_context('spawn').Pool(
            1, _initialize_likelihood_worker, (config,)) as pool:
        actual = inference_setup._cost_function_batch(params, pool)
    assert len(calls) == 1
    np.testing.assert_array_equal(calls[0]['kprod'], [2., 4., 1.])
    assert np.isneginf(actual[[1, 2]]).all()
    expected = [inference_setup.cost_function(row) for row in params]
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


def test_stochastic_batched_inference(inference_setup, tmp_path):
    inference_setup.sim_type = 'stochastic'
    inference_setup.N_simulations = 2
    inference_setup.parallel = True
    inference_setup.nwalkers = 8
    inference_setup.nsteps = 3
    inference_setup.init_seed = np.array([2., .3]) * (
        1 + .05*np.random.default_rng(2026).normal(size=(8, 2)))
    sampler = inference_setup.run_mcmc(
        parallel_method='batched', n_processes=2, N_simulations=2, norm_order=1,
        timepoints=inference_setup.timepoints, printout=False, progess=False,
        filename_csv=str(tmp_path / 'stochastic.csv'),
        filename_txt=str(tmp_path / 'stochastic.txt'))
    assert sampler.get_chain().shape == (3, 8, 2)
    assert np.isfinite(sampler.get_log_prob()).all()
    config = inference_setup._get_likelihood_worker_config()
    assert config['likelihood_kwargs']['N_simulations'] == 2
    assert config['norm_order'] == 1


@pytest.mark.parametrize('parallel, method', [(True, 'unknown'), (False, 'batched')])
def test_invalid_parallel_method(inference_setup, parallel, method):
    inference_setup.parallel = parallel
    with pytest.raises(ValueError, match='parallel_method'):
        inference_setup.run_emcee(parallel_method=method)
