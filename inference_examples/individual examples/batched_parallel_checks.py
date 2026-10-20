"""Importable checks and timing instrumentation for the example notebook.

Worker functions live in this module so multiprocessing can import them
under Windows spawning. Instrumentation counts actual setup and likelihood
calls; it adds no artificial simulation or setup delays.
"""

import multiprocessing
import os
import time

import numpy as np
import pandas as pd

from bioscrape.inference_setup import InferenceSetup
from bioscrape.parallel_inference import (
    _initialize_likelihood_worker, _evaluate_likelihood)
from bioscrape.pid_interfaces import DeterministicInference
from bioscrape.simulator import py_simulate_model
from bioscrape.types import Model

_setup_calls = 0
_setup_seconds = 0.
_likelihood_calls = 0


class _CountedDeterministicInference(DeterministicInference):
    def get_log_likelihood(self, params):
        global _likelihood_calls
        _likelihood_calls += 1
        return super().get_log_likelihood(params)


class _CountedInferenceSetup(InferenceSetup):
    def setup_cost_function(self, **kwargs):
        super().setup_cost_function(**kwargs)
        self.pid_interface.__class__ = _CountedDeterministicInference
        self._likelihood_worker_config['interface_type'] = _CountedDeterministicInference

    def __setstate__(self, state):
        global _setup_calls, _setup_seconds
        start = time.perf_counter()
        super().__setstate__(state)
        _setup_seconds += time.perf_counter()-start
        _setup_calls += 1


def coupled_prior(params):
    """An example joint log prior, with rejection and a finite contribution."""
    if params['kprod'] > 18*params['kdeg']:
        return -np.inf
    return -.01*(params['kprod']-5.)**2


def make_setup(n_conditions=6, n_replicates=3, n_timepoints=51, counted=False):
    """Create a two-species CRN and synthetic multi-condition observations."""
    model = Model(
        species=['X', 'Y'],
        reactions=[([], ['X'], 'massaction', {'k': 'kprod'}),
                   (['X'], [], 'massaction', {'k': 'kdeg'}),
                   (['X'], ['X', 'Y'], 'massaction', {'k': 'ktl'}),
                   (['Y'], [], 'massaction', {'k': 'kloss'})],
        parameters={'kprod': 5., 'kdeg': .5, 'ktl': .6, 'kloss': .2},
        initial_condition_dict={'X': 0., 'Y': 0.})
    times = np.linspace(0, 5, n_timepoints)
    data, initial, conditions = [], [], []
    rng = np.random.default_rng(2026)
    for dose in np.linspace(.3, 1.2, n_conditions):
        state, parameters = {'X': float(dose), 'Y': 0.}, {'ktl': float(dose)}
        model.set_species(state)
        model.set_params(parameters)
        result = py_simulate_model(times, Model=model)
        for _ in range(n_replicates):
            data.append(pd.DataFrame({'time': times,
                'X': result.X.to_numpy()+rng.normal(0, .05, len(times)),
                'Y': result.Y.to_numpy()+rng.normal(0, .05, len(times))}))
            initial.append(state.copy())
            conditions.append(parameters.copy())
    model.set_species({'X': 0., 'Y': 0.})
    model.set_params({'ktl': .6})
    setup_class = _CountedInferenceSetup if counted else InferenceSetup
    return setup_class(
        Model=model, params_to_estimate=['kprod', 'kdeg'], exp_data=data,
        measurements=['X', 'Y'], time_column='time',
        initial_conditions=initial, parameter_conditions=conditions,
        prior={'kprod': ['log-uniform', .1, 20.],
               'kdeg': ['log-uniform', .01, 2.]},
        custom_joint_prior=coupled_prior)


def proposal_batches(n_batches=8, batch_size=64, outside_fraction=.5, seed=2026):
    """Make fixed joint tuples, including deliberately out-of-bounds proposals."""
    rng = np.random.default_rng(seed)
    batches = rng.uniform([1., .1], [12., 1.2], size=(n_batches, batch_size, 2))
    n_outside = int(round(batch_size*outside_fraction))
    for batch in batches:
        outside = rng.choice(batch_size, n_outside, replace=False)
        batch[outside, 0] = 25.
    return batches


def _reset_counters():
    global _setup_calls, _setup_seconds, _likelihood_calls
    _setup_calls, _setup_seconds, _likelihood_calls = 0, 0., 0


def _initialize_counted_pool(queue):
    _reset_counters()
    queue.put((os.getpid(), _setup_calls, _setup_seconds, _likelihood_calls))


def _initialize_counted_batched(config, queue):
    global _setup_calls, _setup_seconds
    _reset_counters()
    start = time.perf_counter()
    _initialize_likelihood_worker(config)
    _setup_seconds = time.perf_counter()-start
    _setup_calls = 1
    queue.put((os.getpid(), _setup_calls, _setup_seconds, _likelihood_calls))


def _report(value):
    return (value, os.getpid(), _setup_calls, _setup_seconds, _likelihood_calls)


def _pool_task(request):
    setup, params = request
    return _report(setup.cost_function(params))


def _batched_task(params):
    return _report(_evaluate_likelihood(params))


class _RecordingPool:
    def __init__(self, pool):
        self.pool = pool
        self.reports = []

    def map(self, function, params, chunksize):
        assert function is _evaluate_likelihood
        reports = self.pool.map(_batched_task, params, chunksize=chunksize)
        self.reports.extend(reports)
        return [row[0] for row in reports]


def benchmark_batch_methods(batches, n_processes=2, repeats=3):
    """Replay identical tuples and count worker setups for both methods.

    The pool replay uses the default multiprocessing chunk size and sends
    the inference object with each chunk, as the existing emcee pool route
    does. The batched replay calls the production batch evaluator and worker
    functions. Pool startup is measured separately from repeated batches.
    Worker setup seconds are summed CPU wall durations across processes.
    """
    rows = []
    outputs = {}
    context = multiprocessing.get_context('spawn')
    for repeat in range(repeats):
        methods = ('pool', 'batched') if repeat % 2 == 0 else ('batched', 'pool')
        for method in methods:
            setup = make_setup(counted=True)
            ready = context.Queue()
            start = time.perf_counter()
            if method == 'pool':
                pool = context.Pool(n_processes, _initialize_counted_pool, (ready,))
            else:
                pool = context.Pool(n_processes, _initialize_counted_batched,
                    (setup._get_likelihood_worker_config(), ready))
            try:
                worker_ready = [ready.get(timeout=60) for _ in range(n_processes)]
                startup_seconds = time.perf_counter()-start
                recorder = _RecordingPool(pool)
                all_results, reports = [], []
                evaluation_start = time.perf_counter()
                for batch in batches:
                    if method == 'pool':
                        returned = pool.map(_pool_task, [(setup, row) for row in batch])
                        reports.extend(returned)
                        all_results.append(np.array([row[0] for row in returned]))
                    else:
                        all_results.append(setup._cost_function_batch(batch, recorder))
                evaluation_seconds = time.perf_counter()-evaluation_start
                if method == 'batched':
                    reports = recorder.reports
                latest = {pid: (count, seconds, calls)
                          for pid, count, seconds, calls in worker_ready}
                for value, pid, count, seconds, likelihood_calls in reports:
                    latest[pid] = (count, seconds, likelihood_calls)
                rows.append(dict(method=method, repeat=repeat+1,
                    processes=n_processes, proposals=int(np.prod(batches.shape[:2])),
                    likelihood_calls=sum(row[2] for row in latest.values()),
                    worker_setups=sum(row[0] for row in latest.values()),
                    summed_worker_setup_seconds=sum(row[1] for row in latest.values()),
                    startup_seconds=startup_seconds,
                    batch_seconds=evaluation_seconds,
                    startup_plus_batch_seconds=startup_seconds+evaluation_seconds))
                outputs[(method, repeat)] = np.array(all_results)
            finally:
                pool.terminate()
                pool.join()
                ready.close()
                ready.join_thread()
        np.testing.assert_allclose(outputs[('pool', repeat)],
                                   outputs[('batched', repeat)], rtol=1e-12, atol=1e-12)
    result = pd.DataFrame(rows)
    for repeat, group in result.groupby('repeat'):
        assert group.likelihood_calls.nunique() == 1
        assert group.loc[group.method == 'batched', 'worker_setups'].item() == n_processes
        assert group.loc[group.method == 'pool', 'worker_setups'].item() > n_processes
    return result


