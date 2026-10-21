"""Contracts for scalar and batched custom joint-prior evaluation."""

import numpy as np
import pytest

from bioscrape.pid_interfaces import (
    PIDInterface, DeterministicInference, StochasticInference)
from bioscrape.types import Model


@pytest.fixture
def prior_interface():
    model = Model(species=['X'], parameters={'a': .3, 'b': .5})
    return PIDInterface(['a', 'b'], model, {
        'a': ['uniform', .01, 2.], 'b': ['uniform', .01, 2.]})


@pytest.mark.parametrize('spec', [
    ['uniform', .01, 2.], ['gaussian', .4, .1, 'positive'],
    ['exponential', 2.], ['gamma', 2., 4.], ['beta', 2., 4.],
    ['log-uniform', .01, 2.], ['log-gaussian', -1., 1.],
    ['custom', lambda name, value: -.5*value**2],
])
def test_batch_joint_matches_scalar_for_all_parameter_priors(prior_interface, spec):
    interface = prior_interface
    interface.prior = {name: spec for name in interface.params_to_estimate}
    calls = []

    def scalar(params):
        return -.7 if params['a'] < .7 else -np.inf

    def batch(params):
        calls.append(params)
        return np.where(params['a'] < .7, -.7, -np.inf)

    interface.custom_joint_prior = scalar
    interface.custom_joint_prior_batch = batch
    points = np.array([[.3, .5], [.8, .4], [.2, .6], [.3, .5]])
    expected = np.array([interface.get_log_prior(row) for row in points])
    np.testing.assert_allclose(interface.get_log_prior_batch(points), expected)
    assert len(calls) == 1
    np.testing.assert_array_equal(calls[0]['a'], points[:, 0])


def test_batch_callback_replaces_scalar_callback_in_batch(prior_interface):
    def unused_scalar(params):
        raise AssertionError('The scalar joint prior must not run in a batch.')

    prior_interface.custom_joint_prior = unused_scalar
    prior_interface.custom_joint_prior_batch = lambda p: np.full(len(p['a']), -.7)
    points = np.array([[.3, .5], [.4, .6]])
    expected = 2*np.log(1/1.99) - .7
    np.testing.assert_allclose(prior_interface.get_log_prior_batch(points), expected)


def test_parameter_rejection_and_empty_batch_skip_joint_callback(prior_interface):
    received = []

    def callback(params):
        received.append(params)
        return np.zeros(len(params['a']))

    prior_interface.custom_joint_prior_batch = callback
    points = np.array([[4., .5], [.3, .5], [.4, -1.], [.2, .6]])
    result = prior_interface.get_log_prior_batch(points)
    np.testing.assert_array_equal(received[0]['a'], [.3, .2])
    np.testing.assert_array_equal(received[0]['b'], [.5, .6])
    assert np.isneginf(result[[0, 2]]).all()
    assert np.isfinite(result[[1, 3]]).all()
    prior_interface.get_log_prior_batch(points[[0, 2]])
    assert prior_interface.get_log_prior_batch(np.empty((0, 2))).shape == (0,)
    assert prior_interface.get_log_prior([4., .5]) == -np.inf
    assert len(received) == 1


@pytest.mark.parametrize('interface_type, attribute', [
    (DeterministicInference, 'LL_det'), (StochasticInference, 'LL_stoch')])
def test_batch_only_callback_is_used_in_scalar_inference(
        prior_interface, interface_type, attribute):
    received = []

    def callback(params):
        received.append(params)
        return np.where(params['a'] < .7, -.4, -np.inf)

    interface = interface_type(
        prior_interface.params_to_estimate, prior_interface.M,
        prior_interface.prior, custom_joint_prior_batch=callback)

    class Likelihood:
        calls = 0

        def set_init_params(self, params):
            pass

        def py_log_likelihood(self):
            self.calls += 1
            return 3.

    likelihood = Likelihood()
    setattr(interface, attribute, likelihood)
    np.testing.assert_allclose(
        interface.get_likelihood_function([.3, .5]), 2*np.log(1/1.99)-.4+3.)
    assert interface.get_likelihood_function([.8, .5]) == -np.inf
    assert likelihood.calls == 1
    assert all(params['a'].shape == (1,) for params in received)


def test_log_space_batch_callback_receives_model_values(prior_interface):
    received = []
    prior_interface.log_space_parameters = True
    prior_interface.custom_joint_prior_batch = (
        lambda p: received.append(p) or -.25*p['a'])
    points = np.array([[.3, .5], [.4, .6]])
    expected = np.array([prior_interface.get_log_prior(row)
                         for row in np.log(points)])
    actual = prior_interface.get_log_prior_batch(np.log(points))
    np.testing.assert_allclose(actual, expected)
    np.testing.assert_allclose(received[-1]['a'], points[:, 0])
    np.testing.assert_allclose(received[-1]['b'], points[:, 1])


def test_nonfinite_joint_values_reject_individual_rows(prior_interface):
    prior_interface.custom_joint_prior_batch = lambda p: np.array([
        0., -np.inf, np.inf, np.nan, -.7])
    result = prior_interface.get_log_prior_batch(np.tile([.3, .5], (5, 1)))
    assert np.isneginf(result[1:4]).all()
    np.testing.assert_allclose(result[[0, 4]], 2*np.log(1/1.99)+np.array([0., -.7]))


@pytest.mark.parametrize('output', [0., [0.], [[0.], [0.]]])
def test_wrong_joint_output_shape_raises(prior_interface, output):
    prior_interface.custom_joint_prior_batch = lambda p: output
    with pytest.raises(ValueError, match='one value per proposal'):
        prior_interface.get_log_prior_batch([[.3, .5], [.4, .6]])


@pytest.mark.parametrize('points', [[.3, .5], [[.3]], np.zeros((1, 2, 1))])
def test_wrong_input_shape_raises(prior_interface, points):
    with pytest.raises(ValueError, match='proposals, parameters'):
        prior_interface.get_log_prior_batch(points)


def test_batch_callback_exception_propagates(prior_interface):
    def callback(params):
        raise RuntimeError('custom prior failure')

    prior_interface.custom_joint_prior_batch = callback
    with pytest.raises(RuntimeError, match='custom prior failure'):
        prior_interface.get_log_prior_batch([[.3, .5]])


def test_existing_check_prior_override_is_preserved_without_batch_hook(prior_interface):
    class OverriddenPrior(PIDInterface):
        def check_prior(self, params_dict):
            return -params_dict['a']**2

    interface = OverriddenPrior(['a', 'b'], prior_interface.M, prior_interface.prior)
    np.testing.assert_allclose(
        interface.get_log_prior_batch([[.3, .5], [.4, .6]]), [-.09, -.16])
