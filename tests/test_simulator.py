import pytest
import test_utils
import numpy as np

from bioscrape.types import Model
from bioscrape.simulator import py_simulate_model

#This test confirms that SSAResult.py_empirical_distribution runs
#  without crashing and returns a sensible distribution. The cdef
#  helper it calls (empirical_distribution) declared its index_ar/
#  max_counts arrays with an np.int32_t buffer type but allocated them
#  with dtype=np.int_ (platform-dependent, 64-bit "long" on most
#  modern 64-bit platforms), so NumPy's buffer protocol rejected the
#  mismatch with "ValueError: Buffer dtype mismatch, expected
#  'int32_t' but got 'long'" on essentially any 64-bit platform. Fixed
#  by allocating with dtype=np.int32 to match the declared buffer
#  type.
def test_empirical_distribution():
    test_utils.set_seed(1234)

    M = Model()
    M.create_reaction(reactants=[], products=['X'],
                      propensity_type='massaction',
                      propensity_param_dict={'k': 'k_prod'})
    M.create_reaction(reactants=['X'], products=[],
                      propensity_type='massaction',
                      propensity_param_dict={'k': 'k_deg'})
    M.set_parameter('k_prod', 10.0)
    M.set_parameter('k_deg', 0.5)
    M.set_species({'X': 0})

    timepoints = np.linspace(0, 2000, 20000)
    result = py_simulate_model(timepoints, Model=M, stochastic=True,
                               return_dataframe=False)
    dist = result.py_empirical_distribution(start_time=200, species=['X'],
                                            Model=M)

    assert np.isclose(dist.sum(), 1.0)

    # Steady state for this birth-death model is Poisson(k_prod/k_deg = 20),
    # so both the mean and the variance should be close to 20.
    xs = np.arange(len(dist))
    mean = (xs * dist).sum()
    var = ((xs - mean)**2 * dist).sum()
    assert np.isclose(mean, 20, atol=2)
    assert np.isclose(var, 20, atol=4)
