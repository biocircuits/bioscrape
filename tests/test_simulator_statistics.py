import numpy as np

from bioscrape.simulator import SSAResult


def test_empirical_distribution_uses_matching_integer_buffer_dtype():
    timepoints = np.array([0.0, 1.0, 2.0, 3.0])
    results = np.array([[0.0], [1.0], [1.0], [2.0]])
    result = SSAResult(timepoints, results)

    distribution = result.py_empirical_distribution(species=[0])

    assert np.allclose(distribution, np.array([0.25, 0.5, 0.25]))
