import bioscrape.random as bsr
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

def test_poisson_fit(mu, n_draws):
	bs_samples = np.array([bsr.py_poisson_rnd(mu) for _ in range(n_draws)])
	bs_bin_edges = [i-0.5 for i in range(bs_samples.max()+1)]
	bs_hist, bs_bin_edges = np.histogram(bs_samples, bins = bs_bin_edges)
	bs_bin_centers = [(bs_bin_edges[i]+bs_bin_edges[i+1])/2 for i in range(len(bs_bin_edges)-1)]

	poisson_support = range(bs_samples.max()+1)
	true_dist = scipy.stats.poisson.pmf(poisson_support, mu)

	print(f"len(bs_samples): {len(bs_samples)}")
	print(f"len(bs_hist): {len(bs_hist)}")
	print(f"len(bs_bin_edges): {len(bs_bin_edges)}")
	print(f"len(bs_bin_centers): {len(bs_bin_centers)}")
	print(f"len(true_dist): {len(true_dist)}")

	plt.plot(poisson_support, true_dist, color = 'black', lw = 2, label = "True Poisson")
	plt.plot(bs_bin_centers, bs_hist/bs_hist.sum(), color = 'red', lw = 1, label = "Sampled Poisson")
	plt.xlabel("Value")
	plt.ylabel("Frequency")
	plt.legend()
	plt.tight_layout()
	plt.show()

test_poisson_fit(1.5, 1000)
test_poisson_fit(5, 1000)
test_poisson_fit(80, 1000)