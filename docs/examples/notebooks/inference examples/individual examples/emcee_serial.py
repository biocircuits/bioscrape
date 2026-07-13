# Run: python emcee_series.py to benchmark emcee (serially) on CPU
# Import packages
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import emcee
import time

import logging
logging.basicConfig(filename='emcee_serial.log', level=logging.INFO)
logging.info("Instance run at {}".format(time.asctime()))
# Define log-prior, log-likelihood, and log-probability functions
def log_prior(theta):
    m, b = theta
    # Try priors
    # if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < log_f < 1.0:
    if -500.0 < m < 500 and -1100.0 < b < 1100.0:
        return 0.0
    return -np.inf

def log_likelihood(theta, x, y):
    m, b = theta
    model = m * x + b
    cost = - np.linalg.norm(y - model)
    return cost

def log_probability(theta, x, y):
    # make the log-prob slow to benchmark serial vs parallel
    # so that overheads are compensated.
    t = time.time() + np.random.uniform(0.005, 0.008)
    while True:
        if time.time() >= t:
            break
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y)

np.random.seed(123)

# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534 # noise model
logging.info("True parameters are m = {}, b = {}, f = {}".format(m_true, b_true, f_true))
# Generate some synthetic data from the model.
N = 500
x = np.sort(10 * np.random.rand(N))
y = m_true * x + b_true
y += np.abs(f_true * y) * np.random.randn(N)

x0 = np.linspace(0, 10, N)

# Inference params
nwalkers = 32
nsteps = 1000
pos = 1*np.zeros(2) + 1e-4 * np.random.randn(nwalkers, 2)
nwalkers, ndim = pos.shape

# Usual emcee call (Not parallel)
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,
                                args=(x, y))
start_time = time.time()
sampler.run_mcmc(pos, nsteps, progress=True)
elapsed = time.time() - start_time
print("Time without parallelization: {:.2f} seconds".format(elapsed))
logging.info("Time without parallelization: {:.2f} seconds".format(elapsed))
fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")
plt.savefig("emcee_serial.png")


# Get autocorrelation time
# tau = sampler.get_autocorr_time()
# print("Autocorrelation time is", tau)


flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

import corner

fig = corner.corner(
    flat_samples, labels=labels, truths=[m_true, b_true]
)
plt.savefig("corner_plot_parallel.png")

mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
print("The 16th 50th and 84th percentiles of the estimated params are", mcmc)
print("Results logged to emcee_serial.log")