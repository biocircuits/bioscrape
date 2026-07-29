.. currentmodule:: bioscrape.analysis

Sensitivity Analysis
=====================

Sensitivity analysis answers a simple question: if a parameter's value
changed slightly, how much would the model's trajectories change?
This is useful on its own, for finding which parameters most strongly
drive a species' behavior over time, but it also has a more practical
role as a precursor to :doc:`parameter inference <parameter_inference>`.
Not every parameter that appears in a model can actually be estimated
from a given measurement -- if a species' trajectory barely changes as
a parameter is varied, no amount of data on that trajectory will pin
the parameter down, and running expensive MCMC sampling to try is
wasted effort. Checking which parameters a measured output is
sensitive to, before fitting, tells you which parameters are worth
estimating from that data and which should be fixed or determined
separately.

Bioscrape can compute the local sensitivity of a model's trajectories
to its parameters, i.e. the coefficients
:math:`s_{ij} = \partial x_i / \partial p_j` at each simulated time
point, where :math:`x_i` is a state (species) and :math:`p_j` is a
parameter.

Computing Sensitivities
-------------------------

The user-facing entry point is `~bioscrape.analysis.py_sensitivity_analysis`::

    from bioscrape.analysis import py_sensitivity_analysis
    import numpy as np

    timepoints = np.linspace(0, 100, 100)
    S = py_sensitivity_analysis(M, timepoints, normalize=False)

When ``normalize=True``, each coefficient is divided by :math:`x_i /
p_j`, giving a dimensionless relative sensitivity.

Two lower-level helpers are also available for computing sensitivities
at a single point rather than along a full trajectory:

- `~bioscrape.analysis.py_get_jacobian`: the Jacobian :math:`\partial
  f/\partial x` of the model's right-hand side at a given state
- `~bioscrape.analysis.py_get_sensitivity_to_parameter`: the
  sensitivity :math:`\partial f/\partial p` to a single named
  parameter at a given state

Both `~py_sensitivity_analysis` and the two helper functions are
implemented on top of `~bioscrape.analysis.SensitivityAnalysis`, a
`~bioscrape.types.Model` subclass that adds finite-difference
computation of Jacobians and parameter sensitivities using a
deterministic simulation of the model.

Checking Sensitivities Against an Analytical Steady State
------------------------------------------------------------

For a simple enough model, the computed sensitivities can be checked
directly against an analytical result. For the birth-death reaction
below, ``X`` reaches a deterministic steady state of ``k_prod /
k_deg``, and its steady-state sensitivities to the two rate parameters
have simple closed forms: :math:`\partial X/\partial k_{prod} =
1/k_{deg}` and :math:`\partial X/\partial k_{deg} =
-k_{prod}/k_{deg}^2`. Running `~py_sensitivity_analysis` and letting
the trajectory settle shows the computed sensitivities converging to
exactly these values:

.. plot::

    from bioscrape.types import Model
    from bioscrape.analysis import py_sensitivity_analysis
    import numpy as np
    import matplotlib.pyplot as plt

    M = Model()
    M.create_reaction(reactants=[], products=['X'],
                      propensity_type='massaction',
                      propensity_param_dict={'k': 'k_prod'})
    M.create_reaction(reactants=['X'], products=[],
                      propensity_type='massaction',
                      propensity_param_dict={'k': 'k_deg'})
    k_prod, k_deg = 10.0, 0.5
    M.set_parameter('k_prod', k_prod)
    M.set_parameter('k_deg', k_deg)
    M.set_species({'X': 0})

    timepoints = np.linspace(0, 40, 200)
    # S has shape (len(timepoints), len(parameters), len(species)); the
    # parameter order matches M.get_parameter_dictionary(), here
    # [k_prod, k_deg], and there is a single species, X.
    S = py_sensitivity_analysis(M, timepoints, normalize=False)

    plt.plot(timepoints, S[:, 0, 0], label=r'$\partial X/\partial k_{prod}$')
    plt.axhline(1 / k_deg, color='0.5', ls='--')
    plt.plot(timepoints, S[:, 1, 0], label=r'$\partial X/\partial k_{deg}$')
    plt.axhline(-k_prod / k_deg**2, color='0.5', ls='--')
    plt.xlabel('Time')
    plt.ylabel('Sensitivity')
    plt.legend()

For a multi-species example -- the sensitivities of a genetic toggle
switch's four species to all nine of its parameters -- see the
`Sensitivity Analysis Using Bioscrape
<https://github.com/biocircuits/bioscrape/blob/master/examples/Sensitivity%20Analysis%20using%20Bioscrape.ipynb>`__
notebook.

Real-world example
--------------------

The by-parts strategy sketched above -- using sensitivity analysis to
find which parameters a measurement can identify, fitting those,
fixing them, and repeating for the next measured output -- is exactly
how bioscrape's sensitivity and inference tools were used together in
[Pan+23d]_ to characterize a Bxb1 integrase/excisionase DNA
recombination circuit in cell-free extract: sensitivity analysis on
the fluorescence data first identified the small subset of parameters
each reporter's trajectory could actually constrain, and Bayesian
inference (:doc:`parameter_inference`) was then run on only those
parameters, one measured species at a time.

.. [Pan+23d] Pandey A, Rodriguez ML, Poole W, Murray RM (2023)
   Characterization of integrase and excisionase activity in a
   cell-free protein expression system using a modeling and analysis
   pipeline. *ACS Synthetic Biology* 12(2):511-523.
   https://doi.org/10.1021/acssynbio.2c00534
