Simulators
==========

The simplest way to simulate a model is the
`~bioscrape.simulator.py_simulate_model` convenience function, which
selects and configures the appropriate simulator based on its keyword
arguments::

    from bioscrape.simulator import py_simulate_model
    import numpy as np

    timepoints = np.linspace(0, 256, 100)
    result = py_simulate_model(timepoints, Model=M, stochastic=True)

By default `~bioscrape.simulator.py_simulate_model` returns a pandas
``DataFrame`` indexed by species and parameter name; pass
``return_dataframe=False`` to get the underlying result object
instead. Different simulation modes can be toggled using its keyword
arguments:

``stochastic = True / False`` Switches between Stochastic (SSA) and
deterministic (ODE) simulation. Defaults to deterministic.

``delay = False / True`` adds delay to stochastic simulation.
Deterministic simulation with delay is not supported. Using delay
requires instantiating reactions with :doc:`delay <delays>`. Defaults
to False.

``volume = False / True`` adds volume parameter for SSA simulations.
For single cell simulation with volume, use the :doc:`Bioscrape
Lineages <lineage_package>` package.

For examples of how to actually do simulations, check `this
notebook <https://github.com/biocircuits/bioscrape/blob/master/examples/Basic%20Examples%20-%20START%20HERE.ipynb>`__.
Advanced users and developers can check out `this
notebook <https://github.com/biocircuits/bioscrape/blob/master/examples/Advanced%20Examples%20-%20Direct%20Simulator%20Instantiation.ipynb>`__
for more details on running simulations more directly.

Deterministic vs. Stochastic Simulation
----------------------------------------

The two agree on average but not trajectory-by-trajectory: for a
birth-death model where a species ``X`` is produced at a constant
rate and degrades in a first-order reaction, a single stochastic
trajectory fluctuates around the deterministic solution rather than
following it exactly, with fluctuations that are largest (relative to
the mean) when molecule counts are small:

.. plot::

    from bioscrape.types import Model
    from bioscrape.simulator import py_simulate_model
    import numpy as np
    import matplotlib.pyplot as plt
    import bioscrape.random
    bioscrape.random.py_seed_random(0)

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

    timepoints = np.linspace(0, 40, 200)
    det = py_simulate_model(timepoints, Model=M, stochastic=False)
    sto = py_simulate_model(timepoints, Model=M, stochastic=True)

    plt.plot(timepoints, sto['X'], color='0.6', label='Stochastic (SSA)')
    plt.plot(timepoints, det['X'], lw=2, label='Deterministic (ODE)')
    plt.xlabel('Time')
    plt.ylabel('X')
    plt.legend()

Delayed Reactions
------------------

Delay changes the shape of the dynamics, not just their timing:
without delay, a reaction's products can appear as soon as the
reaction fires, but with delay, nothing changes until the first
delayed products complete, after which the delayed trajectory
resembles the undelayed one shifted forward in time. For a delayed
transcription model (see :doc:`delays`) where mRNA appears a fixed
time after each firing rather than immediately, comparing simulations
of the same model with and without delay makes this visible directly:

.. plot::

    from bioscrape.types import Model
    from bioscrape.simulator import py_simulate_model
    import numpy as np
    import matplotlib.pyplot as plt
    import bioscrape.random

    M = Model()
    M.create_reaction(
        reactants=[], products=[],
        propensity_type='massaction', propensity_param_dict={'k': 'beta'},
        delay_type='fixed', delay_reactants=[], delay_products=['mRNA'],
        delay_param_dict={'delay': 'tx_delay'})
    M.create_reaction(
        reactants=['mRNA'], products=[],
        propensity_type='massaction', propensity_param_dict={'k': 'delta'})
    M.set_parameter('beta', 2)
    M.set_parameter('delta', 0.2)
    M.set_parameter('tx_delay', 10)
    M.set_species({'mRNA': 0})

    timepoints = np.linspace(0, 50, 200)

    # Re-seeding before each run means both share the same sequence of
    # reaction firings, isolating the effect of the delay itself.
    bioscrape.random.py_seed_random(0)
    no_delay = py_simulate_model(timepoints, Model=M, stochastic=True, delay=False)
    bioscrape.random.py_seed_random(0)
    with_delay = py_simulate_model(timepoints, Model=M, stochastic=True, delay=True)

    plt.plot(timepoints, no_delay['mRNA'], label='No delay')
    plt.plot(timepoints, with_delay['mRNA'], label='Fixed delay = 10')
    plt.xlabel('Time')
    plt.ylabel('mRNA')
    plt.legend()

Volume-Dependent Rates
------------------------

When a volume is present (``volume=True``), species are counts of
molecules but the reaction rates depend on their *concentrations*, so
bioscrape rescales the propensities by the current volume. This
happens automatically; you write the same rate laws as in a
volumeless simulation, and the following adjustments are applied
under the hood:

-  A *zeroth-order* (constitutive) mass-action rate scales with the
   volume: a fixed concentration of product is made per unit time, so
   the number of molecules produced grows with the cell.
-  A *first-order* (unimolecular) mass-action rate is unchanged: it
   depends only on the count of the single reactant.
-  A *second-order* (bimolecular) mass-action rate scales as
   ``1 / volume``: two molecules are less likely to find each other in
   a larger volume.
-  More generally, an n-th order mass-action rate scales as
   ``1 / volume ** (n - 1)``.
-  In a Hill propensity, the input species count is divided by the
   volume, i.e. the rate is computed from the species' *concentration*
   rather than its count, because it is the concentration that sets
   the equilibrium binding (for example of a transcription factor to a
   promoter).

Because these conversions are automatic, a rate constant fit or chosen
for a volumeless simulation may need to be re-scaled when the same
model is simulated with a volume, and vice versa. Single-cell
simulation with volume and cell division uses the :doc:`Bioscrape
Lineages <lineage_package>` package.

Checking a Stochastic Simulation
-----------------------------------

For a simple enough model, the stochastic simulator's output can be
checked directly against a known analytical result, which is a useful
way to build confidence in a new model or a modified reaction network.
The birth-death model above is a textbook example: at steady state,
its species count ``X`` follows a Poisson distribution with mean and
variance both equal to ``k_prod / k_deg``. Simulating it stochastically
for long enough, and comparing the distribution of ``X`` values (after
discarding an initial transient) against that Poisson distribution,
confirms that bioscrape's SSA implementation reproduces the expected
statistics:

.. plot::

    from bioscrape.types import Model
    from bioscrape.simulator import py_simulate_model
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import poisson
    import bioscrape.random
    bioscrape.random.py_seed_random(0)

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
    result = py_simulate_model(timepoints, Model=M, stochastic=True)
    X = result['X'].values

    # Discard the transient before the distribution reaches steady state.
    X_steady_state = X[timepoints > 200].astype(int)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.plot(timepoints[timepoints < 100], X[timepoints < 100])
    ax1.set_xlabel('Time')
    ax1.set_ylabel('X')
    ax1.set_title('Stochastic trajectory')

    counts = np.bincount(X_steady_state)
    xs = np.arange(len(counts))
    ax2.bar(xs, counts / counts.sum(), color='0.7', label='Empirical')
    ax2.plot(xs, poisson.pmf(xs, mu=10.0 / 0.5), lw=2, label='Poisson(20)')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Probability')
    ax2.set_title('Steady-state distribution')
    ax2.legend()

    plt.tight_layout()

`~bioscrape.simulator.SSAResult` provides some of this analysis
directly, including `~bioscrape.simulator.SSAResult.py_first_moment`
and `~bioscrape.simulator.SSAResult.py_correlations` for computing
means and correlations from a simulated trajectory without the manual
post-processing shown above. For a more complete worked example,
including statistics of a multi-species model, see the *Fast
Statistics for SSA Trajectories* notebook.

Built-in Simulators
---------------------

Simulations can be performed deterministically (without delay or cell
division, i.e. regular deterministic ODE's for a system of reactions),
or they can be performed stochastically with or without delay and with
or without cell growth and division using the :doc:`Bioscrape
Lineages <lineage_package>` package. The concrete simulator classes,
selected automatically by `~bioscrape.simulator.py_simulate_model`
based on its ``stochastic``/``delay``/``volume`` arguments, are:

-  `~bioscrape.simulator.DeterministicSimulator`
-  `~bioscrape.simulator.SSASimulator`
-  `~bioscrape.simulator.DelaySSASimulator`
-  `~bioscrape.simulator.VolumeSSASimulator`
-  `~bioscrape.simulator.DelayVolumeSSASimulator`

Simulation Interfaces
------------------------

`~bioscrape.simulator.CSimInterface` and its subclasses
(`~bioscrape.simulator.ModelCSimInterface`,
`~bioscrape.simulator.SafeModelCSimInterface`) adapt a model into the
low-level form used internally by the simulators. Most users will not
need to construct these directly --
`~bioscrape.simulator.py_simulate_model` creates one automatically --
but they can be reused across repeated simulations of the same model
for performance (the ``Interface`` keyword argument), or replaced with
a `~bioscrape.simulator.SafeModelCSimInterface` (via ``safe=True``) to
get warnings about ill-conditioned situations such as negative
propensities.
