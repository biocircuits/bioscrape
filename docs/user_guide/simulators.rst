Simulators
==========

For examples of how to actually do simulations, check `this
notebook <https://github.com/biocircuits/bioscrape/blob/master/examples/Basic%20Examples%20-%20START%20HERE.ipynb>`__.
Briefly, different simulation modes can be easily toggled using the
following keywords in the ``py_SimulateModel(timepoints, Model = ...)``
wrapper:

``stochastic = True / False`` Switches between Stochastic (SSA) and
determinstic (ODE simulation). Defaults to deterministic.

``delay = False / True`` adds delay to stochastic simulation.
Deterministic simulation with delay is not supported. Using delay
requires instantiating reactions with :doc:`delay <delays>`. Defaults
to False.

``volume = False / True`` adds volume parameter for SSA simulations. For
single cell simulation with volume, use the (lineages
package)[Lineage-Package]

Advanced users and developers can check out `this
notebook <https://github.com/biocircuits/bioscrape/blob/master/examples/Advanced%20Examples%20-%20Direct%20Simulator%20Instantiation.ipynb>`__
for more details on running simulations more directly.

Built-in Simulators
-------------------

Simulations can be performed deterministically (without delay or cell
division, i.e. regular deterministic ODE's for a system of reactions),
or they can be performed stochastically with or without delay and with
or without cell growth and division using the :doc:`Bioscrape
Lineages <lineage_package>` package.
