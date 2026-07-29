Bioscrape - Biological Modeling, Simulation, and Inference
==========================================================

Bioscrape is a Python package for deterministic, stochastic, delayed,
and cell-lineage simulation of biological chemical reaction network
models. It also provides parameter inference tools that can be applied
to estimate model parameters from experimental data.

.. rubric:: Features

- Fast Cython simulators for deterministic, stochastic, delayed, and
  cell-lineage models.
- SBML support for loading, saving, and exchanging biological chemical
  reaction network models.
- Python APIs for constructing models, propensities, delays, rules, and
  simulation workflows.
- Bayesian and deterministic parameter inference interfaces built on the
  Bioscrape simulation stack.

.. rubric:: Links

- Source code: https://github.com/biocircuits/bioscrape
- Bug reports: https://github.com/biocircuits/bioscrape/issues
- Mailing list: `SBTools Google Group <https://groups.google.com/g/sbtools/>`_
- Citation details: :doc:`reference/paper`

.. toctree::
   :caption: User Guide
   :maxdepth: 1
   :numbered: 2

   user_guide/overview
   user_guide/installation
   user_guide/examples
   user_guide/simulators
   user_guide/propensities
   user_guide/delays
   user_guide/sensitivity
   user_guide/parameter_inference
   user_guide/sbml_support
   user_guide/lineage_package

.. toctree::
   :caption: Tutorial Examples
   :maxdepth: 1

   Tutorial Examples <examples/index>

.. toctree::
   :caption: Reference Manual
   :maxdepth: 1

   api/index
   reference/paper
   developer/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
