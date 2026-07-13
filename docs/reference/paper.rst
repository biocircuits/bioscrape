Paper and Citation
==================

The canonical software paper for bioscrape is `Fast and flexible simulation and parameter estimation for synthetic biology using bioscrape <https://doi.org/10.21105/joss.05057>`_, published in the Journal of Open Source Software.

Citation
--------

Pandey et al., (2023). Fast and flexible simulation and parameter estimation for synthetic biology using bioscrape. Journal of Open Source Software, 8(83), 5057, https://doi.org/10.21105/joss.05057

BibTeX
------

.. code-block:: bibtex

   @article{Pandey2023,
     doi = {10.21105/joss.05057},
     url = {https://doi.org/10.21105/joss.05057},
     year = {2023},
     publisher = {The Open Journal},
     volume = {8},
     number = {83},
     pages = {5057},
     author = {Pandey, Ayush and Poole, William and Swaminathan, Anandh and Hsiao, Victoria and Murray, Richard M.},
     title = {Fast and flexible simulation and parameter estimation for synthetic biology using bioscrape},
     journal = {Journal of Open Source Software}
   }

Feature Summary
---------------

The paper describes bioscrape as a Cython-based simulator for chemical reaction network models with deterministic, stochastic, delayed-reaction, and growing/dividing-cell simulation capabilities. It also describes the Python API for constructing models programmatically or from SBML, and the parameter inference interface built around repeated forward simulations and MCMC sampling.

The main features highlighted in the paper are:

* Cython-based simulation performance for parameter estimation and stochastic simulation workflows.
* Deterministic simulations, stochastic simulations, stochastic simulation of delayed chemical reaction networks, and simulation of growing and dividing lineages of single cells.
* Markov Chain Monte Carlo based inference tools for estimating parameter distributions from biological circuit data, including common biological data types such as time-series fluorescence and flow cytometry data.
* Local sensitivity analysis for studying model sensitivity to parameters over time.

Source and License
------------------

The JOSS paper page states that the work is licensed under the Creative Commons Attribution 4.0 International License.
