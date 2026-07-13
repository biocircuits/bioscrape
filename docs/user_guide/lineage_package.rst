Cell Lineage Simulation
=======================

This documentation covers the cell-population simulation capabilities of
Bioscrape as defined in the bioscrape.lineage package.

Overview:
---------

bioscrape.lineage is a wrapper around the original bioscrape.types and
bioscrape.simulator packages that contains its own sets of objects and
simulators (largely sub-classed from existing objects in bioscrape)
designed to flexibly simulate growing, dividing (including stochastic
partitioning of molecular species), and dying cell populations. Each
individual cell is represented as a bioscrape CRN model with the ability
to grow, divide, and die. Additionally, this package allows large
ensembles of cells to be simulated efficiently offers flexible
statistical wrappers (in development) to ensure that these simulations
are of manageable size.

Simulators
~~~~~~~~~~

bioscrape.lineage supports SSA simulation and SSA with Delay Simulation
of single cells. Cell lineages can be simulated via the following
methods:

-  Lineage Simulation: Beginning with an initial list of cells, the
   ensemble is simulated and all offspring and their time-course data is
   returned.
-  Cell Propagation: Cells are propagated for a set amount of time. All
   cells states (living, dead optional) at the end of the simulation are
   returned. Time-course trajectories not stored.
-  Single Cell Lineage: Beginning with a single cell, this cell is
   simulated with only one daughter cell being saved (chosen randomly)
   each division, resulting in a time-course set of data which tracks
   only a single cell at any time.

Rules
~~~~~

Cell growth, death, and division be modeled as Rule objects, which are
evaluated every simulation step. There are 3 overarching Rule types:
DivisionRules, VolumeRules, DeathRules. A number of different
sub-classes exist

Events
~~~~~~

Alternatively, cell growth, death, and division can be modeled as Event
Objects which fire according to an internal
:doc:`propensity <propensities>`. There are 3 overarching Event types:
DivisionEvents, VolumeEvents, and DeathEvents.

Volume Splitters
~~~~~~~~~~~~~~~~

As cells grow and divide, their internal species are re-partitioned
according to VolumeSplitter objects which are associated with every
DivisionRule and DivisionEvent.
