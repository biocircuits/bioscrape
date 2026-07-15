import pytest

#The lineage extension is optional and separately built
#  (python setup.py build_ext --inplace lineage), so skip this whole
#  module if it isn't available rather than erroring.
lineage = pytest.importorskip("bioscrape.lineage")

import numpy as np
from bioscrape.lineage import (LineageModel, LineageVolumeSplitter,
    LineageCSimInterface, LineageSSASimulator)


def _make_minimal_lineage_model():
    species = ["X"]
    x0 = {"X": 10}
    rxn = [[], ["X"], "massaction", {"k": 1.0}]
    M = LineageModel(species = species, reactions = [rxn],
                      initial_condition_dict = x0)
    vsplit = LineageVolumeSplitter(M)
    M.create_division_rule("deltaV", {"threshold": 1.0}, vsplit)
    M.create_volume_event("linear volume", {"growth_rate": 0.1},
                           "massaction", {"k": 0.1, "species": ""})
    M.py_initialize()
    return M


#This test confirms that py_initialize_single_cell_interface doesn't
#  raise AttributeError. It previously called
#  self.initialize_single_cell_interface(...), but the actual cdef
#  method was misspelled intialize_single_cell_interface -- since
#  fixed by correcting the cdef method's name (and all its other call
#  sites) to match, rather than matching the call site to the typo.
def test_py_initialize_single_cell_interface():
    M = _make_minimal_lineage_model()
    interface = LineageCSimInterface(M)
    interface.py_set_initial_time(0.0)
    simulator = LineageSSASimulator()
    simulator.py_initialize_single_cell_interface(interface)
