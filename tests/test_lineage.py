import pytest

#The lineage extension is optional and separately built
#  (python setup.py build_ext --inplace lineage), so skip this whole
#  module if it isn't available rather than erroring.
lineage = pytest.importorskip("bioscrape.lineage")

import numpy as np
from bioscrape.lineage import (LineageModel, LineageVolumeSplitter,
    LineageCSimInterface, LineageSSASimulator,
    InteractingLineageSSASimulator, LineageVolumeCellState,
    py_SimulateSingleCell)
from bioscrape.simulator import ModelCSimInterface, VolumeSSASimulator
from bioscrape.types import Model


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


#This test confirms that InteractingLineageSSASimulator.setup_global_
#  volume_simulation doesn't raise AttributeError. It previously did
#  global_species_global_crn_inds.astype(np.int, copy=False) --
#  np.int was removed in NumPy >= 1.24 -- and also discarded the
#  astype() return value instead of assigning it back.
def test_setup_global_volume_simulation():
    global_model = Model(species = ["G"], reactions = [],
                          initial_condition_dict = {"G": 5})
    global_interface = ModelCSimInterface(global_model)
    global_volume_simulator = VolumeSSASimulator()
    inds = np.array([0], dtype = np.int32)

    simulator = InteractingLineageSSASimulator()
    simulator.setup_global_volume_simulation(global_volume_simulator,
                                              global_interface, inds)


#This test confirms that py_SimulateSingleCell actually uses a
#  passed-in initial_cell_state. It previously only handled the
#  initial_cell_state=None case, leaving the local variable v
#  unassigned (and thus raising NameError) for any other value.
def test_simulate_single_cell_initial_cell_state():
    M = _make_minimal_lineage_model()
    interface = LineageCSimInterface(M)
    interface.py_set_initial_time(0.0)
    initial_cell_state = LineageVolumeCellState(v0 = 1, t0 = 0,
                                                 state = interface.py_get_initial_state())
    timepoints = np.linspace(0, 1, 5)

    result = py_SimulateSingleCell(timepoints, interface = interface,
                                    initial_cell_state = initial_cell_state,
                                    return_dataframes = False)
    assert result is not None


#This test confirms that LineageVolumeSplitter's constructor can set
#  up a custom default splitter without crashing. It previously did
#  self.ind2customsplitter == {} (a comparison, typo for =), leaving
#  ind2customsplitter as None and raising TypeError on item
#  assignment the first time a custom default splitter was used.
def test_lineage_volume_splitter_custom_default():
    M = _make_minimal_lineage_model()
    def my_splitter(parent, which):
        return 1.0, 1.0
    vsplit = LineageVolumeSplitter(M,
        options = {"default": "my_custom", "volume": "my_custom"},
        custom_partition_functions = {"my_custom": my_splitter})
    assert vsplit is not None


#This test confirms that LineageVolumeSplitter validates
#  partition_noise is between 0 and 1. It previously compared
#  self.partition_noise (still its cdef zero-default) before the
#  constructor argument was ever assigned, so the check never fired
#  for any input.
def test_lineage_volume_splitter_partition_noise_bounds():
    M = _make_minimal_lineage_model()
    with pytest.raises(ValueError):
        LineageVolumeSplitter(M, partition_noise = 5.0)
    with pytest.raises(ValueError):
        LineageVolumeSplitter(M, partition_noise = -1.0)


#These tests confirm that LineageModel's rules=/events= constructor
#  shortcuts correctly route short-form type aliases (e.g. "linear",
#  "species", "deltaV") to create_volume_rule/create_death_rule/
#  create_division_rule, not just fully-qualified names like
#  "linearvolumerule" containing the literal word "volume"/"death"/
#  "division". Previously, short-form aliases fell through the
#  substring-based routing entirely and were silently forwarded to
#  the base Model constructor as unrecognized rule/event types.
def test_lineage_model_rules_short_form_volume():
    species = ["X"]
    x0 = {"X": 10}
    rxn = [[], ["X"], "massaction", {"k": 1.0}]
    M = LineageModel(species = species, reactions = [rxn],
                      initial_condition_dict = x0,
                      rules = [("linear", {"growth_rate": 0.1})])
    assert M.py_get_num_volume_rules() == 1

def test_lineage_model_rules_short_form_death():
    species = ["X"]
    x0 = {"X": 10}
    rxn = [[], ["X"], "massaction", {"k": 1.0}]
    M = LineageModel(species = species, reactions = [rxn],
                      initial_condition_dict = x0,
                      rules = [("species", {"specie": "X", "threshold": 5})])
    assert M.py_get_num_death_rules() == 1

#Division rules still can't be constructed via the rules= shortcut
#  (create_division_rule requires an explicit VolumeSplitter that
#  the shortcut has no way to supply -- a separate, documented
#  limitation), but a short-form division alias like "deltaV" should
#  now at least be correctly identified as a division rule (raising
#  a clear TypeError about the missing volume_splitter) instead of
#  being silently misrouted to the base Model as an unrecognized
#  rule type.
def test_lineage_model_rules_short_form_division_routes_correctly():
    species = ["X"]
    x0 = {"X": 10}
    rxn = [[], ["X"], "massaction", {"k": 1.0}]
    with pytest.raises(TypeError):
        LineageModel(species = species, reactions = [rxn],
                      initial_condition_dict = x0,
                      rules = [("deltaV", {"threshold": 1.0})])

def test_lineage_model_events_short_form_volume():
    species = ["X"]
    x0 = {"X": 10}
    rxn = [[], ["X"], "massaction", {"k": 1.0}]
    M = LineageModel(species = species, reactions = [rxn],
                      initial_condition_dict = x0,
                      events = [("linear", {"growth_rate": 0.1},
                                 "massaction", {"k": 0.1, "species": ""})])
    assert M.py_get_num_volume_events() == 1
