from bioscrape.simulator import *
from bioscrape.types import *
import numpy as np
import pytest
import sys

timepoints = np.arange(0, 100, .1)

#This test confirms that calling model.py_initialize twice before simulation does not cause bugs
#  This test was written due to a bug found where the c_propensities vector was not cleared between initializations.
def test_initialize():
    rxn1 = (["A", "A"], ["B"], "massaction", {"k":1.0})
    rxn2 = (["B"], [], "massaction", {"k":.1})
    rxn3 = (["A", "B"], ["C"], "massaction", {"k":2.0})

    #create a model which will be automatically initialized
    M1 = Model(species = ["A", "B", "C"], reactions = [rxn1, rxn2, rxn3], initial_condition_dict = {"A":100})
    R1 = py_simulate_model(Model = M1, timepoints = timepoints)

    #Create another model which will be initialied twice
    M2 = Model(species = ["A", "B", "C"], reactions = [rxn2, rxn3], initial_condition_dict = {"A":100})
    M2.create_reaction(rxn1[0], rxn1[1], rxn1[2], rxn1[3]) #this will uninitialize the model
    #Simulating it will cause it to initialize again
    R2 = py_simulate_model(Model = M2, timepoints = timepoints)

    #the models are identical, so R1 and R2 should have the same output
    assert np.allclose(R1["A"], R2["A"])
    assert np.allclose(R1["B"], R2["B"])
    assert np.allclose(R1["C"], R2["C"])

#This test confirms that Schnitz.py_get_dataframe falls back to a numpy
#  array (instead of raising AttributeError) when pandas is unavailable.
def test_schnitz_dataframe_pandas_fallback():
    time = np.array([0.0, 1.0, 2.0])
    data = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    volume = np.array([1.0, 1.0, 1.0])
    sch = Schnitz(time, data, volume)

    #Setting sys.modules["pandas"] = None makes the import system raise
    #  ModuleNotFoundError, regardless of whether the importing code is
    #  Python or Cython (Cython's `import` bypasses builtins.__import__,
    #  so mocking that hook doesn't work here).
    real_pandas = sys.modules.get("pandas", None)
    sys.modules["pandas"] = None
    try:
        with pytest.warns(UserWarning):
            result = sch.py_get_dataframe()
    finally:
        if real_pandas is not None:
            sys.modules["pandas"] = real_pandas
        else:
            del sys.modules["pandas"]

    assert np.array_equal(result, data)

#This test confirms that py_simulate_model(delay=True, volume=True)
#  doesn't raise NotImplementedError. It previously instantiated the
#  abstract DelayVolumeSimulator instead of the concrete
#  DelayVolumeSSASimulator.
def test_simulate_delay_and_volume():
    species = ["A", "B"]
    x0 = {"A": 10, "B": 0}
    rxn = (["A"], [], "massaction", {"k": 1.0}, "fixed", [], ["B"],
           {"delay": 1.0})
    M = Model(species = species, reactions = [rxn],
              initial_condition_dict = x0)
    delay_timepoints = np.linspace(0, 10, 11)

    result = py_simulate_model(delay_timepoints, Model = M,
                                stochastic = True, delay = True,
                                volume = True, return_dataframe = False)

    assert result.py_get_volume().shape == delay_timepoints.shape

