from importlib.metadata import version 
import bioscrape.random
import bioscrape.types
import bioscrape.simulator
import bioscrape.inference

bioscrape.random.py_seed_random()
__version__ = version("Bioscrape")