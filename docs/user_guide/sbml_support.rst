SBML Support
============

In order to use SBML files, please see the final example in the examples
notebook. This example reads in an SBML file and simulates it using
bioscrape. Additionally, please note that SBML compatibility has the
following limitations.

1. Cannot support delays or events when reading in SBML files. Events
   will be ignored and a warning will be printed out.

2. SBML reaction rates must be in a format such that when the reaction
   rates are converted to a string formula, sympy must be able to parse
   the formula. This will work fine for usual PEMDAS rates. This will
   fail for complex function definitions and things like that.

3. Species will be initialized to their initialAmount field when it is
   nonzero. If the initialAmount is zero, then the initialConcentration
   will be used instead.

4. Any information about compartments will be ignored. Bioscrape will
   only parse in information about reactions, species, and rates. If
   your model has more than one compartment or the compartment changes
   in volume, these facts will be ignored.

5. Assignment rules are supported, but any other type of rule will be
   ignored and an associated warning will be printed out.

6. Parameter names must start with a letter and be alphanumeric, same
   for species names. Furthermore, log, exp, abs, heaviside, and other
   associated keywords for functions are not allowed to be variable
   names. When in doubt, just pick something else.
