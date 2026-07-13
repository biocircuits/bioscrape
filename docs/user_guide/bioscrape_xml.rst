XML Specification Language (now deprecated)
===========================================

   Bioscrape's internal XML langauge is now deprecated! Bioscrape now
   uses `Systems Biology Markup
   Language <http://co.mbine.org/standards/sbml>`__ to save CRN models
   except models with Delay. Bioscrape can also load simple SBML models
   consisting of reactions (kinetic laws) with positive rates functions,
   assignment rules, and without compartments. For more information see
   :doc:`SBML Support <sbml_support>` If you stumbled upon this page
   for a reason, it is probably best to jump to the Github issues page
   and let us know why you're here. For posterity's sake, here's the old
   Bioscrape XML description.

Important: Naming Rules
~~~~~~~~~~~~~~~~~~~~~~~

As a rule, names should not contain any strange characters. Basically,
only use alpha numeric characters and underscores in the names of
parameters and species and do not start names with numbers. Do not use
whitespace, colons, semicolons, ^, \|, etc. as part of your species and
parameter names. Additionally, simple functions like abs, log, exp, and
heaviside are reserved keywords for functions, so do not use them as
variable names. Also, the variable t is reserved for time and should not
be used as a variable name. The word volume is reserved for cell volume
and also should not be used as a variable name.

Overview
~~~~~~~~

The XML spec gives a way to specify the reaction related pieces of the
model, i.e. the species, reactions, and rates/propensities and delays
associated with each reaction. The other aspects of the model such a
cell growth and division and type of simulation used
(deterministic/stochastic/etc) all must be specified in your Python
code.

An example XML reaction looks like this.

::

   <reaction text="DNA+RNAP--DNA_RNAP" after="DNA_RNAP--DNA+RNAP+mRNA" >
       <propensity type="massaction" k="k_tx" species="DNA*RNAP" />
       <delay type="fixed" delay="tx_delay" />
   </reaction>

The first part of the reaction tag specifies what happens in the
reaction. The text field specifies what reactants go to what products
immediately. The reactants are split by + signs and the -- separates the
reactants from the products. The ideal delimiter would of course be an
arrow --> but the ">" character is not allowed in XML :-)

The after part of the reaction tag has the same syntax as the text
field, but it contains the stoichiometry that should be executed after
some delay. In this case, when this reaction fires, we will immediately
subtract one DNA and one RNAP and make one DNA_RNAP. After some delay,
we will get rid of that DNA_RNAP and make one mRNA, one RNAP, and one
mRNA. This represents a transcription with delay that models binding of
a single polymerase to each copy of DNA at a time.

The propensity field specifies the reaction rate/propensity. The type
field tells you the type of rate and based on the type of propensity,
there are different additional parameters required. To see all the
possible propensities, please see below. In this case, the propensity is
a mass action propensity, so we have to specify what species are
involved as well as the parameter involved. In this case, the propensity
will have the form k_tx \* DNA \* RNA_P.

The delay field specifies the delay. The type field again tells you the
type of delay. In this case, we have a fixed delay, and the
delay="tx_delay" tells us that parameter tx_delay contains the value of
that fixed delay.

Details on different kinds of :doc:`propensities <propensities>` and
:doc:`delays <delays>` can be found in the propensity functions and
delay functions wiki pages, respectively.

General Equations:
~~~~~~~~~~~~~~~~~~

General propensities and assignment rules make use of general algebraic
equations. The kinds of expressions and functions supported are detailed
below.

Variables in General Equations
------------------------------

The following gives time.

::

   t

The following gives cell volume.

::

   volume

Variables beginning in a "_" (such as \_k) refer to parameters in the
XML model.

Variables which are beginning with a letter (and may contain "*", and
numerals but cannot start with "*" or numerals) refer to species in the
XML model.

Supported Functions in General Equations
----------------------------------------

In addition to the usual parentheses, exponents, mutiplication/division,
addition/subtraction, the following functions are possible.

The following computes the natural logarithm.

::

   log(x)

The following computes the exponential.

::

   exp(x)

The following computes the absolute value. Note the case!! Must be
uppercase A.

::

   Abs(x)

The following computes the maximum or minimum of a bunch of values. Note
the case!! Must be uppercase M.

::

   Max(x,y,z)
   Min(x,y,z)

The two following statements both compute the Heaviside step function.
This is 1 if x is positive and 0 otherwise. Case does not matter here.

::

   heaviside(x)
   Heaviside(x)
