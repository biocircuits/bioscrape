Propensities
============

Types of propensities
---------------------

Below the different built in propensity types are outlined. These
propensities can also be called using the model API
Model.create_propensity and Model.create_reaction functions. These
functions both take the argument propensity_type which is a string of
the type attribute of the propensity and propensity_param_dict which is
a dictionary containing all other attributes. For example,
propensity_type for the mass action propensity below is "massaction" and
the propensity_param_dict would be ``{"k":"k1", "species":"X*Y*Z"}``.

Reaction with Mass-Action Propensity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

       reaction = Reaction(list_of_reactants, list_of_products, "massaction", {"k":"k1"})

A mass action propensity requires one constant parameter k for each
irreversible reaction step. Optionally, a "species" keyword can be
passed in which changes the massaction propensity, ``species = "X*Y"``.
If a species should be raised to a power, then just do something like
X\ *X*\ Y. The final rate will then be k1\ *X*\ X*Y for this example. If
you want a constitutive rate, then ignore the species key.

If there is volume, then the rate will scale appropriately with the
number of species. For example, a constitutive reaction's rate will
scale proportionally with volume, a unimolecular reaction's rate will be
independent of volume, and a bimolecular reaction's rate will scale as 1
/ volume.

Reaction with Positive Hill Propensity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

       reaction = Reaction(list_of_reactants, list_of_products, "hillpositive", {"s1"="X" "k"="k1" "K"="K1" "n"="n1"})

A positive hill function propensity with input species s1, maximal rate
k, cooperativity n, centered at K. For example, the above propensity
would result in the general form rate(X) = k \* X^n1 / (K^n1 + X^n1).

Reaction with Negative Hill Propensity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

       reaction = Reaction(list_of_reactants, list_of_products, "hillnegative", {"s1"="X" "k"="k1" "K"="K1" "n"="n1"})

A negative hill function propensity with input species s1, maximal rate
k, cooperativity n, centered at K. For example, the above propensity
would result in the general form rate(X) = k \* 1 / (K^n1 + X^n1).

Reaction with Proportional Positive Hill Propensity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

       reaction = Reaction(list_of_reactants, list_of_products, "proportionalhillpositive", {"d"="Y" "s1"="X" "k"="k1" "K"="K1" "n"="n1"})

A positive hill function propensity proportional to species d, with
saturating input species s1, maximal rate k, cooperativity n, centered
at K. For example, the above propensity would result in the general form
rate(X) = k \* Y \* X^n1 / (K^n1 + X^n1).

Proportional Negative Hill Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

       reaction = Reaction(list_of_reactants, list_of_products, "proportionalhillnegative", {"d"="Y" "s1"="X" "k"="k1" "K"="K1" "n"="n1"})

A negative hill function propensity proportional to species d, with
desaturating input species s1, maximal rate k, cooperativity n, centered
at K. For example, the above propensity would result in the general form
rate(X) = k \* Y / (K^n1 + X^n1).

General Propensities
~~~~~~~~~~~~~~~~~~~~

If these propensity types don't work for you, you can also create any
custom combination of species, parameters, numbers, and volume that you
can assemble together with standard parentheses, exponents,
addition/subtraction and multiplication/division. The general
propensities also can handle functions of time as well as step
functions, log, exp, Min, and Max (note capitalization). However, using
this type of propensity might be somewhat slower, so always use the
massaction type when you can if you want fast simulations, especially
with stochastic and cell lineage simulations. Here are a few examples of
general propensity reactions:

::

       reaction = Reaction(list_of_reactants, list_of_products, "general", {"rate":" 0.15 + X / (X+_Kx) * Y^_n / (Y^_n+_Ky^_n)"})
       reaction = Reaction(list_of_reactants, list_of_products, "general", {"rate":" log(exp(heaviside(t-4)*(t-4)))*beta"})
       reaction = Reaction(list_of_reactants, list_of_products, "general", {"rate":" (x+1)**2-x**2-2*x-1+t*alpha"})
       reaction = Reaction(list_of_reactants, list_of_products, "general", {"rate":" Min(v1*X, v2*Y)"})

The rate field contains the entire form of the propensity. In this case,
we have a small leak of 0.15 added to an AND gate on two species X and
Y, where both X and Y must be high for this reaction to fire. The top
line specifies three parameters Kx, Ky, and n, and creates a sum of 0.15
with a product of Hill functions in X and Y where the Hill coefficient
for X is 1, but the Hill coefficient for Y is a parameter n. The next
two lines show the more complex expressions that are possible.

This type of rate can be super general, but it can be about 50% slower
than the above rates, so try and use it only when you really need to use
it.

The keyword volume refers to the volume. This term is evaluated to the
volume if you're doing a volume based simulation or 1.0 otherwise. For
example,

::

       volume * X

would be just X when doing a volume less simulation or the volume
multiplied by X if you're doing a simulation with volume.

Also, the keyword t refers to time. An expression like

::

       heaviside(t-60.0)*beta

would enforce a reaction to have no rate until time reaches 60, after
which it has a constant rate of beta. The heaviside function allows you
to essentially model events within the reaction rate.
