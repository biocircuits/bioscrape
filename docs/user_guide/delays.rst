Delays in Bioscrape
===================

In stochastic simulations, bioscrape supports delay. In a delay
reaction, delay inputs/outputs are consumed/produced after some amount
of delay. Reactions may have a mix of delay and non-delay inputs and
outputs. Bioscrape innately supports a number of delay-types:

Note that delays and delay inputs/outputs will be ignored if a model
with delays is simulated deterministically.

Types of delays
---------------

1. fixed: constant delay with parameter "delay".
2. Gaussian: gaussian distributed delay with parameters "mean" and
   "std".
3. Gamma: gamma distributed delay with shape parameter "k" and scale
   parameter "theta".

Delays for each reaction can be specified in a Bioscrape model.

No Delay
--------

By default, all reactions are created without a delay.

Fixed Delay
-----------

The reaction below creates a product "B" from "A" with a ``massaction``
propensity after a fixed delay.

::

   rxn = (["A"], ["A"], "massaction", {"k":2},
          "fixed", [], ["T"], {"delay":"k_delay"})

If the delay is a fixed constant, then set the delay to a fixed constant
specified by the parameter "k_delay" in the model.

Gaussian Delay
--------------

The reaction below creates a product "T" from "G" with a
``proportionalhillpositive`` propensity after a delay that has a
Gaussian distribution.

::

   rxn = (["G"], ["G"], "proportionalhillpositive", {"d":"G", "s1":"I", "k":"ktx", "K":"KI", "n":"n"},
          "gaussian", [], ["T"], {"mean":10.0, "std":1.0})

If the delay is Gaussian, then it is distributed according to a normal
distribution with mean parameter mean and standard deviation parameter
std.

Gamma Distributed Delay
-----------------------

The following line creates a reaction that forms a product "X" from "T"
with a ``hillpositive`` propensity after a delay that is gamma
distributed.

::

   rxn2d = (["T"], ["T"], "hillpositive", {"s1":"T", "k":"ktl", "K":"KR", "n":1},
           "gamma", [], ["X"], {"k":10.0, "theta":3.0})

If the delay is distributed with a gamma distribution as is commonly the
case in biological circuits, then you can specify that the delay type is
gamma with "k" and "theta" being parameters that contain the k and theta
parameters of the gamma distribution.

Advanced: The Delay Queue
-------------------------

To do simulations with delay, you may also create a delay queue object
that stores reactions that have been queued up to happen in the future.
This is part of the initial state because the initial state is both the
current species level as well as any reactions that have been queued up
to occur. Currently, the only type of delay queue object available is
called an ArrayDelayQueue and one can be created with the code

::

   q = ArrayDelayQueue.setup_queue(num_reactions, num_timepoints, dt)

This will create a delay queue that is capable of storing events up to
num_timepoints*dt time units into the future with a temporal resolution
of dt. Making dt smaller makes your timing more accurate, but will make
your code slower. It will especially slow things down when you simulate
cell division. You should select a suitable dt, and then select
num_timepoints to be such that the maximum possible storeable delay is
big enough for your purposes. The num_reactions parameter should be
picked based on your system of reactions. If you want to add reactions
to the queue, you can do so using the following code.

::

   q.py_add_reaction(time, rxn_id, amount)

The time is the time that the reaction occurs, the reaction ID is the
index of that reaction in the stoichiometric matrix, and the amount is
the number of times you want that reaction to fire at the specified
time. This would typically be 1.0 since the reaction would usually fire
once in a stochastic simulation, but in some cases you might want it to
fire multiple times. If you want to set up a delay simulation with
reactions queued up as part of the initial state, you can create a queue
and use py_add_reaction() to pre-queue reactions that will occur.
