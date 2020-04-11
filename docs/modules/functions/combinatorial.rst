Combinatorial
=============

This module implements various combinatorial functions.

bell
----

.. autoclass:: diofant.functions.combinatorial.numbers.bell
   :members:

bernoulli
---------

.. autoclass:: diofant.functions.combinatorial.numbers.bernoulli
   :members:

binomial
--------

.. autoclass:: diofant.functions.combinatorial.factorials.binomial
   :members:

catalan
-------

.. autoclass:: diofant.functions.combinatorial.numbers.catalan
   :members:


euler
-----

.. autoclass:: diofant.functions.combinatorial.numbers.euler
   :members:


factorial
---------

.. autoclass:: diofant.functions.combinatorial.factorials.factorial
   :members:

subfactorial
------------

.. autoclass:: diofant.functions.combinatorial.factorials.subfactorial
   :members:

factorial2 / double factorial
-----------------------------

.. autoclass:: diofant.functions.combinatorial.factorials.factorial2
   :members:


FallingFactorial
----------------

.. autoclass:: diofant.functions.combinatorial.factorials.FallingFactorial
   :members:

fibonacci
---------

.. autoclass:: diofant.functions.combinatorial.numbers.fibonacci
   :members:

harmonic
--------

.. autoclass:: diofant.functions.combinatorial.numbers.harmonic
   :members:


lucas
-----

.. autoclass:: diofant.functions.combinatorial.numbers.lucas
   :members:


RisingFactorial
---------------

.. autoclass:: diofant.functions.combinatorial.factorials.RisingFactorial
   :members:

stirling
--------

.. autofunction:: diofant.functions.combinatorial.numbers.stirling

Enumeration
===========

Three functions are available. Each of them attempts to efficiently compute
a given combinatorial quantity for a given set or multiset which can be
entered as an integer, sequence or multiset (dictionary with
elements as keys and multiplicities as values). The ``k`` parameter indicates
the number of elements to pick (or the number of partitions to make). When
``k`` is None, the sum of the enumeration for all ``k`` (from 0 through the
number of items represented by ``n``) is returned. A ``replacement`` parameter
is recognized for combinations and permutations; this indicates that any item
may appear with multiplicity as high as the number of items in the original
set.

.. autofunction:: diofant.functions.combinatorial.numbers.nC

.. autofunction:: diofant.functions.combinatorial.numbers.nP

.. autofunction:: diofant.functions.combinatorial.numbers.nT

Note that the integer for ``n`` indicates *identical* items for ``nT`` but
indicates ``n`` *different* items for ``nC`` and ``nP``.
