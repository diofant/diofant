=========
Iterables
=========

variations
----------

variations(seq, n) Returns all the variations of the list of size n.

Has an optional third argument. Must be a boolean value and makes the method
return the variations with repetition if set to True, or the variations
without repetition if set to False.

Examples::
    >>> from diofant.utilities.iterables import variations
    >>> list(variations([1, 2, 3], 2))
    [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
    >>> list(variations([1, 2, 3], 2, True))
    [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]


partitions
----------

Although the combinatorics module contains Partition and IntegerPartition
classes for investigation and manipulation of partitions, there are a few
functions to generate partitions that can be used as low-level tools for
routines:  ``partitions`` and ``multiset_partitions``. The former gives
integer partitions, and the latter gives enumerated partitions of elements.

partitions::

    >>> from diofant.utilities.iterables import partitions
    >>> [p.copy() for s, p in partitions(7, m=2, size=True) if s == 2]
    [{1: 1, 6: 1}, {2: 1, 5: 1}, {3: 1, 4: 1}]

multiset_partitions::

    >>> from diofant.utilities.iterables import multiset_partitions
    >>> list(multiset_partitions(3, 2))
    [[[0, 1], [2]], [[0, 2], [1]], [[0], [1, 2]]]
    >>> list(multiset_partitions([1, 1, 1, 2], 2))
    [[[1, 1, 1], [2]], [[1, 1, 2], [1]], [[1, 1], [1, 2]]]

Docstring
=========

.. automodule:: diofant.utilities.iterables
   :members:
