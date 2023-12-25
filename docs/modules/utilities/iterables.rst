=========
Iterables
=========

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
