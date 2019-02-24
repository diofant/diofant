Logic
=====

.. module:: diofant.logic

Introduction
------------

The logic module for Diofant allows to form and manipulate logic expressions
using symbolic and Boolean values.

Forming logical expressions
---------------------------

You can build Boolean expressions with the standard python operators ``&``
(:class:`~diofant.logic.boolalg.And`), ``|`` (:class:`~diofant.logic.boolalg.Or`),
``~`` (:class:`~diofant.logic.boolalg.Not`)::

    >>> y | (x & y)
    y | (x & y)
    >>> x | y
    x | y
    >>> ~x
    ~x

You can also form implications with ``>>`` and ``<<``::

    >>> x >> y
    Implies(x, y)
    >>> x << y
    Implies(y, x)

Like most types in Diofant, Boolean expressions inherit from
:class:`~diofant.core.basic.Basic`::

    >>> (y & x).subs({x: True, y: True})
    true
    >>> (x | y).atoms()
    {x, y}

The logic module also includes the following functions to derive boolean expressions
from their truth tables-

.. autofunction:: diofant.logic.boolalg.SOPform

.. autofunction:: diofant.logic.boolalg.POSform

Boolean functions
-----------------------

.. autoclass:: diofant.logic.boolalg.BooleanTrue

.. autoclass:: diofant.logic.boolalg.BooleanFalse

.. autoclass:: diofant.logic.boolalg.And

.. autoclass:: diofant.logic.boolalg.Or

.. autoclass:: diofant.logic.boolalg.Not

.. autoclass:: diofant.logic.boolalg.Xor

.. autoclass:: diofant.logic.boolalg.Nand

.. autoclass:: diofant.logic.boolalg.Nor

.. autoclass:: diofant.logic.boolalg.Implies

.. autoclass:: diofant.logic.boolalg.Equivalent

.. autoclass:: diofant.logic.boolalg.ITE

The following functions can be used to handle Conjunctive and Disjunctive Normal
forms-

.. autofunction:: diofant.logic.boolalg.to_cnf

.. autofunction:: diofant.logic.boolalg.to_dnf

.. autofunction:: diofant.logic.boolalg.is_cnf

.. autofunction:: diofant.logic.boolalg.is_dnf

Simplification and equivalence-testing
--------------------------------------

.. autofunction:: diofant.logic.boolalg.simplify_logic

Diofant's simplify() function can also be used to simplify logic expressions to their
simplest forms.

.. autofunction:: diofant.logic.boolalg.bool_map

Inference
---------

.. module: diofant.logic.inference

This module implements some inference routines in propositional logic.

The function satisfiable will test that a given Boolean expression is satisfiable,
that is, you can assign values to the variables to make the sentence `True`.

For example, the expression ``x & ~x`` is not satisfiable, since there are no
values for ``x`` that make this sentence ``True``. On the other hand, ``(x
| y) & (x | ~y) & (~x | y)`` is satisfiable with both ``x`` and ``y`` being
``True``.

    >>> satisfiable(x & ~x)
    False
    >>> satisfiable((x | y) & (x | ~y) & (~x | y))
    {x: True, y: True}

As you see, when a sentence is satisfiable, it returns a model that makes that
sentence ``True``. If it is not satisfiable it will return ``False``.

.. autofunction:: diofant.logic.inference.satisfiable

.. TODO: write about CNF file format
