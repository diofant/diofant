Core
====

.. automodule:: diofant.core

sympify
-------
.. module:: diofant.core.sympify

sympify
^^^^^^^
.. autofunction:: sympify

assumptions
-----------

.. automodule:: diofant.core.assumptions
    :members:

cache
-------
.. module:: diofant.core.cache

cacheit
^^^^^^^
.. autofunction:: cacheit

basic
-----
.. automodule:: diofant.core.basic
   :members:

core
----
.. module:: diofant.core.core

singleton
---------
.. automodule:: diofant.core.singleton
   :members:

evaluate
--------

.. automodule:: diofant.core.evaluate
   :members:

expr
----
.. module:: diofant.core.expr

Expr
----
.. autoclass:: Expr
   :members:

AtomicExpr
----------
.. autoclass:: AtomicExpr
   :members:

symbol
------
.. module:: diofant.core.symbol

Symbol
^^^^^^
.. autoclass:: Symbol
   :members:

Wild
^^^^
.. autoclass:: Wild
   :members:

Dummy
^^^^^
.. autoclass:: Dummy
   :members:

symbols
^^^^^^^
.. autofunction:: symbols

var
^^^
.. autofunction:: var

numbers
-------
.. module:: diofant.core.numbers

Number
^^^^^^
.. autoclass:: Number
   :members:

Float
^^^^^
.. autoclass:: Float
   :members:

Rational
^^^^^^^^
.. autoclass:: Rational
   :members:

Integer
^^^^^^^
.. autoclass:: Integer
   :members:

NumberSymbol
^^^^^^^^^^^^
.. autoclass:: NumberSymbol
   :members:

mod_inverse
^^^^^^^^^^^
.. autofunction:: mod_inverse

integer_digits
^^^^^^^^^^^^^^
.. autofunction:: integer_digits

Zero
^^^^

.. autoclass:: Zero
   :members:

One
^^^

.. autoclass:: One
   :members:

NegativeOne
^^^^^^^^^^^

.. autoclass:: NegativeOne
   :members:

Half
^^^^

.. autoclass:: Half
   :members:

NaN
^^^

.. autoclass:: NaN
   :members:

Infinity
^^^^^^^^

.. autoclass:: Infinity
   :members:

NegativeInfinity
^^^^^^^^^^^^^^^^

.. autoclass:: NegativeInfinity
   :members:

ComplexInfinity
^^^^^^^^^^^^^^^

.. autoclass:: ComplexInfinity
   :members:

Exp1
^^^^

.. autoclass:: Exp1
   :members:

ImaginaryUnit
^^^^^^^^^^^^^

.. autoclass:: ImaginaryUnit
   :members:

Pi
^^

.. autoclass:: Pi
   :members:

EulerGamma
^^^^^^^^^^

.. autoclass:: EulerGamma
   :members:

Catalan
^^^^^^^

.. autoclass:: Catalan
   :members:

GoldenRatio
^^^^^^^^^^^

.. autoclass:: GoldenRatio
   :members:

power
-----
.. module:: diofant.core.power

Pow
^^^
.. autoclass:: Pow
   :members:

integer_nthroot
^^^^^^^^^^^^^^^
.. autofunction:: integer_nthroot

mul
---
.. module:: diofant.core.mul

Mul
^^^
.. autoclass:: Mul
   :members:

add
---
.. module:: diofant.core.add

Add
^^^
.. autoclass:: Add
   :members:

mod
---
.. module:: diofant.core.mod

Mod
^^^
.. autoclass:: Mod
   :members:

relational
----------
.. module:: diofant.core.relational

Rel
^^^
.. autoclass:: Rel
   :members:

Eq
^^
.. autoclass:: Eq
   :members:

Ne
^^
.. autoclass:: Ne
   :members:

Lt
^^
.. autoclass:: Lt
   :members:

Le
^^
.. autoclass:: Le
   :members:

Gt
^^
.. autoclass:: Gt
   :members:

Ge
^^
.. autoclass:: Ge
   :members:

Relational
^^^^^^^^^^
.. autoclass:: Relational
   :members:

Equality
^^^^^^^^
.. autoclass:: Equality
   :members:

GreaterThan
^^^^^^^^^^^
.. autoclass:: GreaterThan
   :members:

LessThan
^^^^^^^^
.. autoclass:: LessThan
   :members:

Unequality
^^^^^^^^^^
.. autoclass:: Unequality
   :members:

StrictGreaterThan
^^^^^^^^^^^^^^^^^
.. autoclass:: StrictGreaterThan
   :members:

StrictLessThan
^^^^^^^^^^^^^^
.. autoclass:: StrictLessThan
   :members:

multidimensional
----------------
.. module:: diofant.core.multidimensional

vectorize
^^^^^^^^^
.. autoclass:: vectorize
   :members:

function
--------
.. module:: diofant.core.function

Lambda
^^^^^^
.. autoclass:: Lambda
   :members:

WildFunction
^^^^^^^^^^^^
.. autoclass:: WildFunction
   :members:

Derivative
^^^^^^^^^^
.. autoclass:: Derivative
   :members:

diff
^^^^
.. autofunction:: diff

FunctionClass
^^^^^^^^^^^^^
.. autoclass:: FunctionClass
   :members:

Function
^^^^^^^^
.. autoclass:: Function
   :members:

.. note:: Not all functions are the same

   Diofant defines many functions (like ``cos`` and ``factorial``). It also
   allows the user to create generic functions which act as argument
   holders. Such functions are created just like symbols:

   >>> f = Function('f')
   >>> f(2) + f(x)
   f(2) + f(x)

   If you want to see which functions appear in an expression you can use
   the atoms method:

   >>> e = (f(x) + cos(x) + 2)
   >>> e.atoms(Function)
   {f(x), cos(x)}

   If you just want the function you defined, not Diofant functions, the
   thing to search for is AppliedUndef:

   >>> from diofant.core.function import AppliedUndef
   >>> e.atoms(AppliedUndef)
   {f(x)}

Subs
^^^^
.. autoclass:: Subs
   :members:

expand
^^^^^^
.. autofunction:: expand

PoleError
^^^^^^^^^
.. autoclass:: PoleError
   :members:

count_ops
^^^^^^^^^
.. autofunction:: count_ops

expand_mul
^^^^^^^^^^
.. autofunction:: expand_mul

expand_log
^^^^^^^^^^
.. autofunction:: expand_log

expand_func
^^^^^^^^^^^
.. autofunction:: expand_func

expand_trig
^^^^^^^^^^^
.. autofunction:: expand_trig

expand_complex
^^^^^^^^^^^^^^
.. autofunction:: expand_complex

expand_multinomial
^^^^^^^^^^^^^^^^^^
.. autofunction:: expand_multinomial

expand_power_exp
^^^^^^^^^^^^^^^^
.. autofunction:: expand_power_exp

expand_power_base
^^^^^^^^^^^^^^^^^
.. autofunction:: expand_power_base

nfloat
^^^^^^
.. autofunction:: nfloat

evalf
-----
.. module:: diofant.core.evalf

.. autoclass:: EvalfMixin
   :members:

PrecisionExhausted
^^^^^^^^^^^^^^^^^^
.. autoclass:: PrecisionExhausted
   :members:

N
^
.. autofunction:: N

containers
----------
.. module:: diofant.core.containers

Tuple
^^^^^
.. autoclass:: Tuple
   :members:

Dict
^^^^
.. autoclass:: Dict
   :members:

compatibility
-------------
.. automodule:: diofant.core.compatibility
   :members:

exprtools
---------
.. module:: diofant.core.exprtools

gcd_terms
^^^^^^^^^
.. autofunction:: gcd_terms

factor_terms
^^^^^^^^^^^^
.. autofunction:: factor_terms
