=========
 Gotchas
=========

..
    >>> init_printing(pretty_print=True, use_unicode=True)

Lets recall again, that Diofant is nothing more than a Python library,
like :mod:`numpy` or even the Python standard library module
:mod:`sys`.  What this means is that Diofant does not add anything to
the Python language.  Limitations that are inherent in the language
are also inherent in Diofant.

In this section we are trying to collect some things that could
surprise newcomers.

Numbers
=======

To begin with, it should be clear for you, that if you type a numeric
literal --- it will create a Python number of type :class:`int` or
:class:`float`.

Diofant uses its own classes for numbers, for example
:class:`~diofant.core.numbers.Integer` instead of :class:`int`.  In
most cases, Python numeric types will be correctly coersed to Diofant
numbers during expression construction.

    >>> 3 + x**2
     2
    x  + 3
    >>> type(_ - x**2)
    <class 'diofant.core.numbers.Integer'>

But if you use some arithmetic operators between two numerical
literals, Python will evaluate such expression before Diofant has a
chance to get to them.

    >>> x**(3/2)
     1.5
    x

.. tip::

   While working in the IPython console, you could use
   :class:`~diofant.interactive.session.IntegerDivisionWrapper` AST
   transformer to wrap all integer divisions with
   :class:`~diofant.core.numbers.Rational` automatically.

The universal solution is using correct Diofant numeric class to
construct numbers explicitly.  For example,
:class:`~diofant.core.numbers.Rational` in the above example

    >>> x**Rational(3, 2)
     3/2
    x

Equality
========

You may think that ``==``, which is used for equality testing in
Python, is used for Diofant to test mathematical equality.  This is
not quite correct either.  Let us see what happens when we use ``==``.

    >>> (x + 1)**2 == x**2 + 2*x + 1
    False

But, `(x + 1)^2` *does* equal `x^2 + 2x + 1`. What is going on here?

In Diofant, ``==`` represents structural equality testing and `(x +
1)^2` and `x^2 + 2x + 1` are not the same in this sense.  One is the
power and the other is the addition of three terms.

There is a separate class, called
:class:`~diofant.core.relational.Eq`, which can be used to create a
symbolic equation

    >>> Eq((x + 1)**2 - x**2, 2*x + 1)
       2          2
    - x  + (x + 1)  = 2⋅x + 1

It is not always return a :class:`bool` object, like ``==``, but you
may use some simplification methods to prove (or disprove) equation.

    >>> expand(_)
    true

Naming of Functions
===================

Diofant uses different names for some mathematical functions than most
computer algebra systems.  In particular, the inverse trigonometric
functions use the python names of
:func:`~diofant.functions.elementary.trigonometric.asin`,
:func:`~diofant.functions.elementary.trigonometric.acos` and so on
instead of ``arcsin`` and ``arccos``.
