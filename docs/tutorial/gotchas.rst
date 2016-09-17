=====================
 Gotchas and Pitfalls
=====================

To begin, we should make something about Diofant clear.  Diofant is nothing more
than a Python library, like ``NumPy``, ``Django``, or even modules in the
Python standard library ``sys`` or ``re``.  What this means is that Diofant does
not add anything to the Python language.  Limitations that are inherent in the
Python language are also inherent in Diofant.  It also means that Diofant tries to
use Python idioms whenever possible, making programming with Diofant easy for
those already familiar with programming with Python.

For example, implicit multiplication (like ``3x`` or ``3 x``) is not
allowed in Python, and thus not allowed in Diofant: to multiply ``3``
and ``x``, you must type ``3*x`` with the ``*``.  Also, to
raise something to a power, use ``**``, not ``^`` (logical exclusive
or in Python) as many computer algebra systems use.  Parentheses
``()`` change operator precedence as you would normally expect.

.. _tutorial-gotchas-symbols:

Variables and Symbols
=====================

One consequence of this fact is that Diofant can be used in any environment
where Python is available.  We just import it, like we would any other
library:

    >>> from diofant import *

This imports all the functions and classes from Diofant into our interactive
Python session.  Now, suppose we start to do a computation.

    >>> x + 1
    Traceback (most recent call last):
    ...
    NameError: name 'x' is not defined

Oops! What happened here?  We tried to use the variable ``x``, but it tells us
that ``x`` is not defined.  In Python, variables have no meaning until they
are defined.  Diofant is no different.  Unlike many symbolic manipulation
systems you may have used, in Diofant, variables are not defined automatically.
To define variables, we must use ``symbols``.

    >>> x = symbols('x')
    >>> x + 1
    x + 1

``symbols`` takes a string of variable names separated by spaces or commas,
and creates Symbols out of them.  We can then assign these to variable names.
Later, we will investigate some convenient ways we can work around this issue.
For now, let us just define the most common variable names, ``x``, ``y``, and
``z``, for use through the rest of this section

    >>> x, y, z = symbols('x y z')

As a final note, we note that the name of a Symbol and the name of the
variable it is assigned to need not have anything to do with one another.

    >>> a, b = symbols('b a')
    >>> a
    b
    >>> b
    a

Here we have done the very confusing thing of assigning a Symbol with the name
``a`` to the variable ``b``, and a Symbol of the name ``b`` to the variable
``a``.  Now the Python variable named ``a`` points to the Diofant Symbol named
``b``, and visa versa.  How confusing.  We could have also done something like

    >>> crazy = symbols('unrelated')
    >>> crazy + 1
    unrelated + 1

This also shows that Symbols can have names longer than one character if we
want.

Usually, the best practice is to assign Symbols to Python variables of the
same name, although there are exceptions:  Symbol names can contain characters
that are not allowed in Python variable names, or may just want to avoid
typing long names by assigning Symbols with long names to single letter Python
variables.

To avoid confusion, throughout this tutorial, Symbol names and Python variable
names will always coincide.  Furthermore, the word "Symbol" will refer to a
Diofant Symbol and the word "variable" will refer to a Python variable.

Finally, let us be sure we understand the difference between Diofant Symbols and
Python variables.  Consider the following::

  x = symbols('x')
  expr = x + 1
  x = 2
  print(expr)

What do you think the output of this code will be?  If you thought ``3``,
you're wrong.  Let's see what really happens

    >>> x = symbols('x')
    >>> expr = x + 1
    >>> x = 2
    >>> print(expr)
    x + 1

Changing ``x`` to ``2`` had no effect on ``expr``.  This is because ``x = 2``
changes the Python variable ``x`` to ``2``, but has no effect on the Diofant
Symbol ``x``, which was what we used in creating ``expr``.  When we created
``expr``, the Python variable ``x`` was a Symbol.  After we created, it, we
changed the Python variable ``x`` to 2.  But ``expr`` remains the same.  This
behavior is not unique to Diofant.  All Python programs work this way: if a
variable is changed, expressions that were already created with that variable
do not change automatically.  For example

    >>> x = 'abc'
    >>> expr = x + 'def'
    >>> expr
    'abcdef'
    >>> x = 'ABC'
    >>> expr
    'abcdef'


.. tip::

   To change the value of a Symbol in an expression, use ``subs``

     >>> x = symbols('x')
     >>> expr = x + 1
     >>> expr.subs(x, 2)
     3

In this example, if we want to know what ``expr`` is with the new value of
``x``, we need to reevaluate the code that created ``expr``, namely, ``expr =
x + 1``.  This can be complicated if several lines created ``expr``.  One
advantage of using a symbolic computation system like Diofant is that we can
build a symbolic representation for ``expr``, and then substitute ``x`` with
values.  The correct way to do this in Diofant is to use ``subs``, which will be
discussed in more detail later.

    >>> x = symbols('x')
    >>> expr = x + 1
    >>> expr.subs(x, 2)
    3

You can also import common symbol names from ``diofant.abc`` module.

    >>> from diofant.abc import w
    >>> w
    w
    >>> import diofant
    >>> dir(diofant.abc)
    ['A', 'B', ..., 'zeta']

If you want control over the assumptions of the variables, use
:func:`~diofant.core.symbol.Symbol` and :func:`~diofant.core.symbol.symbols`.

Lastly, it is recommended that you not use :class:`I <diofant.core.numbers.ImaginaryUnit>`,
:class:`E <diofant.core.numbers.Exp1>`, :class:`~diofant.core.singleton.S`,
:func:`~diofant.core.evalf.N`, or :class:`O <diofant.series.order.Order>`,
for variable or symbol names, as those
are used for the imaginary unit (:math:`i`), the base of the natural
logarithm (:math:`e`), the :func:`~diofant.core.sympify.sympify` function (see :ref:`Symbolic
Expressions<symbolic-expressions>` below), numeric evaluation (:func:`~diofant.core.evalf.N`
is equivalent to :ref:`evalf()<evalf-label>` ),
the `big O <http://en.wikipedia.org/wiki/Big_O_notation>`_ order symbol
(as in :math:`O(n\log{n})`).  You can use the
mnemonic ``QCOSINE`` to remember what Symbols are defined by default in Diofant.
Or better yet, always use lowercase letters for Symbol names.  Python will
not prevent you from overriding default Diofant names or functions, so be
careful.

    >>> cos(pi)  # cos and pi are a built-in diofant names.
    -1
    >>> pi = 3   # Notice that there is no warning for overriding pi.
    >>> cos(pi)
    cos(3)
    >>> def cos(x):  # No warning for overriding built-in functions either.
    ...     return 5*x
    ...
    >>> cos(pi)
    15
    >>> from diofant import cos  # reimport to restore normal behavior

To get a full list of all default names in Diofant do:

    >>> import diofant
    >>> dir(diofant)
    ['AbelianGroup', 'Abs', 'Add', ..., 'zoo']

If you have `IPython <http://ipython.org/>`_ installed and
use :command:`ipython`, instead of bare Python shell, you can also press
the TAB key to get a list of all built-in names and to autocomplete.

.. _tutorial_gotchas_equals:

Equals signs
============

Another very important consequence of the fact that Diofant does not extend
Python syntax is that ``=`` does not represent equality in Diofant.  Rather it
is Python variable assignment.  This is hard-coded into the Python language,
and Diofant makes no attempts to change that.

You may think, however, that ``==``, which is used for equality testing in
Python, is used for Diofant as equality.  This is not quite correct either.  Let
us see what happens when we use ``==``.

    >>> x + 1 == 4
    False

Instead of treating ``x + 1 == 4`` symbolically, we just got ``False``.  In
Diofant, ``==`` represents exact structural equality testing.  This means that
``a == b`` means that we are *asking* if `a = b`.  We always get a ``bool`` as
the result of ``==``.  There is a separate object, called ``Eq``, which can be
used to create symbolic equalities

    >>> Eq(x + 1, 4)
    Eq(x + 1, 4)

There is one additional caveat about ``==`` as well.  Suppose we want to know
if `(x + 1)^2 = x^2 + 2x + 1`.  We might try something like this:

    >>> (x + 1)**2 == x**2 + 2*x + 1
    False

We got ``False`` again. However, `(x + 1)^2` *does* equal `x^2 + 2x + 1`. What
is going on here?  Did we find a bug in Diofant, or is it just not powerful
enough to recognize this basic algebraic fact?

Recall from above that ``==`` represents *exact* structural equality testing.
"Exact" here means that two expressions will compare equal with ``==`` only if
they are exactly equal structurally.  Here, `(x + 1)^2` and `x^2 + 2x + 1` are
not the same symbolically. One is the power of an addition of two terms, and
the other is the addition of three terms.

It turns out that when using Diofant as a library, having ``==`` test for exact
symbolic equality is far more useful than having it represent symbolic
equality, or having it test for mathematical equality.  However, as a new
user, you will probably care more about the latter two.  We have already seen
an alternative to representing equalities symbolically, ``Eq``.  To test if
two things are equal, it is best to recall the basic fact that if `a = b`,
then `a - b = 0`.  Thus, the best way to check if `a = b` is to take `a - b`
and simplify it, and see if it goes to 0.  We will learn :ref:`later
<tutorial-simplify>` that the function to do this is called ``simplify``. This
method is not infallible---in fact, it can be `theoretically proven
<http://en.wikipedia.org/wiki/Richardson%27s_theorem>`_ that it is impossible
to determine if two symbolic expressions are identically equal in
general---but for most common expressions, it works quite well.

    >>> a = (x + 1)**2
    >>> b = x**2 + 2*x + 1
    >>> simplify(a - b)
    0
    >>> c = x**2 - 2*x + 1
    >>> simplify(a - c)
    4*x

There is also a method called ``equals`` that tests if two expressions are
equal by evaluating them numerically at random points.

    >>> a = cos(x)**2 - sin(x)**2
    >>> b = cos(2*x)
    >>> a.equals(b)
    True


.. _symbolic-expressions:

Symbolic Expressions
====================

.. _python-vs-diofant-numbers:

Python numbers vs. Diofant Numbers
----------------------------------

Diofant uses its own classes for integers, rational numbers, and floating
point numbers instead of the default Python `int` and `float`
types because it allows for more control.  But you have to be careful.
If you type an expression that just has numbers in it, it will default
to a Python expression.  Use the :func:`diofant.core.sympify.sympify` function, or just
:func:`S <diofant.core.sympify.sympify>`, to ensure that something is a Diofant expression.

    >>> 6.2  # Python float. Notice the floating point accuracy problems.
    6.2
    >>> type(6.2)
    <class 'float'>
    >>> Float(6.2)  # Diofant Float has no such problems because of arbitrary precision.
    6.20000000000000

If you include numbers in a Diofant expression, they will be sympified
automatically, but there is one gotcha you should be aware of.  If you
do ``<number>/<number>`` inside of a Diofant expression, Python will
evaluate the two numbers before Diofant has a chance to get
to them.  The solution is to :func:`~diofant.core.sympify.sympify` one of the
numbers, or use :class:`~diofant.core.numbers.Rational`.

    >>> x**(1/2)
    x**0.5
    >>> x**Rational(1, 2)  # use the Rational class
    sqrt(x)

With a power of ``1/2`` you can also use ``sqrt`` shorthand:

    >>> sqrt(x) == x**Rational(1, 2)
    True

If the two integers are not directly separated by a division sign then
you don't have to worry about this problem:

    >>> x**(2*x/3)
    x**(2*x/3)

.. note::

    A common mistake is copying an expression that is printed and reusing
    it.  If the expression has a :class:`~diofant.core.numbers.Rational`
    (i.e., ``<number>/<number>``) in it, you will not get the same result,
    obtaining the Python result for the division rather than a Diofant
    Rational.

    >>> x = Symbol('x')
    >>> solve(7*x - 22, x)
    [22/7]
    >>> 22/7  # After copy and paste we get a float
    3.142857142857143
    >>> # One solution is to just assign the expression to a variable
    >>> # if we need to use it again.
    >>> a = solve(7*x - 22, x)
    >>> a
    [22/7]

    The other solution is using the Rational class:

    >>> Rational(22, 7)
    22/7

    >>> 1/2   # With division imported it evaluates to a python float
    0.5
    >>> 1//2  # You can still achieve integer division with //
    0

    But be careful: you will now receive floats where you might have desired
    a Rational:

    >>> x**(1/2)
    x**0.5

:class:`~diofant.core.numbers.Rational` only works for number/number and is only meant for
rational numbers.  If you want a fraction with symbols or expressions in
it, just use ``/``.  If you do number/expression or expression/number,
then the number will automatically be converted into a Diofant Number.
You only need to be careful with number/number.

    >>> Rational(2, x)
    Traceback (most recent call last):
    ...
    TypeError: invalid input: x
    >>> 2/x
    2/x

Evaluating Expressions with Floats and Rationals
------------------------------------------------

Diofant keeps track of the precision of ``Float`` objects. The default precision is
15 digits. When an expression involving a ``Float`` is evaluated, the result
will be expressed to 15 digits of precision but those digits (depending
on the numbers involved with the calculation) may not all be significant.

The first issue to keep in mind is how the ``Float`` is created: it is created
with a value and a precision. The precision indicates how precise of a value
to use when that ``Float`` (or an expression it appears in) is evaluated.

The values can be given as strings, integers, floats, or rationals.

    - strings and integers are interpreted as exact

    >>> Float(100)
    100.000000000000
    >>> Float('100', 5)
    100.00

    - to have the precision match the number of digits, the null string
      can be used for the precision

    >>> Float(100, '')
    100.
    >>> Float('12.34')
    12.3400000000000
    >>> Float('12.34', '')
    12.34

    >>> s, r = [Float(j, 3) for j in ('0.25', Rational(1, 7))]
    >>> for f in [s, r]:
    ...     print(f)
    0.250
    0.143

Next, notice that each of those values looks correct to 3 digits. But if we try
to evaluate them to 20 digits, a difference will become apparent:

    The 0.25 (with precision of 3) represents a number that has a non-repeating
    binary decimal; 1/7 is repeating in binary and decimal -- it cannot be
    represented accurately too far past those first 3 digits (the correct
    decimal is a repeating 142857):

    >>> s.n(20)
    0.25000000000000000000
    >>> r.n(20)
    0.14285278320312500000

    It is important to realize that although a Float is being displayed in
    decimal at aritrary precision, it is actually stored in binary. Once the
    Float is created, its binary information is set at the given precision.
    The accuracy of that value cannot be subsequently changed; so 1/7, at a
    precision of 3 digits, can be padded with binary zeros, but these will
    not make it a more accurate value of 1/7.

If inexact, low-precision numbers are involved in a calculation with
with higher precision values, the evalf engine will increase the precision
of the low precision values and inexact results will be obtained. This is
feature of calculations with limited precision:

    >>> Float('0.1', 10) + Float('0.1', 3)
    0.2000061035

Although the ``evalf`` engine tried to maintain 10 digits of precision (since
that was the highest precision represented) the 3-digit precision used
limits the accuracy to about 4 digits -- not all the digits you see
are significant. evalf doesn't try to keep track of the number of
significant digits.

That very simple expression involving the addition of two numbers with
different precisions will hopefully be instructive in helping you
understand why more complicated expressions (like trig expressions that
may not be simplified) will not evaluate to an exact zero even though,
with the right simplification, they should be zero. Consider this
unsimplified trig identity, multiplied by a big number:

    >>> big = 12345678901234567890
    >>> big_trig_identity = big*cos(x)**2 + big*sin(x)**2 - big*1
    >>> abs(big_trig_identity.subs(x, .1).n(2)) > 1000
    true

When the `\cos` and `\sin` terms were evaluated to 15 digits of precision and
multiplied by the big number, they gave a large number that was only
precise to 15 digits (approximately) and when the 20 digit big number
was subtracted the result was not zero.

There are three things that will help you obtain more precise numerical
values for expressions:

    1) Pass the desired substitutions with the call to evaluate. By doing
    the subs first, the ``Float`` values can not be updated as necessary. By
    passing the desired substitutions with the call to evalf the ability
    to re-evaluate as necessary is gained and the results are impressively
    better:

    >>> big_trig_identity.n(2, {x: 0.1})
    -0.e-91

    2) Use Rationals, not Floats. During the evaluation process, the
    Rational can be computed to an arbitrary precision while the Float,
    once created -- at a default of 15 digits -- cannot. Compare the
    value of ``-1.4e+3`` above with the nearly zero value obtained when
    replacing x with a Rational representing 1/10 -- before the call
    to evaluate:

    >>> big_trig_identity.subs(x, Rational(1, 10)).n(2)
    0.e-91

    3) Try to simplify the expression. In this case, Diofant will recognize
    the trig identity and simplify it to zero so you don't even have to
    evaluate it numerically:

    >>> big_trig_identity.simplify()
    0


.. _Immutability-of-Expressions:

Immutability of Expressions
---------------------------

Expressions in Diofant are immutable, and cannot be modified by an in-place
operation.  This means that a function will always return an object, and the
original expression will not be modified. The following example snippet
demonstrates how this works::

	def main():
	    var('x y a b')
	    expr = 3*x + 4*y
	    print('original =', expr)
	    expr_modified = expr.subs({x: a, y: b})
	    print('modified =', expr_modified)

	if __name__ == "__main__":
	    main()

The output shows that the :func:`~diofant.core.basic.Basic.subs` function has replaced variable
``x`` with variable ``a``, and variable ``y`` with variable ``b``::

	original = 3*x + 4*y
	modified = 3*a + 4*b

The :func:`~diofant.core.basic.Basic.subs` function does not modify the original expression `expr``.
Rather, a modified copy of the expression is returned. This returned object
is stored in the variable ``expr_modified``. Note that unlike C/C++ and
other high-level languages, Python does not require you to declare a variable
before it is used.

Inverse Trig Functions
----------------------

Diofant uses different names for some functions than most computer algebra
systems.  In particular, the inverse trig functions use the python names
of :func:`~diofant.functions.elementary.trigonometric.asin`,
:func:`~diofant.functions.elementary.trigonometric.acos` and
so on instead of the usual ``arcsin``
and ``arccos``.  Use the methods described in the section
:ref:`Variables and Symbols <tutorial-gotchas-symbols>`
above to see the names of all Diofant functions.
