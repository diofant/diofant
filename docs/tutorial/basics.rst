========
 Basics
========

..
   >>> init_printing(pretty_print=True, use_unicode=True)

Here we discuss some of the most basic aspects of expression
manipulation in Diofant.

Assumptions
===========

The assumptions system allows users to declare certain mathematical
properties on symbols, such as being positive, imaginary or integer.

By default, all symbols are complex valued.  This assumption makes it
easier to treat mathematical problems in full generality.

    >>> sqrt(x**2)
       ____
      ╱  2
    ╲╱  x


Yet obviously we can simplify above expression if some additional
mathematical properties on ``x`` are assumed.  This is where
assumptions system come into play.

Assumptions are set on :class:`~diofant.core.symbol.Symbol` objects
when they are created. For instance, we can create a symbol that is
assumed to be positive.

    >>> p = symbols('p', positive=True)

And then, certain simplifications will be possible:

    >>> sqrt(p**2)
    p

The assumptions system additionally has deductive capabilities.  You
might check assumptions on any expression with ``is_assumption``
attributes, like :attr:`~diofant.core.expr.Expr.is_positive`.

    >>> p.is_positive
    True
    >>> (1 + p).is_positive
    True
    >>> (-p).is_positive
    False

.. note::

   ``False`` is returned also if certain assumption doesn't make sense
   for given object.

In a three-valued logic, used by system, ``None`` represents the
"unknown" case.

    >>> (p - 1).is_positive is None
    True

Substitution
============

One of the most common things you might want to do with a mathematical
expression is substitution with :meth:`~diofant.core.basic.Basic.subs`
method.  It replaces all instances of something in an expression with
something else.

    >>> expr = cos(x) + 1
    >>> expr.subs({x: y})
    cos(y) + 1
    >>> expr
    cos(x) + 1

We see that performing substitution leaves original expression
``expr`` unchanged.

.. note::

   Almost all Diofant expressions are immutable.  No function (or
   method) will change them in-place.

Use several method calls to perform a sequence of substitutions in
same variable:

    >>> x**y
     y
    x
    >>> _.subs({y: x**y}).subs({y: x**x})
     ⎛ ⎛ x⎞⎞
     ⎜ ⎝x ⎠⎟
     ⎝x    ⎠
    x

Use flag ``simultaneous`` to do all substitutions at once.

    >>> (x - y).subs({x: y, y: x})
    0
    >>> (x - y).subs({x: y, y: x}, simultaneous=True)
    -x + y

Numerics
========

To evaluate a numerical expression into a floating point number with
arbitrary precision, use :meth:`~diofant.core.evalf.EvalfMixin.evalf`.
By default, 15 digits of precision are used.

    >>> expr = sqrt(8)
    >>> expr.evalf()
    2.82842712474619

But you can change that.  Let's compute the first 70 digits of `\pi`.

    >>> pi.evalf(70)
    3.141592653589793238462643383279502884197169399375105820974944592307816

Complex numbers are supported:

    >>> (1/(pi + I)).evalf()
    0.289025482222236 - 0.0919996683503752⋅ⅈ

If the expression contains symbols or for some other reason cannot be evaluated
numerically, calling :meth:`~diofant.core.evalf.EvalfMixin.evalf` returns the
original expression or a partially evaluated expression.

    >>> (pi*x**2 + x/3).evalf()
                      2
    3.14159265358979⋅x  + 0.333333333333333⋅x

You can also use the standard Python functions :class:`float` and
:class:`complex` to convert symbolic expressions to regular Python numbers:

    >>> float(pi)
    3.141592653589793
    >>> complex(pi + E*I)
    (3.141592653589793+2.718281828459045j)

Sometimes there are roundoff errors smaller than the desired precision
that remain after an expression is evaluated.  Such numbers can be
removed by setting the ``chop`` flag.

    >>> one = cos(1)**2 + sin(1)**2
    >>> (one - 1).evalf(strict=False)
    -0.e-146
    >>> (one - 1).evalf(chop=True)
    0

Discussed above method is not effective enough if you intend to
evaluate an expression at many points, there are better ways,
especially if you only care about machine precision.

Substitution may be used to evaluate an expression for some floating point
number

    >>> expr = sin(x)/x
    >>> expr.subs({x: 0.1})
    0.998334166468282

but this method is slow.

The easiest way to convert an expression to the form that can be numerically
evaluated with libraries like :mod:`numpy` or the standard library :mod:`math`
module --- use the :func:`~diofant.utilities.lambdify.lambdify` function.

    >>> f = lambdify(x, expr, 'math')
    >>> f(0.1)
    0.9983341664682815

Using the :mod:`numpy` library gives the generated function access to powerful
vectorized ufuncs that are backed by compiled C code.

    >>> f = lambdify(x, expr, 'numpy')
    >>> f(range(1, 5))
    [ 0.84147098  0.45464871  0.04704    -0.18920062]
