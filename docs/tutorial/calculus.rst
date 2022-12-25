==========
 Calculus
==========

..
    >>> init_printing(pretty_print=True, use_unicode=True)

This section covers how to do basic calculus tasks such as derivatives,
integrals, limits, and series expansions in Diofant.

Derivatives
===========

To take derivatives, use the :func:`~diofant.core.function.diff`
function.

    >>> diff(cos(x))
    -sin(x)
    >>> diff(exp(x**2), x)
       ⎛ 2⎞
       ⎝x ⎠
    2⋅ℯ    ⋅x

:func:`~diofant.core.function.diff` can take multiple derivatives at
once.  To take multiple derivatives, pass the variable as many times
as you wish to differentiate, or pass a tuple (variable and order).
For example, both of the following find the third derivative of `x^4`.

    >>> diff(x**4, x, x, x)
    24⋅x
    >>> diff(x**4, (x, 3))
    24⋅x

You can also take derivatives with respect to many variables at once.  Just
pass each derivative in order, using the same syntax as for single variable
derivatives.  For example, each of the following will compute
`\frac{\partial^7}{\partial x\partial y^2\partial z^4} e^{x y z}`.

    >>> expr = exp(x*y*z)
    >>> diff(expr, x, y, y, z, z, z, z)
     x⋅y⋅z  3  2 ⎛ 3  3  3       2  2  2                ⎞
    ℯ     ⋅x ⋅y ⋅⎝x ⋅y ⋅z  + 14⋅x ⋅y ⋅z  + 52⋅x⋅y⋅z + 48⎠
    >>> diff(expr, x, (y, 2), (z, 4))
     x⋅y⋅z  3  2 ⎛ 3  3  3       2  2  2                ⎞
    ℯ     ⋅x ⋅y ⋅⎝x ⋅y ⋅z  + 14⋅x ⋅y ⋅z  + 52⋅x⋅y⋅z + 48⎠
    >>> diff(expr, x, y, y, (z, 4))
     x⋅y⋅z  3  2 ⎛ 3  3  3       2  2  2                ⎞
    ℯ     ⋅x ⋅y ⋅⎝x ⋅y ⋅z  + 14⋅x ⋅y ⋅z  + 52⋅x⋅y⋅z + 48⎠

:func:`~diofant.core.function.diff` can also be called as a method
:meth:`~diofant.core.expr.Expr.diff`.  The two ways of calling
:func:`~diofant.core.function.diff` are exactly the same, and are
provided only for convenience.

    >>> expr.diff(x, y, y, (z, 4))
     x⋅y⋅z  3  2 ⎛ 3  3  3       2  2  2                ⎞
    ℯ     ⋅x ⋅y ⋅⎝x ⋅y ⋅z  + 14⋅x ⋅y ⋅z  + 52⋅x⋅y⋅z + 48⎠

To create an unevaluated derivative, use the
:class:`~diofant.core.function.Derivative` class.  It has the same
syntax as :func:`~diofant.core.function.diff`.

    >>> Derivative(expr, x, y, y, (z, 4))
         7
        ∂     ⎛ x⋅y⋅z⎞
    ──────────⎝ℯ     ⎠
      4   2
    ∂z  ∂y  ∂x

Such unevaluated objects also used when Diofant does not know how to compute
the derivative of an expression.

To evaluate an unevaluated derivative, use the
:meth:`~diofant.core.basic.Basic.doit` method.

    >>> _.doit()
     x⋅y⋅z  3  2 ⎛ 3  3  3       2  2  2                ⎞
    ℯ     ⋅x ⋅y ⋅⎝x ⋅y ⋅z  + 14⋅x ⋅y ⋅z  + 52⋅x⋅y⋅z + 48⎠

Integrals
=========

To compute an integral, use the :func:`~diofant.integrals.integrals.integrate`
function.  There are two kinds of integrals, definite and indefinite.  To
compute an indefinite integral, do

    >>> integrate(cos(x))
    sin(x)

.. note::

    For indefinite integrals, Diofant does not include the constant of
    integration.

For example, to compute a definite integral

.. math::

   \int_0^\infty e^{-x}\,dx,

we would do

    >>> integrate(exp(-x), (x, 0, oo))
    1

.. tip::

    `\infty` in Diofant is ``oo`` (that's the lowercase letter "oh" twice).

As with indefinite integrals, you can pass multiple limit tuples to perform a
multiple integral.  For example, to compute

.. math::

   \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} e^{- x^{2} - y^{2}}\, dx\, dy,

do

    >>> integrate(exp(-x**2 - y**2), (x, -oo, oo), (y, -oo, oo))
    π

If :func:`~diofant.integrals.integrals.integrate` is unable to compute an
integral, it returns an unevaluated
:class:`~diofant.integrals.integrals.Integral` object.

    >>> integrate(x**x)
    ⌠
    ⎮  x
    ⎮ x  dx
    ⌡
    >>> print(_)
    Integral(x**x, x)

As with :class:`~diofant.core.function.Derivative`, you can create an
unevaluated integral directly.  To later evaluate this integral, call
:meth:`~diofant.integrals.integrals.Integral.doit`.

    >>> Integral(log(x)**2)
    ⌠
    ⎮    2
    ⎮ log (x) dx
    ⌡
    >>> _.doit()
             2
    x⋅log (x) - 2⋅x⋅log(x) + 2⋅x

:func:`~diofant.integrals.integrals.integrate` uses powerful algorithms that
are always improving to compute both definite and indefinite integrals,
including a partial implementation of the `Risch algorithm
<https://en.wikipedia.org/wiki/Risch_algorithm>`_

    >>> Integral((x**4 + x**2*exp(x) - x**2 - 2*x*exp(x) - 2*x -
    ...           exp(x))*exp(x)/((x - 1)**2*(x + 1)**2*(exp(x) + 1)))
    ⌠
    ⎮  x ⎛ x  2      x      x    4    2      ⎞
    ⎮ ℯ ⋅⎝ℯ ⋅x  - 2⋅ℯ ⋅x - ℯ  + x  - x  - 2⋅x⎠
    ⎮ ──────────────────────────────────────── dx
    ⎮        ⎛ x    ⎞        2        2
    ⎮        ⎝ℯ  + 1⎠⋅(x - 1) ⋅(x + 1)
    ⌡
    >>> _.doit()
        x
      ℯ         ⎛ x    ⎞
    ────── + log⎝ℯ  + 1⎠
     2
    x  - 1

and an algorithm using `Meijer G-functions
<https://en.wikipedia.org/wiki/Meijer_g-function>`_ that is useful for computing
integrals in terms of special functions, especially definite integrals

    >>> Integral(sin(x**2))
    ⌠
    ⎮    ⎛ 2⎞
    ⎮ sin⎝x ⎠ dx
    ⌡
    >>> _.doit()
                          ⎛  ___  ⎞
        ___   ___         ⎜╲╱ 2 ⋅x⎟
    3⋅╲╱ 2 ⋅╲╱ π ⋅fresnels⎜───────⎟⋅Γ(3/4)
                          ⎜   ___ ⎟
                          ⎝ ╲╱ π  ⎠
    ──────────────────────────────────────
                   8⋅Γ(7/4)


    >>> Integral(x**y*exp(-x), (x, 0, oo))
    ∞
    ⌠
    ⎮  -x  y
    ⎮ ℯ  ⋅x  dx
    ⌡
    0
    >>> _.doit()
    ⎧ Γ(y + 1)    for -re(y) < 1
    ⎪
    ⎪∞
    ⎪⌠
    ⎨⎮  -x  y
    ⎪⎮ ℯ  ⋅x  dx    otherwise
    ⎪⌡
    ⎪0
    ⎩

This last example returned a
:class:`~diofant.functions.elementary.piecewise.Piecewise` expression because
the integral does not converge unless `\Re(y) > 1.`

Sums and Products
=================

Much like integrals, there are
:func:`~diofant.concrete.summations.summation` and
:func:`~diofant.concrete.products.product` to compute sums and
products respectively.

    >>> summation(2**x, (x, 0, y - 1))
     y
    2  - 1
    >>> product(z, (x, 1, y))
     y
    z

Unevaluated sums and products are represented by
:class:`~diofant.concrete.summations.Sum` and
:class:`~diofant.concrete.products.Product` classes.

    >>> Sum(1, (x, 1, z))
     z
    ___
    ╲
     ╲   1
     ╱
    ╱
    ‾‾‾
    x = 1
    >>> _.doit()
    z

Limits
======

Diofant can compute symbolic limits with the
:func:`~diofant.calculus.limits.limit` function.  To compute a directional limit

.. math::

   \lim_{x\to 0^+}\frac{\sin x}{x},

do

    >>> limit(sin(x)/x, x, 0)
    1

:func:`~diofant.calculus.limits.limit` should be used instead of
:meth:`~diofant.core.basic.Basic.subs` whenever the point of evaluation is a
singularity.  Even though Diofant has objects to represent `\infty`, using them
for evaluation is not reliable because they do not keep track of things like
rate of growth.  Also, things like `\infty - \infty` and
`\frac{\infty}{\infty}` return `\mathrm{nan}` (not-a-number).  For example

    >>> expr = x**2/exp(x)
    >>> expr.subs({x: oo})
    nan
    >>> limit(expr, x, oo)
    0

Like :class:`~diofant.core.function.Derivative` and
:class:`~diofant.integrals.integrals.Integral`,
:func:`~diofant.calculus.limits.limit` has an unevaluated counterpart,
:class:`~diofant.calculus.limits.Limit`.  To evaluate it, use
:meth:`~diofant.calculus.limits.Limit.doit`.

    >>> Limit((cos(x) - 1)/x, x, 0)
         cos(x) - 1
     lim ──────────
    x─→0⁺    x
    >>> _.doit()
    0

To change side, pass ``'-'`` as a third argument to
:func:`~diofant.calculus.limits.limit`.  For example, to compute

.. math::

   \lim_{x\to 0^-}\frac{1}{x},

do

    >>> limit(1/x, x, 0, dir=1)
    -∞

You can also evaluate bidirectional limit

    >>> limit(sin(x)/x, x, 0, dir=Reals)
    1
    >>> limit(1/x, x, 0, dir=Reals)
    Traceback (most recent call last):
    ...
    PoleError: left and right limits for expression 1/x at point x=0 seems to be not equal

Series Expansion
================

Diofant can compute asymptotic series expansions of functions around a point.

    >>> exp(sin(x)).series(x, 0, 4)
             2
            x     ⎛ 4⎞
    1 + x + ── + O⎝x ⎠
            2

The `O\left (x^4\right )` term, an instance of :class:`~diofant.calculus.order.O`
at the end represents the Landau order term at `x=0` (not to be confused with
big O notation used in computer science, which generally represents the Landau
order term at `x=\infty`).  Order terms can be created and manipulated outside
of :meth:`~diofant.core.expr.Expr.series`.

    >>> x + x**3 + x**6 + O(x**4)
         3    ⎛ 4⎞
    x + x  + O⎝x ⎠
    >>> x*O(1)
    O(x)

If you do not want the order term, use the
:meth:`~diofant.core.expr.Expr.removeO` method.

    >>> exp(x).series(x, 0, 3).removeO()
     2
    x
    ── + x + 1
    2

The :class:`~diofant.calculus.order.O` notation supports arbitrary limit points:

    >>> exp(x - 1).series(x, x0=1)
           2          3          4          5
    (x - 1)    (x - 1)    (x - 1)    (x - 1)         ⎛       6       ⎞
    ──────── + ──────── + ──────── + ──────── + x + O⎝(x - 1) ; x → 1⎠
       2          6          24        120
