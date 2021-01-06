================
 Simplification
================

..
    >>> init_printing(pretty_print=True, use_unicode=True)

The generic way to do *nontrivial* simplifications of expressions is
calling :func:`~diofant.simplify.simplify.simplify` function.

    >>> simplify(sin(x)**2 + cos(x)**2)
    1
    >>> simplify((x**3 + x**2 - x - 1)/(x**2 + 2*x + 1))
    x - 1
    >>> simplify(gamma(x)/gamma(x - 2))
    (x - 2)⋅(x - 1)

There are also more directed simplification functions.  These apply
very specific rules to the input expression and are typically able to
make guarantees about the output.  For instance, the
:func:`~diofant.polys.polytools.factor` function, given a polynomial
with rational coefficients in several variables, is guaranteed to
produce a factorization into irreducible factors.

The :func:`~diofant.simplify.simplify.simplify` function applies
almost all available in Diofant such specific simplification rules in
some heuristics sequence to produce the simplest result.

.. tip::

   The optional ``measure`` keyword argument for
   :func:`~diofant.simplify.simplify.simplify` lets the user specify
   the Python function used to determine how "simple" an expression
   is.  The default is :func:`~diofant.core.function.count_ops`, which
   returns the total number of operations in the expression.

That is why it is usually slow.  But more important pitfail is that
sometimes :func:`~diofant.simplify.simplify.simplify` doesn't
"simplify" how you might expect, if, for example, it miss some
transformation or apply it too early or too late.  Lets look on an
example

    >>> simplify(x**2 + 2*x + 1)
     2
    x  + 2⋅x + 1
    >>> factor(_)
           2
    (x + 1)

Obviously, the factored form is more "simple", as it has less
arithmetic operations.

The function :func:`~diofant.simplify.simplify.simplify` is best when
used interactively, when you just want to whittle down an expression
to a simpler form.  You may then choose to apply specific functions
once you see what :func:`~diofant.simplify.simplify.simplify` returns,
to get a more precise result.  It is also useful when you have no idea
what form an expression will take, and you need a catchall function to
simplify it.

Rational Functions
==================

:func:`~diofant.core.function.expand` is one of the most common
simplification functions in Diofant.  Although it has a lot of scopes,
for now, we will consider its function in expanding polynomial
expressions.

    >>> expand((x + 1)**2)
     2
    x  + 2⋅x + 1
    >>> expand((x + 2)*(x - 3))
     2
    x  - x - 6

Given a polynomial, :func:`~diofant.core.function.expand` will put it
into a canonical form of a sum of monomials with help of more directed
expansion methods, namely
:func:`~diofant.core.function.expand_multinomial` and
:func:`~diofant.core.function.expand_mul`.

:func:`~diofant.core.function.expand` may not sound like a
simplification function.  After all, by its very name, it makes
expressions bigger, not smaller.  Usually this is the case, but often
an expression will become smaller upon calling
:func:`~diofant.core.function.expand` on it due to cancellation.

    >>> expand((x + 1)*(x - 2) - (x - 1)*x)
    -2

Function :func:`~diofant.polys.polytools.factor` takes a multivariate
polynomial with rational coefficients and factors it into irreducible
factors.

    >>> factor(x**3 - x**2 + x - 1)
            ⎛ 2    ⎞
    (x - 1)⋅⎝x  + 1⎠
    >>> factor(x**2*z + 4*x*y*z + 4*y**2*z)
               2
    z⋅(x + 2⋅y)

For polynomials, :func:`~diofant.polys.polytools.factor` is the
opposite of :func:`~diofant.core.function.expand`.

.. note::

   The input to :func:`~diofant.polys.polytools.factor` and
   :func:`~diofant.core.function.expand` need not be polynomials in
   the strict sense.  They will intelligently factor or expand any
   kind of expression (though, for example, the factors may not be
   irreducible if the input is no longer a polynomial over the
   rationals).

       >>> expand((cos(x) + sin(x))**2)
          2                           2
       sin (x) + 2⋅sin(x)⋅cos(x) + cos (x)
       >>> factor(_)
                        2
       (sin(x) + cos(x))

:func:`~diofant.simplify.radsimp.collect` collects common powers of a
term in an expression.

    >>> x*y + x - 3 + 2*x**2 - z*x**2 + x**3
     3    2        2
    x  - x ⋅z + 2⋅x  + x⋅y + x - 3
    >>> collect(_, x)
     3    2
    x  + x ⋅(-z + 2) + x⋅(y + 1) - 3

:func:`~diofant.simplify.radsimp.collect` is particularly useful in
conjunction with the :meth:`~diofant.core.expr.Expr.coeff` method.

    >>> _.coeff(x, 2)
    -z + 2

:func:`~diofant.polys.polytools.cancel` will take any rational
function and put it into the standard canonical form, `p/q`, where `p`
and `q` are expanded polynomials with no common factors.

    >>> 1/x + (3*x/2 - 2)/(x - 4)
    3⋅x
    ─── - 2
     2        1
    ─────── + ─
     x - 4    x
    >>> cancel(_)
       2
    3⋅x  - 2⋅x - 8
    ──────────────
         2
      2⋅x  - 8⋅x

    >>> expr = (x*y**2 - 2*x*y*z + x*z**2 + y**2 - 2*y*z + z**2)/(x**2 - 1)
    >>> expr
       2                2    2            2
    x⋅y  - 2⋅x⋅y⋅z + x⋅z  + y  - 2⋅y⋅z + z
    ───────────────────────────────────────
                      2
                     x  - 1
    >>> cancel(_)
     2            2
    y  - 2⋅y⋅z + z
    ───────────────
         x - 1

.. note::

   Since :func:`~diofant.polys.polytools.factor` will completely
   factorize both the numerator and the denominator of an expression,
   it can also be used to do the same thing:

       >>> factor(expr)
              2
       (y - z)
       ────────
        x - 1

   However, it's less efficient if you are only interested in making
   sure that the expression is in canceled form.

:func:`~diofant.polys.partfrac.apart` performs a `partial fraction
decomposition
<https://en.wikipedia.org/wiki/Partial_fraction_decomposition>`_ on a
rational function.

    >>> (4*x**3 + 21*x**2 + 10*x + 12)/(x**4 + 5*x**3 + 5*x**2 + 4*x)
       3       2
    4⋅x  + 21⋅x  + 10⋅x + 12
    ────────────────────────
      4      3      2
     x  + 5⋅x  + 5⋅x  + 4⋅x
    >>> apart(_)
     2⋅x - 1       1     3
    ────────── - ───── + ─
     2           x + 4   x
    x  + x + 1

Trigonometric Functions
=======================

To simplify expressions using trigonometric identities, use
:func:`~diofant.simplify.trigsimp.trigsimp` function.

    >>> trigsimp(sin(x)**2 + cos(x)**2)
    1
    >>> trigsimp(sin(x)**4 - 2*cos(x)**2*sin(x)**2 + cos(x)**4)
    cos(4⋅x)   1
    ──────── + ─
       2       2
    >>> trigsimp(sin(x)*tan(x)/sec(x))
       2
    sin (x)

It also works with hyperbolic functions.

    >>> trigsimp(cosh(x)**2 + sinh(x)**2)
    cosh(2⋅x)
    >>> trigsimp(sinh(x)/tanh(x))
    cosh(x)

Much like :func:`~diofant.simplify.simplify.simplify` function,
:func:`~diofant.simplify.trigsimp.trigsimp` applies various
trigonometric identities to the input expression, and then uses a
heuristic to return the "best" one.

To expand trigonometric functions, that is, apply the sum or double
angle identities, use :func:`~diofant.core.function.expand_trig`
function.

    >>> expand_trig(sin(x + y))
    sin(x)⋅cos(y) + sin(y)⋅cos(x)
    >>> expand_trig(tan(2*x))
       2⋅tan(x)
    ─────────────
         2
    - tan (x) + 1

Powers and Logarithms
=====================

:func:`~diofant.simplify.powsimp.powdenest` function applies identity
`(x^a)^b = x^{a b}`, from left to right, if assumptions allow.

    >>> a, b = symbols('a b', real=True)
    >>> p = symbols('p', positive=True)
    >>> powdenest((p**a)**b)
     a⋅b
    p

:func:`~diofant.simplify.powsimp.powsimp` function reduces expression
by combining powers with similar bases and exponent.

   >>> powsimp(z**x*z**y)
     x + y
    z

Again, as for :func:`~diofant.simplify.powsimp.powdenest` above, for
the identity `x^a y^a = (x y)^a`, that combine bases, we should be
careful about assumptions.

   >>> q = symbols('q', positive=True)
   >>> powsimp(p**a*q**a)
        a
   (p⋅q)

In general, this identity doesn't hold.  For example, if `x = y = -1`
and `a = 1/2`.

:func:`~diofant.core.function.expand_power_exp` and
:func:`~diofant.core.function.expand_power_base` functions do reverse
of :func:`~diofant.simplify.powsimp.powsimp`.

    >>> expand_power_exp(x**(y + z))
     y  z
    x ⋅x
    >>> expand_power_base((p*q)**a)
     a  a
    p ⋅q

Logarithms have similar issues as powers.  There are two main
identities

1. `\log{(xy)} = \log{(x)} + \log{(y)}`
2. `\log{(x^n)} = n\log{(x)}`

Neither identity is true for arbitrary complex `x` and `y`, due to the
branch cut in the complex plane for the complex logarithm.

To apply above identities from left to right, use
:func:`~diofant.core.function.expand_log`.  As for powers, the
identities will not be applied unless they are valid with given set of
assumptions for symbols.

    >>> expand_log(log(p*q))
    log(p) + log(q)
    >>> expand_log(log(p/q))
    log(p) - log(q)
    >>> expand_log(log(p**2))
    2⋅log(p)
    >>> expand_log(log(p**a))
    a⋅log(p)
    >>> expand_log(log(x*y))
    log(x⋅y)

To apply identities from right to left, i.e. do reverse of
:func:`~diofant.core.function.expand_log`, use
:func:`~diofant.simplify.simplify.logcombine` function.

    >>> logcombine(log(p) + log(q))
    log(p⋅q)
    >>> logcombine(a*log(p))
       ⎛ a⎞
    log⎝p ⎠
    >>> logcombine(a*log(z))
    a⋅log(z)

Special Functions
=================

Diofant implements dozens of :ref:`special functions
<special-functions>`, ranging from functions in combinatorics to
mathematical physics.

To expand special functions in terms of some identities, use
:func:`~diofant.core.function.expand_func`.  For example the `gamma
function <https://en.wikipedia.org/wiki/Gamma_function>`_
:class:`~diofant.functions.special.gamma_functions.gamma` can be
expanded as

    >>> expand_func(gamma(x + 3))
    x⋅(x + 1)⋅(x + 2)⋅Γ(x)

This method also can help if you would like to rewrite the generalized
hypergeometric function
:class:`~diofant.functions.special.hyper.hyper` or the Meijer
G-function :class:`~diofant.functions.special.hyper.meijerg` in terms
of more standard functions.

    >>> expand_func(hyper([1, 1], [2], z))
    -log(-z + 1)
    ─────────────
         z
    >>> meijerg([[1], [1]], [[1], []], -z)
    ╭─╮1, 1 ⎛1  1 │   ⎞
    │╶┐     ⎜     │ -z⎟
    ╰─╯2, 1 ⎝1    │   ⎠
    >>> expand_func(_)
    z ___
    ╲╱ ℯ

Another type of expand rule is expanding complex valued expressions and putting
them into a normal form. For this :func:`~diofant.core.function.expand_complex`
is used.   Note that it will always perform arithmetic expand to obtain the
desired normal form.

    >>> expand_complex(x + I*y)
    ⅈ⋅(re(y) + im(x)) + re(x) - im(y)

The same behavior can be obtained by using
:meth:`~diofant.core.expr.Expr.as_real_imag` method.

    >>> (x + I*y).as_real_imag()
    (re(x) - im(y), re(y) + im(x))

To simplify combinatorial expressions, involving
:class:`~diofant.functions.combinatorial.factorials.factorial`,
:class:`~diofant.functions.combinatorial.factorials.binomial` or
:class:`~diofant.functions.special.gamma_functions.gamma` --- use
:func:`~diofant.simplify.combsimp.combsimp` function.

    >>> combsimp(factorial(x)/factorial(x - 3))
    x⋅(x - 2)⋅(x - 1)
    >>> combsimp(binomial(x + 1, y + 1)/binomial(x, y))
    x + 1
    ─────
    y + 1
    >>> combsimp(gamma(x)*gamma(1 - x))
       π
    ────────
    sin(π⋅x)

CSE
===

Before evaluating a large expression, it is often useful to identify common
subexpressions, collect them and evaluate them at once. This is called common
subexpression elimination (CSE) and implemented in the
:func:`~diofant.simplify.cse_main.cse` function.

    >>> cse(sqrt(sin(x)))
    ⎛    ⎡  ________⎤⎞
    ⎝[], ⎣╲╱ sin(x) ⎦⎠

    >>> cse(sqrt(sin(x) + 5)*sqrt(sin(x) + 4))
    ⎛                ⎡  ________   ________⎤⎞
    ⎝[(x₀, sin(x))], ⎣╲╱ x₀ + 4 ⋅╲╱ x₀ + 5 ⎦⎠

    >>> cse(sqrt(sin(x + 1) + 5 + cos(y))*sqrt(sin(x + 1) + 4 + cos(y)))
    ⎛                             ⎡  ________   ________⎤⎞
    ⎝[(x₀, sin(x + 1) + cos(y))], ⎣╲╱ x₀ + 4 ⋅╲╱ x₀ + 5 ⎦⎠

    >>> cse((x - y)*(z - y) + sqrt((x - y)*(z - y)))
    ⎛                                     ⎡  ____     ⎤⎞
    ⎝[(x₀, -y), (x₁, (x + x₀)⋅(x₀ + z))], ⎣╲╱ x₁  + x₁⎦⎠

Optimizations to be performed before and after common subexpressions
elimination can be passed in the``optimizations`` optional argument.

    >>> cse((x - y)*(z - y) + sqrt((x - y)*(z - y)), optimizations='basic')
    ⎛                          ⎡  ____     ⎤⎞
    ⎝[(x₀, -(x - y)⋅(y - z))], ⎣╲╱ x₀  + x₀⎦⎠

However, these optimizations can be very slow for large expressions. Moreover,
if speed is a concern, one can pass the option ``order='none'``. Order of
terms will then be dependent on hashing algorithm implementation, but speed
will be greatly improved.
