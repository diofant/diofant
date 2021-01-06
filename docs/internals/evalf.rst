.. _evalf-label:

Numerical evaluation
====================

Floating-point numbers
----------------------

Floating-point numbers in Diofant are instances of the class ``Float``. A ``Float``
can be created with a custom precision as second argument:

    >>> Float(0.1)
    0.100000000000000
    >>> Float(0.1, 10)
    0.1000000000
    >>> Float(0.125, 30)
    0.125000000000000000000000000000
    >>> Float(0.1, 30)
    0.100000000000000005551115123126

As the last example shows, some Python floats are only accurate to about 15
digits as inputs, while others (those that have a denominator that is a
power of 2, like 0.125 = 1/8) are exact. To create a ``Float`` from a
high-precision decimal number, it is better to pass a string, ``Rational``,
or ``evalf`` a ``Rational``:

    >>> Float('0.1', 30)
    0.100000000000000000000000000000
    >>> Float(Rational(1, 10), 30)
    0.100000000000000000000000000000
    >>> Rational(1, 10).evalf(30)
    0.100000000000000000000000000000


The precision of a number determines 1) the precision to use when performing
arithmetic with the number, and 2) the number of digits to display when printing
the number. When two numbers with different precision are used together in an
arithmetic operation, the higher of the precisions is used for the result. The
product of 0.1 +/- 0.001 and 3.1415 +/- 0.0001 has an uncertainty of about 0.003
and yet 5 digits of precision are shown.

    >>> Float(0.1, 3)*Float(3.1415, 5)
    0.31417

So the displayed precision should not be used as a model of error propagation or
significance arithmetic; rather, this scheme is employed to ensure stability of
numerical algorithms.

Function :func:`~diofant.core.evalf.N` (or
:meth:`~diofant.core.evalf.EvalfMixin.evalf` method) can be used to
change the precision of existing floating-point numbers:

    >>> N(3.5)
    3.50000000000000
    >>> N(3.5, 5)
    3.5000

However, you can "increase" precision of the
:class:`~diofant.core.numbers.Float` number only with it's class
constructor:

    >>> Float(3.5, 30)
    3.50000000000000000000000000000


Accuracy and error handling
---------------------------

When the input to ``N`` or ``evalf`` is a complicated expression, numerical
error propagation becomes a concern. As an example, consider the 100'th
Fibonacci number and the excellent (but not exact) approximation `\varphi^{100} / \sqrt{5}`
where `\varphi` is the golden ratio. With ordinary floating-point arithmetic,
subtracting these numbers from each other erroneously results in a complete
cancellation:

    >>> a, b = GoldenRatio**1000/sqrt(5), fibonacci(1000)
    >>> float(a)
    4.3466557686937455e+208
    >>> float(b)
    4.3466557686937455e+208
    >>> float(a) - float(b)
    0.0

``N`` and ``evalf`` keep track of errors and automatically increase the
precision used internally in order to obtain a correct result:

    >>> N(fibonacci(100) - GoldenRatio**100/sqrt(5))
    -5.64613129282185e-22


Unfortunately, numerical evaluation cannot tell an expression that is exactly
zero apart from one that is merely very small. The working precision is
therefore capped, by default to around 100 digits. If we try with the 1000'th
Fibonacci number, the following happens:

    >>> N(fibonacci(1000) - GoldenRatio**1000/sqrt(5))
    Traceback (most recent call last):
    ...
    PrecisionExhausted: ...

The exception indicates that ``N`` failed to achieve full accuracy.  To force a
higher working precision, the ``maxn`` keyword argument can be used:

    >>> N(fibonacci(1000) - GoldenRatio**1000/sqrt(5), maxn=500)
    -4.60123853010113e-210


Normally, ``maxn`` can be set very high (thousands of digits), but be aware that
this may cause significant slowdown in extreme cases.

Also, you can set ``strict`` keyword argument to ``False`` to obtain imprecise
answer instead of exception.  For example, if we add a term so that the
Fibonacci approximation becomes exact (the full form of Binet's formula), we
get an expression that is exactly zero, but ``N`` does not know this:

    >>> f = fibonacci(100) - (GoldenRatio**100 - (GoldenRatio-1)**100)/sqrt(5)
    >>> N(f, strict=False)
    0.e-126
    >>> N(f, maxn=1000, strict=False)
    0.e-1336


In situations where such cancellations are known to occur, the ``chop`` options
is useful. This basically replaces very small numbers in the real or
imaginary portions of a number with exact zeros:

    >>> N(f, chop=True)
    0
    >>> N(3 + I*f, chop=True)
    3.00000000000000


In situations where you wish to remove meaningless digits, re-evaluation or
the use of the ``round`` method are useful:

    >>> Float('.1')*Float('.12345')
    0.012297
    >>> ans = _
    >>> N(ans, 1)
    0.01
    >>> ans.round(2)
    0.01


If you are dealing with a numeric expression that contains no floats, it
can be evaluated to arbitrary precision. To round the result relative to
a given decimal, the round method is useful:

    >>> v = 10*pi + cos(1)
    >>> N(v)
    31.9562288417661
    >>> v.round(3)
    31.956


Sums and integrals
------------------

Sums (in particular, infinite series) and integrals can be used like regular
closed-form expressions, and support arbitrary-precision evaluation:

    >>> Sum(1/n**n, (n, 1, oo)).evalf()
    1.29128599706266
    >>> Integral(x**(-x), (x, 0, 1)).evalf()
    1.29128599706266
    >>> Sum(1/n**n, (n, 1, oo)).evalf(50)
    1.2912859970626635404072825905956005414986193682745
    >>> Integral(x**(-x), (x, 0, 1)).evalf(50)
    1.2912859970626635404072825905956005414986193682745
    >>> (Integral(exp(-x**2), (x, -oo, oo)) ** 2).evalf(30)
    3.14159265358979323846264338328


By default, the tanh-sinh quadrature algorithm is used to evaluate integrals.
This algorithm is very efficient and robust for smooth integrands (and even
integrals with endpoint singularities), but may struggle with integrals that
are highly oscillatory or have mid-interval discontinuities. In many cases,
``evalf``/``N`` will correctly estimate the error. With the following integral,
the result is accurate but only good to four digits:

    >>> f = abs(sin(x))
    >>> Integral(abs(sin(x)), (x, 0, 4)).evalf()
    Traceback (most recent call last):
    ...
    PrecisionExhausted: ...


It is better to split this integral into two pieces:

    >>> (Integral(f, (x, 0, pi)) + Integral(f, (x, pi, 4))).evalf()
    2.34635637913639


A similar example is the following oscillatory integral:


    >>> Integral(sin(x)/x**2, (x, 1, oo)).evalf()
    Traceback (most recent call last):
    ...
    PrecisionExhausted: ...


It can be dealt with much more efficiently by telling ``evalf`` or ``N`` to
use an oscillatory quadrature algorithm:

    >>> Integral(sin(x)/x**2, (x, 1, oo)).evalf(quad='osc')
    0.504067061906928
    >>> Integral(sin(x)/x**2, (x, 1, oo)).evalf(20, quad='osc')
    0.50406706190692837199


Oscillatory quadrature requires an integrand containing a factor cos(ax+b) or
sin(ax+b). Note that many other oscillatory integrals can be transformed to
this form with a change of variables:

    >>> init_printing(pretty_print=True, use_unicode=False,
    ...               wrap_line=False, no_global=True)
    >>> intgrl = Integral(sin(1/x), (x, 0, 1)).transform(x, 1/x)
    >>> intgrl
     oo
      /
     |
     |  sin(x)
     |  ------ dx
     |     2
     |    x
     |
    /
    1
    >>> N(intgrl, quad='osc')
    0.504067061906928


Infinite series use direct summation if the series converges quickly enough.
Otherwise, extrapolation methods (generally the Euler-Maclaurin formula but
also Richardson extrapolation) are used to speed up convergence. This allows
high-precision evaluation of slowly convergent series:

    >>> Sum(1/k**2, (k, 1, oo)).evalf(strict=False)
    1.64493406684823
    >>> zeta(2).evalf()
    1.64493406684823
    >>> Sum(1/k-log(1+1/k), (k, 1, oo)).evalf()
    0.577215664901533
    >>> Sum(1/k-log(1+1/k), (k, 1, oo)).evalf(50)
    0.57721566490153286060651209008240243104215933593992
    >>> EulerGamma.evalf(50)
    0.57721566490153286060651209008240243104215933593992


The Euler-Maclaurin formula is also used for finite series, allowing them to
be approximated quickly without evaluating all terms:

    >>> Sum(1/k, (k, 10000000, 20000000)).evalf()
    0.693147255559946


Note that ``evalf`` makes some assumptions that are not always optimal. For
fine-tuned control over numerical summation, it might be worthwhile to manually
use the method ``Sum.euler_maclaurin``.

Special optimizations are used for rational hypergeometric series (where the
term is a product of polynomials, powers, factorials, binomial coefficients and
the like). ``N``/``evalf`` sum series of this type very rapidly to high
precision. For example, this Ramanujan formula for pi can be summed to 10,000
digits in a fraction of a second with a simple command:

    >>> f = factorial
    >>> R = 9801/sqrt(8)/Sum(f(4*n)*(1103+26390*n)/f(n)**4/396**(4*n),
    ...                      (n, 0, oo))
    >>> N(R, 10000, strict=False)
    3.141592653589793238462643383279502884197169399375105820974944592307...

Numerical simplification
------------------------

The function ``nsimplify`` attempts to find a formula that is numerically equal
to the given input. This feature can be used to guess an exact formula for an
approximate floating-point input, or to guess a simpler formula for a
complicated symbolic input. The algorithm used by ``nsimplify`` is capable of
identifying simple fractions, simple algebraic expressions, linear combinations
of given constants, and certain elementary functional transformations of any of
the preceding.

Optionally, ``nsimplify`` can be passed a list of constants to include (e.g. pi)
and a minimum numerical tolerance. Here are some elementary examples:

    >>> nsimplify(0.1)
    1/10
    >>> nsimplify(6.28, [pi], tolerance=0.01)
    2*pi
    >>> nsimplify(pi, tolerance=0.01)
    22/7
    >>> nsimplify(pi, tolerance=0.001)
    355
    ---
    113
    >>> nsimplify(0.33333, tolerance=1e-4)
    1/3
    >>> nsimplify(2.0**(1/3.), tolerance=0.001)
    635
    ---
    504
    >>> nsimplify(2.0**(1/3.), tolerance=0.001, full=True)
    3 ___
    \/ 2


Here are several more advanced examples:

    >>> nsimplify(Float('0.130198866629986772369127970337', 30), [pi, E])
        1
    ----------
    5*pi
    ---- + 2*E
     7
    >>> nsimplify(cos(atan('1/3')))
        ____
    3*\/ 10
    --------
       10
    >>> nsimplify(4/(1+sqrt(5)), [GoldenRatio])
    -2 + 2*GoldenRatio
    >>> nsimplify(2 + exp(2*atan('1/4')*I))
    49   8*I
    -- + ---
    17    17
    >>> nsimplify((1/(exp(3*pi*I/5)+1)))
               ___________
              /   ___
    1        /  \/ 5    1
    - - I*  /   ----- + -
    2     \/      10    4
    >>> nsimplify(I**I, [pi])
     -pi
     ----
      2
    E
    >>> nsimplify(Sum(1/n**2, (n, 1, oo)), [pi])
      2
    pi
    ---
     6
    >>> nsimplify(gamma('1/4')*gamma('3/4'), [pi])
      ___
    \/ 2 *pi

uFuncify
--------

While NumPy operations are very efficient for vectorized data they sometimes
incur unnecessary costs when chained together. Consider the following operation

.. code:: python

    x = get_numpy_array(...)
    y = sin(x)/x

The operators ``sin`` and ``/`` call routines that execute tight for loops in
``C``. The resulting computation looks something like this

.. code:: c

    for(int i = 0; i < n; i++)
    {
        temp[i] = sin(x[i]);
    }
    for(int i = i; i < n; i++)
    {
        y[i] = temp[i] / x[i];
    }

This is slightly sub-optimal because

1.  We allocate an extra ``temp`` array
2.  We walk over ``x`` memory twice when once would have been sufficient

A better solution would fuse both element-wise operations into a single for loop

.. code:: c

    for(int i = i; i < n; i++)
    {
        y[i] = sin(x[i]) / x[i];
    }

Statically compiled projects like NumPy are unable to take advantage of such
optimizations. Fortunately, Diofant is able to generate efficient low-level C
or Fortran code. It can then depend on projects like ``Cython`` or ``f2py`` to
compile and reconnect that code back up to Python. Fortunately this process is
well automated and a Diofant user wishing to make use of this code generation
should call the ``ufuncify`` function

    >>> expr = sin(x)/x

    >>> from diofant.utilities.autowrap import ufuncify
    >>> f = ufuncify((x,), expr)

This function ``f`` consumes and returns a NumPy array. Generally ``ufuncify``
performs at least as well as ``lambdify``. If the expression is complicated
then ``ufuncify`` often significantly outperforms the NumPy backed solution.
Jensen has a good `blog post <https://ojensen.wordpress.com/2010/08/10/fast-ufunc-ish-hydrogen-solutions/>`_
on this topic.
