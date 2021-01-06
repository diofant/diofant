=============
 Polynomials
=============

..
    >>> init_printing(pretty_print=True, use_unicode=True)

We show here some functions, that provide different algorithms dealing
with polynomials in the form of Diofant expression.

Note, that coefficients of a polynomial can be elements of any
commutative ring, this ring is called the ``domain`` of the polynomial
ring and may be specified as a keyword parameter for functions.
Polynomial generators also can be specified via an arbitrary number of
arguments after required arguments of functions.

Division
========

The function :func:`~diofant.polys.polytools.div` provides
division of polynomials with remainder.
That is, for polynomials ``f`` and ``g``, it computes ``q`` and ``r``, such
that `f = g \cdot q + r` and `\deg(r) < q`. For polynomials in one variables
with coefficients in a field, say, the rational numbers, ``q`` and ``r`` are
uniquely defined this way

    >>> f, g = 5*x**2 + 10*x + 3, 2*x + 2

    >>> div(f, g)
    ⎛5⋅x   5    ⎞
    ⎜─── + ─, -2⎟
    ⎝ 2    2    ⎠
    >>> expand(_[0]*g + _[1])
       2
    5⋅x  + 10⋅x + 3

As you can see, ``q`` has a non-integer coefficient. If you want to do division
only in the ring of polynomials with integer coefficients, you can specify an
additional parameter

    >>> div(f, g, field=False)
    ⎛      2    ⎞
    ⎝5, 5⋅x  - 7⎠

But be warned, that this ring is no longer Euclidean and that the degree of the
remainder doesn't need to be smaller than that of ``f``. Since 2 doesn't divide 5,
`2 x` doesn't divide `5 x^2`, even if the degree is smaller. But

    >>> g = 5*x + 1

    >>> div(f, g, field=False)
    (x, 9⋅x + 3)
    >>> expand(_[0]*g + _[1])
       2
    5⋅x  + 10⋅x + 3

This also works for polynomials with multiple variables

    >>> div(x*y + y*z, 3*x + 3*z)
    ⎛y   ⎞
    ⎜─, 0⎟
    ⎝3   ⎠

GCD and LCM
===========

With division, there is also the computation of the greatest common divisor and
the least common multiple.

When the polynomials have integer coefficients, the contents' gcd is also
considered

    >>> gcd((12*x + 12)*x, 16*x**2)
    4⋅x

But if the polynomials have rational coefficients, then the returned polynomial is
monic

    >>> gcd(3*x**2/2, 9*x/4)
    x

Symbolic exponents are supported

    >>> gcd(2*x**(n + 4) - x**(n + 2), 4*x**(n + 1) + 3*x**n)
     n
    x

It also works with multiple variables. In this case, the variables are ordered
alphabetically, be default, which has influence on the leading coefficient

    >>> gcd(x*y/2 + y**2, 3*x + 6*y)
    x + 2⋅y

The lcm is connected with the gcd and one can be computed using the other

    >>> f, g = x*y**2 + x**2*y, x**2*y**2

    >>> gcd(f, g)
    x⋅y
    >>> lcm(f, g)
     3  2    2  3
    x ⋅y  + x ⋅y
    >>> expand(f*g)
     4  3    3  4
    x ⋅y  + x ⋅y
    >>> expand(gcd(f, g, x, y)*lcm(f, g, x, y))
     4  3    3  4
    x ⋅y  + x ⋅y

Square-free factorization
=========================

The square-free factorization of a univariate polynomial is the product of all
factors (not necessarily irreducible) of degree 1, 2 etc

    >>> sqf(2*x**2 + 5*x**3 + 4*x**4 + x**5)
                    2
            ⎛ 2    ⎞
    (x + 2)⋅⎝x  + x⎠

Factorization
=============

Factorization supported over different domains, lets compute one for the
rational field, its algebraic extension or the finite field of order 5

    >>> f = x**4 - 3*x**2 + 1

    >>> factor(f)
    ⎛ 2        ⎞ ⎛ 2        ⎞
    ⎝x  - x - 1⎠⋅⎝x  + x - 1⎠
    >>> factor(f, extension=GoldenRatio)
    (x - φ)⋅(x + φ)⋅(x - 1 + φ)⋅(x - φ + 1)
    >>> factor(f, modulus=5)
           2        2
    (x + 2) ⋅(x + 3)

The finite fields of prime power order are supported

    >>> factor(x**3 + 3*x + 2, modulus=4)
            ⎛ 2        ⎞
    (x + 1)⋅⎝x  + x + 2⎠

You also may use ``gaussian`` keyword to obtain a factorization over
Gaussian rationals

    >>> factor(4*x**4 + 8*x**3 + 77*x**2 + 18*x + 153, gaussian=True)
      ⎛    3⋅ⅈ⎞ ⎛    3⋅ⅈ⎞
    4⋅⎜x - ───⎟⋅⎜x + ───⎟⋅(x + 1 - 4⋅ⅈ)⋅(x + 1 + 4⋅ⅈ)
      ⎝     2 ⎠ ⎝     2 ⎠

Computing with multivariate polynomials over various domains is as simple as in
univariate case.

    >>> factor(x**2 + 4*x*y + 4*y**2)
             2
    (x + 2⋅y)
    >>> factor(x**3 + y**3, extension=sqrt(-3))
            ⎛      ⎛        ___  ⎞⎞ ⎛      ⎛        ___  ⎞⎞
            ⎜      ⎜  1   ╲╱ 3 ⋅ⅈ⎟⎟ ⎜      ⎜  1   ╲╱ 3 ⋅ⅈ⎟⎟
    (x + y)⋅⎜x + y⋅⎜- ─ - ───────⎟⎟⋅⎜x + y⋅⎜- ─ + ───────⎟⎟
            ⎝      ⎝  2      2   ⎠⎠ ⎝      ⎝  2      2   ⎠⎠

Gröbner bases
=============

Buchberger's algorithm is implemented, supporting various monomial orders

    >>> groebner([x**2 + 1, y**4*x + x**3])
                 ⎛⎡ 2       4    ⎤                           ⎞
    GroebnerBasis⎝⎣x  + 1, y  - 1⎦, x, y, domain=ℤ, order=lex⎠


    >>> groebner([x**2 + 1, y**4*x + x**3, x*y*z**3], order=grevlex)
                 ⎛⎡ 4       3   2    ⎤                                  ⎞
    GroebnerBasis⎝⎣y  - 1, z , x  + 1⎦, x, y, z, domain=ℤ, order=grevlex⎠
