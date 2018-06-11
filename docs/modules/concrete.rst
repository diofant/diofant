Concrete Mathematics
====================

Hypergeometric terms
--------------------

The center stage, in recurrence solving and summations, play hypergeometric
terms. Formally these are sequences annihilated by first order linear
recurrence operators. In simple words if we are given term `a(n)` then it is
hypergeometric if its consecutive term ratio is a rational function in `n`.

To check if a sequence is of this type you can use the ``is_hypergeometric``
method which is available in Basic class. Here is simple example involving a
polynomial:

    >>> (n**2 + 1).is_hypergeometric(n)
    True

Of course polynomials are hypergeometric but are there any more complicated
sequences of this type? Here are some trivial examples:

    >>> factorial(n).is_hypergeometric(n)
    True
    >>> binomial(n, k).is_hypergeometric(n)
    True
    >>> rf(n, k).is_hypergeometric(n)
    True
    >>> ff(n, k).is_hypergeometric(n)
    True
    >>> gamma(n).is_hypergeometric(n)
    True
    >>> (2**n).is_hypergeometric(n)
    True

We see that all species used in summations and other parts of concrete
mathematics are hypergeometric. Note also that binomial coefficients and both
rising and falling factorials are hypergeometric in both their arguments:

    >>> binomial(n, k).is_hypergeometric(k)
    True
    >>> rf(n, k).is_hypergeometric(k)
    True
    >>> ff(n, k).is_hypergeometric(k)
    True

To say more, all previously shown examples are valid for integer linear
arguments:

    >>> factorial(2*n).is_hypergeometric(n)
    True
    >>> binomial(3*n+1, k).is_hypergeometric(n)
    True
    >>> rf(n+1, k-1).is_hypergeometric(n)
    True
    >>> ff(n-1, k+1).is_hypergeometric(n)
    True
    >>> gamma(5*n).is_hypergeometric(n)
    True
    >>> (2**(n-7)).is_hypergeometric(n)
    True

However nonlinear arguments make those sequences fail to be hypergeometric:

    >>> factorial(n**2).is_hypergeometric(n)
    False
    >>> (2**(n**3 + 1)).is_hypergeometric(n)
    False

If not only the knowledge of being hypergeometric or not is needed, you can use
``hypersimp()`` function. It will try to simplify combinatorial expression and
if the term given is hypergeometric it will return a quotient of polynomials of
minimal degree. Otherwise is will return `None` to say that sequence is not
hypergeometric:

    >>> hypersimp(factorial(2*n), n)
    2*(n + 1)*(2*n + 1)
    >>> hypersimp(factorial(n**2), n)


Concrete Class Reference
------------------------
.. autoclass:: diofant.concrete.summations.Sum
   :members:

.. autoclass:: diofant.concrete.products.Product
   :members:

.. autoclass:: diofant.concrete.expr_with_limits.ExprWithLimits
   :members:

.. autoclass:: diofant.concrete.expr_with_intlimits.ExprWithIntLimits
   :members:

Concrete Functions Reference
----------------------------

.. autofunction:: diofant.concrete.summations.summation

.. autofunction:: diofant.concrete.products.product

.. autofunction:: diofant.concrete.gosper.gosper_normal

.. autofunction:: diofant.concrete.gosper.gosper_term

.. autofunction:: diofant.concrete.gosper.gosper_sum
