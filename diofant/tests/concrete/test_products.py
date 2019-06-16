import pytest

from diofant import Function
from diofant import KroneckerDelta as Kd
from diofant import (Product, Rational, Sum, Symbol, cos, exp, factorial, log,
                     oo, product, rf, simplify, sqrt, symbols)
from diofant.concrete.expr_with_intlimits import ReorderError


__all__ = ()

a, i, k, n, m, x = symbols('a,i,k,n,m,x', integer=True)
f = Function('f')


def test_karr_convention():
    # Test the Karr product convention that we want to hold.
    # See his paper "Summation in Finite Terms" for a detailed
    # reasoning why we really want exactly this definition.
    # The convention is described for sums on page 309 and
    # essentially in section 1.4, definition 3. For products
    # we can find in analogy:
    #
    # \prod_{m <= i < n} f(i) 'has the obvious meaning'      for m < n
    # \prod_{m <= i < n} f(i) = 0                            for m = n
    # \prod_{m <= i < n} f(i) = 1 / \prod_{n <= i < m} f(i)  for m > n
    #
    # It is important to note that he defines all products with
    # the upper limit being *exclusive*.
    # In contrast, diofant and the usual mathematical notation has:
    #
    # prod_{i = a}^b f(i) = f(a) * f(a+1) * ... * f(b-1) * f(b)
    #
    # with the upper limit *inclusive*. So translating between
    # the two we find that:
    #
    # \prod_{m <= i < n} f(i) = \prod_{i = m}^{n-1} f(i)
    #
    # where we intentionally used two different ways to typeset the
    # products and its limits.

    i = Symbol("i", integer=True)
    k = Symbol("k", integer=True)
    j = Symbol("j", integer=True)

    # A simple example with a concrete factors and symbolic limits.

    # The normal product: m = k and n = k + j and therefore m < n:
    m = k
    n = k + j

    a = m
    b = n - 1
    S1 = Product(i**2, (i, a, b)).doit()

    # The reversed product: m = k + j and n = k and therefore m > n:
    m = k + j
    n = k

    a = m
    b = n - 1
    S2 = Product(i**2, (i, a, b)).doit()

    assert simplify(S1 * S2) == 1

    # Test the empty product: m = k and n = k and therefore m = n:
    m = k
    n = k

    a = m
    b = n - 1
    Sz = Product(i**2, (i, a, b)).doit()

    assert Sz == 1

    # Another example this time with an unspecified factor and
    # numeric limits. (We can not do both tests in the same example.)
    f = Function("f")

    # The normal product with m < n:
    m = 2
    n = 11

    a = m
    b = n - 1
    S1 = Product(f(i), (i, a, b)).doit()

    # The reversed product with m > n:
    m = 11
    n = 2

    a = m
    b = n - 1
    S2 = Product(f(i), (i, a, b)).doit()

    assert simplify(S1 * S2) == 1

    # Test the empty product with m = n:
    m = 5
    n = 5

    a = m
    b = n - 1
    Sz = Product(f(i), (i, a, b)).doit()

    assert Sz == 1


def test_karr_proposition_2a():
    # Test Karr, page 309, proposition 2, part a
    i = Symbol("i", integer=True)
    u = Symbol("u", integer=True)
    v = Symbol("v", integer=True)

    def test_the_product(m, n):
        # g
        g = i**3 + 2*i**2 - 3*i
        # f = Delta g
        f = simplify(g.subs({i: i + 1})/g)
        # The product
        a = m
        b = n - 1
        P = Product(f, (i, a, b)).doit()
        # Test if Product_{m <= i < n} f(i) = g(n) / g(m)
        assert simplify(P/(g.subs({i: n})/g.subs({i: m}))) == 1

    # m < n
    test_the_product(u,   u+v)
    # m = n
    test_the_product(u,   u  )
    # m > n
    test_the_product(u+v, u  )


def test_karr_proposition_2b():
    # Test Karr, page 309, proposition 2, part b
    i = Symbol("i", integer=True)
    u = Symbol("u", integer=True)
    v = Symbol("v", integer=True)
    w = Symbol("w", integer=True)

    def test_the_product(l, n, m):
        # Productmand
        s = i**3
        # First product
        a = l
        b = n - 1
        S1 = Product(s, (i, a, b)).doit()
        # Second product
        a = l
        b = m - 1
        S2 = Product(s, (i, a, b)).doit()
        # Third product
        a = m
        b = n - 1
        S3 = Product(s, (i, a, b)).doit()
        # Test if S1 = S2 * S3 as required
        assert simplify(S1 / (S2 * S3)) == 1

    # l < m < n
    test_the_product(u,     u+v,   u+v+w)
    # l < m = n
    test_the_product(u,     u+v,   u+v  )
    # l < m > n
    test_the_product(u,     u+v+w, v    )
    # l = m < n
    test_the_product(u,     u,     u+v  )
    # l = m = n
    test_the_product(u,     u,     u    )
    # l = m > n
    test_the_product(u+v,   u+v,   u    )
    # l > m < n
    test_the_product(u+v,   u,     u+w  )
    # l > m = n
    test_the_product(u+v,   u,     u    )
    # l > m > n
    test_the_product(u+v+w, u+v,   u    )


def test_simple_products():
    assert product(2, (k, a, n)) == 2**(n - a + 1)
    assert product(k, (k, 1, n)) == factorial(n)
    assert product(k**3, (k, 1, n)) == factorial(n)**3

    assert product(k + 1, (k, 0, n - 1)) == factorial(n)
    assert product(k + 1, (k, a, n - 1)) == rf(1 + a, n - a)

    assert product(cos(k), (k, 0, 5)) == cos(1)*cos(2)*cos(3)*cos(4)*cos(5)
    assert product(cos(k), (k, 3, 5)) == cos(3)*cos(4)*cos(5)
    assert product(cos(k), (k, 1, Rational(5, 2))) != cos(1)*cos(2)

    assert isinstance(product(k**k, (k, 1, n)), Product)

    assert Product(x**k, (k, 1, n)).variables == [k]

    pytest.raises(ValueError, lambda: Product(n))
    pytest.raises(ValueError, lambda: Product(n*k))
    pytest.raises(ValueError, lambda: Product(n, k))
    pytest.raises(ValueError, lambda: Product(n, k, 1))
    pytest.raises(ValueError, lambda: Product(n, k, 1, 10))
    pytest.raises(ValueError, lambda: Product(n, (k, 1)))

    assert product(1, (n, 1, oo)) == 1  # issue sympy/sympy#8301
    assert product(2, (n, 1, oo)) == oo
    assert isinstance(product(-1, (n, 1, oo)), Product)

    assert product(Kd(n, m), (m, 1, 3)) == 0
    assert product(Kd(n, m), (m, 1, 1)) == Kd(n, 1)


def test_multiple_products():
    assert product(x, (n, 1, k), (k, 1, m)) == x**(m**2/2 + m/2)
    assert product(f(n), (
        n, 1, m), (m, 1, k)) == Product(f(n), (n, 1, m), (m, 1, k)).doit()
    assert Product(f(n), (m, 1, k), (n, 1, k)).doit() == \
        Product(Product(f(n), (m, 1, k)), (n, 1, k)).doit() == \
        product(f(n), (m, 1, k), (n, 1, k)) == \
        product(product(f(n), (m, 1, k)), (n, 1, k)) == \
        Product(f(n)**k, (n, 1, k))
    assert Product(
        x, (x, 1, k), (k, 1, n)).doit() == Product(factorial(k), (k, 1, n))

    assert Product(x**k, (n, 1, k), (k, 1, m)).variables == [n, k]


def test_rational_products():
    assert product(1 + 1/k, (k, 1, n)) == rf(2, n)/factorial(n)


def test_special_products():
    # Wallis product
    assert product((4*k)**2 / (4*k**2 - 1), (k, 1, n)) == \
        4**n*factorial(n)**2/rf(Rational(1, 2), n)/rf(Rational(3, 2), n)

    # Euler's product formula for sin
    assert product(1 + a/k**2, (k, 1, n)) == \
        rf(1 - sqrt(-a), n)*rf(1 + sqrt(-a), n)/factorial(n)**2


def test__eval_product():
    # issue sympy/sympy#4809
    a = Function('a')
    assert product(2*a(i), (i, 1, n)) == 2**n * Product(a(i), (i, 1, n))
    # issue sympy/sympy#4810
    assert product(2**i, (i, 1, n)) == 2**(n/2 + n**2/2)

    assert (product((i - 1)*(i**6 + i - 1), (i, n, m)) ==
            rf(n - 1, m - n + 1)*product(i**6 + i - 1, (i, n, m)))

    assert (product(log(i)**2*cos(i)**3, (i, n, m)) ==
            Product(log(i)**2*cos(i)**3, (i, n, m)))


def test_product_pow():
    # issue sympy/sympy#4817
    assert product(2**f(k), (k, 1, n)) == 2**Sum(f(k), (k, 1, n))
    assert product(2**(2*f(k)), (k, 1, n)) == 2**Sum(2*f(k), (k, 1, n))


def test_infinite_product():
    # issue sympy/sympy#5737
    assert isinstance(Product(2**(1/factorial(n)), (n, 0, oo)), Product)


def test_conjugate_transpose():
    p = Product(x**k, (k, 1, 3))
    assert p.adjoint().doit() == p.doit().adjoint()
    assert p.conjugate().doit() == p.doit().conjugate()
    assert p.transpose().doit() == p.doit().transpose()

    A, B = symbols("A B", commutative=False)
    p = Product(A*B**k, (k, 1, 3))
    assert p.adjoint().doit() == p.doit().adjoint()
    assert p.conjugate().doit() == p.doit().conjugate()
    assert p.transpose().doit() == p.doit().transpose()


def test_simplify():
    y, t, b, c = symbols('y, t, b, c', integer=True)

    assert simplify(Product(x*y, (x, n, m), (y, a, k)) *
                    Product(y, (x, n, m), (y, a, k))) == \
        Product(x*y**2, (x, n, m), (y, a, k))
    assert simplify(3 * y * Product(x, (x, n, m)) * Product(x, (x, m + 1, a))) \
        == 3 * y * Product(x, (x, n, a))
    assert simplify(Product(x, (x, k + 1, a)) * Product(x, (x, n, k))) == \
        Product(x, (x, n, a))
    assert simplify(Product(x, (x, k + 1, a)) * Product(x + 1, (x, n, k))) == \
        Product(x, (x, k + 1, a)) * Product(x + 1, (x, n, k))
    assert simplify(Product(x, (t, a, b)) * Product(y, (t, a, b)) *
                    Product(x, (t, b+1, c))) == Product(x*y, (t, a, b)) * \
        Product(x, (t, b+1, c))
    assert simplify(Product(x, (t, a, b)) * Product(x, (t, b+1, c)) *
                    Product(y, (t, a, b))) == Product(x*y, (t, a, b)) * \
        Product(x, (t, b+1, c))


def test_change_index():
    b, y, c, d, z = symbols('b, y, c, d, z', integer=True)

    assert Product(x, (x, a, b)).change_index(x, x + 1, y) == \
        Product(y - 1, (y, a + 1, b + 1))
    assert Product(x**2, (x, a, b)).change_index(x, x - 1) == \
        Product((x + 1)**2, (x, a - 1, b - 1))
    assert Product(x**2, (x, a, b)).change_index(x, -x, y) == \
        Product((-y)**2, (y, -b, -a))
    assert Product(x, (x, a, b)).change_index(x, -x - 1) == \
        Product(-x - 1, (x, - b - 1, -a - 1))
    assert Product(x*y, (x, a, b), (y, c, d)).change_index(x, x - 1, z) == \
        Product((z + 1)*y, (z, a - 1, b - 1), (y, c, d))

    p = Product(x, (x, 1, 5))
    pytest.raises(ValueError, lambda: p.change_index(x, x**2, y))
    pytest.raises(ValueError, lambda: p.change_index(x, 2*x, y))


def test_index():
    p = Product(x, (x, 1, 3), (x, 1, 4))
    pytest.raises(ValueError, lambda: p.index(x))


def test_reorder():
    b, y, c, d, z = symbols('b, y, c, d, z', integer=True)

    assert Product(x*y, (x, a, b), (y, c, d)).reorder((0, 1)) == \
        Product(x*y, (y, c, d), (x, a, b))
    assert Product(x, (x, a, b), (x, c, d)).reorder((0, 1)) == \
        Product(x, (x, c, d), (x, a, b))
    assert Product(x*y + z, (x, a, b), (z, m, n), (y, c, d)).reorder(
        (2, 0), (0, 1)) == Product(x*y + z, (z, m, n), (y, c, d), (x, a, b))
    assert Product(x*y*z, (x, a, b), (y, c, d), (z, m, n)).reorder(
        (0, 1), (1, 2), (0, 2)) == \
        Product(x*y*z, (x, a, b), (z, m, n), (y, c, d))
    assert Product(x*y*z, (x, a, b), (y, c, d), (z, m, n)).reorder(
        (x, y), (y, z), (x, z)) == \
        Product(x*y*z, (x, a, b), (z, m, n), (y, c, d))
    assert Product(x*y, (x, a, b), (y, c, d)).reorder((x, 1)) == \
        Product(x*y, (y, c, d), (x, a, b))
    assert Product(x*y, (x, a, b), (y, c, d)).reorder((y, x)) == \
        Product(x*y, (y, c, d), (x, a, b))

    pytest.raises(ValueError,
                  lambda: Product(x*y, (x, a, b), (y, c, d)).reorder((x,)))


def test_reverse_order():
    x, y, a, b, c, d = symbols('x, y, a, b, c, d', integer=True)

    assert Product(x, (x, 0, 3)).reverse_order(0) == Product(1/x, (x, 4, -1))
    assert Product(x*y, (x, 1, 5), (y, 0, 6)).reverse_order(0, 1) == \
        Product(x*y, (x, 6, 0), (y, 7, -1))
    assert Product(x, (x, 1, 2)).reverse_order(0) == Product(1/x, (x, 3, 0))
    assert Product(x, (x, 1, 3)).reverse_order(0) == Product(1/x, (x, 4, 0))
    assert Product(x, (x, 1, a)).reverse_order(0) == Product(1/x, (x, a + 1, 0))
    assert Product(x, (x, a, 5)).reverse_order(0) == Product(1/x, (x, 6, a - 1))
    assert Product(x, (x, a + 1, a + 5)).reverse_order(0) == \
        Product(1/x, (x, a + 6, a))
    assert Product(x, (x, a + 1, a + 2)).reverse_order(0) == \
        Product(1/x, (x, a + 3, a))
    assert Product(x, (x, a + 1, a + 1)).reverse_order(0) == \
        Product(1/x, (x, a + 2, a))
    assert Product(x, (x, a, b)).reverse_order(0) == Product(1/x, (x, b + 1, a - 1))
    assert Product(x, (x, a, b)).reverse_order(x) == Product(1/x, (x, b + 1, a - 1))
    assert Product(x*y, (x, a, b), (y, 2, 5)).reverse_order(x, 1) == \
        Product(x*y, (x, b + 1, a - 1), (y, 6, 1))
    assert Product(x*y, (x, a, b), (y, 2, 5)).reverse_order(y, x) == \
        Product(x*y, (x, b + 1, a - 1), (y, 6, 1))
    p = Product(x, (x, 0, 3))
    assert p.reverse_order(10) == p


def test_reorder_limit():
    x, y, a, b, c, d = symbols('x, y, a, b, c, d', integer=True)
    pytest.raises(ReorderError,
                  lambda: Product(x**2, (x, a, b), (y, x, d)).reorder_limit(1, 0))


def test_sympyissue_9983():
    p = Product(1 + 1/n**Rational(2, 3), (n, 1, oo))
    assert p == p.doit()


def test_rewrite_Sum():
    assert Product(1 - 1/(4*k**2), (k, 1, oo)).rewrite(Sum) == \
        exp(Sum(log(1 - 1/(4*k**2)), (k, 1, oo)))
