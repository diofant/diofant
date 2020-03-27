from diofant import (And, Eq, Integer, Piecewise, Sum, exp, oo, piecewise_fold,
                     symbols)
from diofant.concrete.delta import deltaproduct as dp
from diofant.concrete.delta import deltasummation as ds
from diofant.functions import KroneckerDelta as Kd


__all__ = ()

i, j, k, l, m = symbols('i j k l m', integer=True, finite=True)
x, y = symbols('x y', commutative=False)


def test_deltaproduct_trivial():
    assert dp(x, (j, 1, 0)) == 1
    assert dp(x, (j, 1, 3)) == x**3
    assert dp(x + y, (j, 1, 3)) == (x + y)**3
    assert dp(x*y, (j, 1, 3)) == (x*y)**3
    assert dp(Kd(i, j), (k, 1, 3)) == Kd(i, j)
    assert dp(x*Kd(i, j), (k, 1, 3)) == x**3*Kd(i, j)
    assert dp(x*y*Kd(i, j), (k, 1, 3)) == (x*y)**3*Kd(i, j)


def test_deltaproduct_basic():
    assert dp(Kd(i, j), (j, 1, 3)) == 0
    assert dp(Kd(i, j), (j, 1, 1)) == Kd(i, 1)
    assert dp(Kd(i, j), (j, 2, 2)) == Kd(i, 2)
    assert dp(Kd(i, j), (j, 3, 3)) == Kd(i, 3)
    assert dp(Kd(i, j), (j, 1, k)) == Kd(i, 1)*Kd(k, 1) + Kd(k, 0)
    assert dp(Kd(i, j), (j, k, 3)) == Kd(i, 3)*Kd(k, 3) + Kd(k, 4)
    assert dp(Kd(i, j), (j, k, l)) == Kd(i, l)*Kd(k, l) + Kd(k, l + 1)
    assert dp(Kd(i, 1), (i, j**2, k**2)) == (Kd(1, j**2)*Kd(j**2, k**2) +
                                             Kd(k**2, j**2 - 1))


def test_deltaproduct_mul_x_kd():
    assert dp(x*Kd(i, j), (j, 1, 3)) == 0
    assert dp(x*Kd(i, j), (j, 1, 1)) == x*Kd(i, 1)
    assert dp(x*Kd(i, j), (j, 2, 2)) == x*Kd(i, 2)
    assert dp(x*Kd(i, j), (j, 3, 3)) == x*Kd(i, 3)
    assert dp(x*Kd(i, j), (j, 1, k)) == x*Kd(i, 1)*Kd(k, 1) + Kd(k, 0)
    assert dp(x*Kd(i, j), (j, k, 3)) == x*Kd(i, 3)*Kd(k, 3) + Kd(k, 4)
    assert dp(x*Kd(i, j), (j, k, l)) == x*Kd(i, l)*Kd(k, l) + Kd(k, l + 1)


def test_deltaproduct_mul_add_x_y_kd():
    assert dp((x + y)*Kd(i, j), (j, 1, 3)) == 0
    assert dp((x + y)*Kd(i, j), (j, 1, 1)) == (x + y)*Kd(i, 1)
    assert dp((x + y)*Kd(i, j), (j, 2, 2)) == (x + y)*Kd(i, 2)
    assert dp((x + y)*Kd(i, j), (j, 3, 3)) == (x + y)*Kd(i, 3)
    assert dp((x + y)*Kd(i, j), (j, 1, k)) == \
        (x + y)*Kd(i, 1)*Kd(k, 1) + Kd(k, 0)
    assert dp((x + y)*Kd(i, j), (j, k, 3)) == \
        (x + y)*Kd(i, 3)*Kd(k, 3) + Kd(k, 4)
    assert dp((x + y)*Kd(i, j), (j, k, l)) == \
        (x + y)*Kd(i, l)*Kd(k, l) + Kd(k, l + 1)


def test_deltaproduct_add_kd_kd():
    assert dp(Kd(i, k) + Kd(j, k), (k, 1, 3)) == 0
    assert dp(Kd(i, k) + Kd(j, k), (k, 1, 1)) == Kd(i, 1) + Kd(j, 1)
    assert dp(Kd(i, k) + Kd(j, k), (k, 2, 2)) == Kd(i, 2) + Kd(j, 2)
    assert dp(Kd(i, k) + Kd(j, k), (k, 3, 3)) == Kd(i, 3) + Kd(j, 3)
    assert dp(Kd(i, k) + Kd(j, k), (k, 1, l)) == Kd(l, 0) + \
        Kd(i, 1)*Kd(l, 1) + Kd(j, 1)*Kd(l, 1) + \
        Kd(i, 1)*Kd(j, 2)*Kd(l, 2) + Kd(j, 1)*Kd(i, 2)*Kd(l, 2)
    assert dp(Kd(i, k) + Kd(j, k), (k, l, 3)) == Kd(l, 4) + \
        Kd(i, 3)*Kd(l, 3) + Kd(j, 3)*Kd(l, 3) + \
        Kd(i, 2)*Kd(j, 3)*Kd(l, 2) + Kd(i, 3)*Kd(j, 2)*Kd(l, 2)
    assert dp(Kd(i, k) + Kd(j, k), (k, l, m)) == Kd(l, m + 1) + \
        Kd(i, m)*Kd(l, m) + Kd(j, m)*Kd(l, m) + \
        Kd(i, m)*Kd(j, m - 1)*Kd(l, m - 1) + Kd(i, m - 1)*Kd(j, m)*Kd(l, m - 1)


def test_deltaproduct_mul_x_add_kd_kd():
    assert dp(x*(Kd(i, k) + Kd(j, k)), (k, 1, 3)) == 0
    assert dp(x*(Kd(i, k) + Kd(j, k)), (k, 1, 1)) == x*(Kd(i, 1) + Kd(j, 1))
    assert dp(x*(Kd(i, k) + Kd(j, k)), (k, 2, 2)) == x*(Kd(i, 2) + Kd(j, 2))
    assert dp(x*(Kd(i, k) + Kd(j, k)), (k, 3, 3)) == x*(Kd(i, 3) + Kd(j, 3))
    assert dp(x*(Kd(i, k) + Kd(j, k)), (k, 1, l)) == Kd(l, 0) + \
        x*Kd(i, 1)*Kd(l, 1) + x*Kd(j, 1)*Kd(l, 1) + \
        x**2*Kd(i, 1)*Kd(j, 2)*Kd(l, 2) + x**2*Kd(j, 1)*Kd(i, 2)*Kd(l, 2)
    assert dp(x*(Kd(i, k) + Kd(j, k)), (k, l, 3)) == Kd(l, 4) + \
        x*Kd(i, 3)*Kd(l, 3) + x*Kd(j, 3)*Kd(l, 3) + \
        x**2*Kd(i, 2)*Kd(j, 3)*Kd(l, 2) + x**2*Kd(i, 3)*Kd(j, 2)*Kd(l, 2)
    assert dp(x*(Kd(i, k) + Kd(j, k)), (k, l, m)) == Kd(l, m + 1) + \
        x*Kd(i, m)*Kd(l, m) + x*Kd(j, m)*Kd(l, m) + \
        x**2*Kd(i, m - 1)*Kd(j, m)*Kd(l, m - 1) + \
        x**2*Kd(i, m)*Kd(j, m - 1)*Kd(l, m - 1)


def test_deltaproduct_mul_add_x_y_add_kd_kd():
    assert dp((x + y)*(Kd(i, k) + Kd(j, k)), (k, 1, 3)) == 0
    assert dp((x + y)*(Kd(i, k) + Kd(j, k)), (k, 1, 1)) == \
        (x + y)*(Kd(i, 1) + Kd(j, 1))
    assert dp((x + y)*(Kd(i, k) + Kd(j, k)), (k, 2, 2)) == \
        (x + y)*(Kd(i, 2) + Kd(j, 2))
    assert dp((x + y)*(Kd(i, k) + Kd(j, k)), (k, 3, 3)) == \
        (x + y)*(Kd(i, 3) + Kd(j, 3))
    assert dp((x + y)*(Kd(i, k) + Kd(j, k)), (k, 1, l)) == Kd(l, 0) + \
        (x + y)*Kd(i, 1)*Kd(l, 1) + (x + y)*Kd(j, 1)*Kd(l, 1) + \
        (x + y)**2*Kd(i, 1)*Kd(j, 2)*Kd(l, 2) + \
        (x + y)**2*Kd(j, 1)*Kd(i, 2)*Kd(l, 2)
    assert dp((x + y)*(Kd(i, k) + Kd(j, k)), (k, l, 3)) == Kd(l, 4) + \
        (x + y)*Kd(i, 3)*Kd(l, 3) + (x + y)*Kd(j, 3)*Kd(l, 3) + \
        (x + y)**2*Kd(i, 2)*Kd(j, 3)*Kd(l, 2) + \
        (x + y)**2*Kd(i, 3)*Kd(j, 2)*Kd(l, 2)
    assert dp((x + y)*(Kd(i, k) + Kd(j, k)), (k, l, m)) == Kd(l, m + 1) + \
        (x + y)*Kd(i, m)*Kd(l, m) + (x + y)*Kd(j, m)*Kd(l, m) + \
        (x + y)**2*Kd(i, m - 1)*Kd(j, m)*Kd(l, m - 1) + \
        (x + y)**2*Kd(i, m)*Kd(j, m - 1)*Kd(l, m - 1)


def test_deltaproduct_add_mul_x_y_mul_x_kd():
    assert dp(x*y + x*Kd(i, j), (j, 1, 3)) == (x*y)**3 + \
        x*(x*y)**2*Kd(i, 1) + (x*y)*x*(x*y)*Kd(i, 2) + (x*y)**2*x*Kd(i, 3)
    assert dp(x*y + x*Kd(i, j), (j, 1, 1)) == x*y + x*Kd(i, 1)
    assert dp(x*y + x*Kd(i, j), (j, 2, 2)) == x*y + x*Kd(i, 2)
    assert dp(x*y + x*Kd(i, j), (j, 3, 3)) == x*y + x*Kd(i, 3)
    assert dp(x*y + x*Kd(i, j), (j, 1, k)) == \
        (x*y)**k + Piecewise(
            ((x*y)**(i - 1)*x*(x*y)**(k - i), And(Integer(1) <= i, i <= k)),
            (0, True)
    )
    assert dp(x*y + x*Kd(i, j), (j, k, 3)) == \
        (x*y)**(-k + 4) + Piecewise(
            ((x*y)**(i - k)*x*(x*y)**(3 - i), And(k <= i, i <= 3)),
            (0, True)
    )
    assert dp(x*y + x*Kd(i, j), (j, k, l)) == \
        (x*y)**(-k + l + 1) + Piecewise(
            ((x*y)**(i - k)*x*(x*y)**(l - i), And(k <= i, i <= l)),
            (0, True)
    )


def test_deltaproduct_mul_x_add_y_kd():
    assert dp(x*(y + Kd(i, j)), (j, 1, 3)) == (x*y)**3 + \
        x*(x*y)**2*Kd(i, 1) + (x*y)*x*(x*y)*Kd(i, 2) + (x*y)**2*x*Kd(i, 3)
    assert dp(x*(y + Kd(i, j)), (j, 1, 1)) == x*(y + Kd(i, 1))
    assert dp(x*(y + Kd(i, j)), (j, 2, 2)) == x*(y + Kd(i, 2))
    assert dp(x*(y + Kd(i, j)), (j, 3, 3)) == x*(y + Kd(i, 3))
    assert dp(x*(y + Kd(i, j)), (j, 1, k)) == \
        (x*y)**k + Piecewise(
            ((x*y)**(i - 1)*x*(x*y)**(k - i), And(Integer(1) <= i, i <= k)),
            (0, True)
    )
    assert dp(x*(y + Kd(i, j)), (j, k, 3)) == \
        (x*y)**(-k + 4) + Piecewise(
            ((x*y)**(i - k)*x*(x*y)**(3 - i), And(k <= i, i <= 3)),
            (0, True)
    )
    assert dp(x*(y + Kd(i, j)), (j, k, l)) == \
        (x*y)**(-k + l + 1) + Piecewise(
            ((x*y)**(i - k)*x*(x*y)**(l - i), And(k <= i, i <= l)),
            (0, True)
    )


def test_deltaproduct_mul_x_add_y_twokd():
    assert dp(x*(y + 2*Kd(i, j)), (j, 1, 3)) == (x*y)**3 + \
        2*x*(x*y)**2*Kd(i, 1) + 2*x*y*x*x*y*Kd(i, 2) + 2*(x*y)**2*x*Kd(i, 3)
    assert dp(x*(y + 2*Kd(i, j)), (j, 1, 1)) == x*(y + 2*Kd(i, 1))
    assert dp(x*(y + 2*Kd(i, j)), (j, 2, 2)) == x*(y + 2*Kd(i, 2))
    assert dp(x*(y + 2*Kd(i, j)), (j, 3, 3)) == x*(y + 2*Kd(i, 3))
    assert dp(x*(y + 2*Kd(i, j)), (j, 1, k)) == \
        (x*y)**k + Piecewise(
            (2*(x*y)**(i - 1)*x*(x*y)**(k - i), And(Integer(1) <= i, i <= k)),
            (0, True)
    )
    assert dp(x*(y + 2*Kd(i, j)), (j, k, 3)) == \
        (x*y)**(-k + 4) + Piecewise(
            (2*(x*y)**(i - k)*x*(x*y)**(3 - i), And(k <= i, i <= 3)),
            (0, True)
    )
    assert dp(x*(y + 2*Kd(i, j)), (j, k, l)) == \
        (x*y)**(-k + l + 1) + Piecewise(
            (2*(x*y)**(i - k)*x*(x*y)**(l - i), And(k <= i, i <= l)),
            (0, True)
    )


def test_deltaproduct_mul_add_x_y_add_y_kd():
    assert dp((x + y)*(y + Kd(i, j)), (j, 1, 3)) == ((x + y)*y)**3 + \
        (x + y)*((x + y)*y)**2*Kd(i, 1) + \
        (x + y)*y*(x + y)**2*y*Kd(i, 2) + \
        ((x + y)*y)**2*(x + y)*Kd(i, 3)
    assert dp((x + y)*(y + Kd(i, j)), (j, 1, 1)) == (x + y)*(y + Kd(i, 1))
    assert dp((x + y)*(y + Kd(i, j)), (j, 2, 2)) == (x + y)*(y + Kd(i, 2))
    assert dp((x + y)*(y + Kd(i, j)), (j, 3, 3)) == (x + y)*(y + Kd(i, 3))
    assert dp((x + y)*(y + Kd(i, j)), (j, 1, k)) == \
        ((x + y)*y)**k + Piecewise(
            (((x + y)*y)**(i - 1)*(x + y)*((x + y)*y)**(k - i),
             And(Integer(1) <= i, i <= k)),
            (0, True)
    )
    assert dp((x + y)*(y + Kd(i, j)), (j, k, 3)) == \
        ((x + y)*y)**(-k + 4) + Piecewise(
            (((x + y)*y)**(i - k)*(x + y)*((x + y)*y)**(3 - i),
             And(k <= i, i <= 3)),
            (0, True)
    )
    assert dp((x + y)*(y + Kd(i, j)), (j, k, l)) == \
        ((x + y)*y)**(-k + l + 1) + Piecewise(
            (((x + y)*y)**(i - k)*(x + y)*((x + y)*y)**(l - i),
             And(k <= i, i <= l)),
            (0, True)
    )


def test_deltaproduct_mul_add_x_kd_add_y_kd():
    assert dp((x + Kd(i, k))*(y + Kd(i, j)), (j, 1, 3)) == \
        Kd(i, 1)*(Kd(i, k) + x)*((Kd(i, k) + x)*y)**2 + \
        Kd(i, 2)*(Kd(i, k) + x)*y*(Kd(i, k) + x)**2*y + \
        Kd(i, 3)*((Kd(i, k) + x)*y)**2*(Kd(i, k) + x) + \
        ((Kd(i, k) + x)*y)**3
    assert dp((x + Kd(i, k))*(y + Kd(i, j)), (j, 1, 1)) == \
        (x + Kd(i, k))*(y + Kd(i, 1))
    assert dp((x + Kd(i, k))*(y + Kd(i, j)), (j, 2, 2)) == \
        (x + Kd(i, k))*(y + Kd(i, 2))
    assert dp((x + Kd(i, k))*(y + Kd(i, j)), (j, 3, 3)) == \
        (x + Kd(i, k))*(y + Kd(i, 3))
    assert dp((x + Kd(i, k))*(y + Kd(i, j)), (j, 1, k)) == \
        ((x + Kd(i, k))*y)**k + Piecewise(
            (((x + Kd(i, k))*y)**(i - 1)*(x + Kd(i, k)) *
             ((x + Kd(i, k))*y)**(-i + k), And(Integer(1) <= i, i <= k)),
            (0, True)
    )
    assert dp((x + Kd(i, k))*(y + Kd(i, j)), (j, k, 3)) == \
        ((x + Kd(i, k))*y)**(4 - k) + Piecewise(
            (((x + Kd(i, k))*y)**(i - k)*(x + Kd(i, k)) *
             ((x + Kd(i, k))*y)**(-i + 3), And(k <= i, i <= 3)),
            (0, True)
    )
    assert dp((x + Kd(i, k))*(y + Kd(i, j)), (j, k, l)) == \
        ((x + Kd(i, k))*y)**(-k + l + 1) + Piecewise(
            (((x + Kd(i, k))*y)**(i - k)*(x + Kd(i, k)) *
             ((x + Kd(i, k))*y)**(-i + l), And(k <= i, i <= l)),
            (0, True)
    )


def test_deltasummation_trivial():
    assert ds(x, (j, 1, 0)) == 0
    assert ds(x, (j, 1, 3)) == 3*x
    assert ds(x + y, (j, 1, 3)) == 3*(x + y)
    assert ds(x*y, (j, 1, 3)) == 3*x*y
    assert ds(Kd(i, j), (k, 1, 3)) == 3*Kd(i, j)
    assert ds(x*Kd(i, j), (k, 1, 3)) == 3*x*Kd(i, j)
    assert ds(x*y*Kd(i, j), (k, 1, 3)) == 3*x*y*Kd(i, j)


def test_deltasummation_basic_numerical():
    n = symbols('n', integer=True, nonzero=True)
    assert ds(Kd(n, 0), (n, 1, 3)) == 0

    # return unevaluated, until it gets implemented
    assert ds(Kd(i**2, j**2), (j, -oo, oo)) == \
        Sum(Kd(i**2, j**2), (j, -oo, oo))

    assert Piecewise((Kd(i, k), And(Integer(1) <= i, i <= 3)), (0, True)) == \
        ds(Kd(i, j)*Kd(j, k), (j, 1, 3)) == \
        ds(Kd(j, k)*Kd(i, j), (j, 1, 3))

    assert ds(Kd(i, k), (k, -oo, oo)) == 1
    assert ds(Kd(i, k), (k, 0, oo)) == Piecewise((1, Integer(0) <= i), (0, True))
    assert ds(Kd(i, k), (k, 1, 3)) == \
        Piecewise((1, And(Integer(1) <= i, i <= 3)), (0, True))
    assert ds(k*Kd(i, j)*Kd(j, k), (k, -oo, oo)) == j*Kd(i, j)
    assert ds(j*Kd(i, j), (j, -oo, oo)) == i
    assert ds(i*Kd(i, j), (i, -oo, oo)) == j
    assert ds(x, (i, 1, 3)) == 3*x
    assert ds((i + j)*Kd(i, j), (j, -oo, oo)) == 2*i


def test_deltasummation_basic_symbolic():
    assert ds(Kd(exp(i), 0), (i, 1, 3)) == 0
    assert ds(Kd(exp(i), 0), (i, -1, 3)) == 0
    assert ds(Kd(exp(i), 1), (i, 0, 3)) == 1
    assert ds(Kd(exp(i), 1), (i, 1, 3)) == 0
    assert ds(Kd(exp(i), 1), (i, -10, 3)) == 1
    assert ds(Kd(i, j), (j, 1, 3)) == \
        Piecewise((1, And(Integer(1) <= i, i <= 3)), (0, True))
    assert ds(Kd(i, j), (j, 1, 1)) == Piecewise((1, Eq(i, 1)), (0, True))
    assert ds(Kd(i, j), (j, 2, 2)) == Piecewise((1, Eq(i, 2)), (0, True))
    assert ds(Kd(i, j), (j, 3, 3)) == Piecewise((1, Eq(i, 3)), (0, True))
    assert ds(Kd(i, j), (j, 1, k)) == \
        Piecewise((1, And(Integer(1) <= i, i <= k)), (0, True))
    assert ds(Kd(i, j), (j, k, 3)) == \
        Piecewise((1, And(k <= i, i <= 3)), (0, True))
    assert ds(Kd(i, j), (j, k, l)) == \
        Piecewise((1, And(k <= i, i <= l)), (0, True))


def test_deltasummation_mul_x_kd():
    assert ds(x*Kd(i, j), (j, 1, 3)) == \
        Piecewise((x, And(Integer(1) <= i, i <= 3)), (0, True))
    assert ds(x*Kd(i, j), (j, 1, 1)) == Piecewise((x, Eq(i, 1)), (0, True))
    assert ds(x*Kd(i, j), (j, 2, 2)) == Piecewise((x, Eq(i, 2)), (0, True))
    assert ds(x*Kd(i, j), (j, 3, 3)) == Piecewise((x, Eq(i, 3)), (0, True))
    assert ds(x*Kd(i, j), (j, 1, k)) == \
        Piecewise((x, And(Integer(1) <= i, i <= k)), (0, True))
    assert ds(x*Kd(i, j), (j, k, 3)) == \
        Piecewise((x, And(k <= i, i <= 3)), (0, True))
    assert ds(x*Kd(i, j), (j, k, l)) == \
        Piecewise((x, And(k <= i, i <= l)), (0, True))


def test_deltasummation_mul_add_x_y_kd():
    assert ds((x + y)*Kd(i, j), (j, 1, 3)) == \
        Piecewise((x + y, And(Integer(1) <= i, i <= 3)), (0, True))
    assert ds((x + y)*Kd(i, j), (j, 1, 1)) == \
        Piecewise((x + y, Eq(i, 1)), (0, True))
    assert ds((x + y)*Kd(i, j), (j, 2, 2)) == \
        Piecewise((x + y, Eq(i, 2)), (0, True))
    assert ds((x + y)*Kd(i, j), (j, 3, 3)) == \
        Piecewise((x + y, Eq(i, 3)), (0, True))
    assert ds((x + y)*Kd(i, j), (j, 1, k)) == \
        Piecewise((x + y, And(Integer(1) <= i, i <= k)), (0, True))
    assert ds((x + y)*Kd(i, j), (j, k, 3)) == \
        Piecewise((x + y, And(k <= i, i <= 3)), (0, True))
    assert ds((x + y)*Kd(i, j), (j, k, l)) == \
        Piecewise((x + y, And(k <= i, i <= l)), (0, True))


def test_deltasummation_add_kd_kd():
    assert ds(Kd(i, k) + Kd(j, k), (k, 1, 3)) == piecewise_fold(
        Piecewise((1, And(Integer(1) <= i, i <= 3)), (0, True)) +
        Piecewise((1, And(Integer(1) <= j, j <= 3)), (0, True)))
    assert ds(Kd(i, k) + Kd(j, k), (k, 1, 1)) == piecewise_fold(
        Piecewise((1, Eq(i, 1)), (0, True)) +
        Piecewise((1, Eq(j, 1)), (0, True)))
    assert ds(Kd(i, k) + Kd(j, k), (k, 2, 2)) == piecewise_fold(
        Piecewise((1, Eq(i, 2)), (0, True)) +
        Piecewise((1, Eq(j, 2)), (0, True)))
    assert ds(Kd(i, k) + Kd(j, k), (k, 3, 3)) == piecewise_fold(
        Piecewise((1, Eq(i, 3)), (0, True)) +
        Piecewise((1, Eq(j, 3)), (0, True)))
    assert ds(Kd(i, k) + Kd(j, k), (k, 1, l)) == piecewise_fold(
        Piecewise((1, And(Integer(1) <= i, i <= l)), (0, True)) +
        Piecewise((1, And(Integer(1) <= j, j <= l)), (0, True)))
    assert ds(Kd(i, k) + Kd(j, k), (k, l, 3)) == piecewise_fold(
        Piecewise((1, And(l <= i, i <= 3)), (0, True)) +
        Piecewise((1, And(l <= j, j <= 3)), (0, True)))
    assert ds(Kd(i, k) + Kd(j, k), (k, l, m)) == piecewise_fold(
        Piecewise((1, And(l <= i, i <= m)), (0, True)) +
        Piecewise((1, And(l <= j, j <= m)), (0, True)))


def test_deltasummation_add_mul_x_kd_kd():
    assert ds(x*Kd(i, k) + Kd(j, k), (k, 1, 3)) == piecewise_fold(
        Piecewise((x, And(Integer(1) <= i, i <= 3)), (0, True)) +
        Piecewise((1, And(Integer(1) <= j, j <= 3)), (0, True)))
    assert ds(x*Kd(i, k) + Kd(j, k), (k, 1, 1)) == piecewise_fold(
        Piecewise((x, Eq(i, 1)), (0, True)) +
        Piecewise((1, Eq(j, 1)), (0, True)))
    assert ds(x*Kd(i, k) + Kd(j, k), (k, 2, 2)) == piecewise_fold(
        Piecewise((x, Eq(i, 2)), (0, True)) +
        Piecewise((1, Eq(j, 2)), (0, True)))
    assert ds(x*Kd(i, k) + Kd(j, k), (k, 3, 3)) == piecewise_fold(
        Piecewise((x, Eq(i, 3)), (0, True)) +
        Piecewise((1, Eq(j, 3)), (0, True)))
    assert ds(x*Kd(i, k) + Kd(j, k), (k, 1, l)) == piecewise_fold(
        Piecewise((x, And(Integer(1) <= i, i <= l)), (0, True)) +
        Piecewise((1, And(Integer(1) <= j, j <= l)), (0, True)))
    assert ds(x*Kd(i, k) + Kd(j, k), (k, l, 3)) == piecewise_fold(
        Piecewise((x, And(l <= i, i <= 3)), (0, True)) +
        Piecewise((1, And(l <= j, j <= 3)), (0, True)))
    assert ds(x*Kd(i, k) + Kd(j, k), (k, l, m)) == piecewise_fold(
        Piecewise((x, And(l <= i, i <= m)), (0, True)) +
        Piecewise((1, And(l <= j, j <= m)), (0, True)))


def test_deltasummation_mul_x_add_kd_kd():
    assert ds(x*(Kd(i, k) + Kd(j, k)), (k, 1, 3)) == piecewise_fold(
        Piecewise((x, And(Integer(1) <= i, i <= 3)), (0, True)) +
        Piecewise((x, And(Integer(1) <= j, j <= 3)), (0, True)))
    assert ds(x*(Kd(i, k) + Kd(j, k)), (k, 1, 1)) == piecewise_fold(
        Piecewise((x, Eq(i, 1)), (0, True)) +
        Piecewise((x, Eq(j, 1)), (0, True)))
    assert ds(x*(Kd(i, k) + Kd(j, k)), (k, 2, 2)) == piecewise_fold(
        Piecewise((x, Eq(i, 2)), (0, True)) +
        Piecewise((x, Eq(j, 2)), (0, True)))
    assert ds(x*(Kd(i, k) + Kd(j, k)), (k, 3, 3)) == piecewise_fold(
        Piecewise((x, Eq(i, 3)), (0, True)) +
        Piecewise((x, Eq(j, 3)), (0, True)))
    assert ds(x*(Kd(i, k) + Kd(j, k)), (k, 1, l)) == piecewise_fold(
        Piecewise((x, And(Integer(1) <= i, i <= l)), (0, True)) +
        Piecewise((x, And(Integer(1) <= j, j <= l)), (0, True)))
    assert ds(x*(Kd(i, k) + Kd(j, k)), (k, l, 3)) == piecewise_fold(
        Piecewise((x, And(l <= i, i <= 3)), (0, True)) +
        Piecewise((x, And(l <= j, j <= 3)), (0, True)))
    assert ds(x*(Kd(i, k) + Kd(j, k)), (k, l, m)) == piecewise_fold(
        Piecewise((x, And(l <= i, i <= m)), (0, True)) +
        Piecewise((x, And(l <= j, j <= m)), (0, True)))


def test_deltasummation_mul_add_x_y_add_kd_kd():
    assert ds((x + y)*(Kd(i, k) + Kd(j, k)), (k, 1, 3)) == piecewise_fold(
        Piecewise((x + y, And(Integer(1) <= i, i <= 3)), (0, True)) +
        Piecewise((x + y, And(Integer(1) <= j, j <= 3)), (0, True)))
    assert ds((x + y)*(Kd(i, k) + Kd(j, k)), (k, 1, 1)) == piecewise_fold(
        Piecewise((x + y, Eq(i, 1)), (0, True)) +
        Piecewise((x + y, Eq(j, 1)), (0, True)))
    assert ds((x + y)*(Kd(i, k) + Kd(j, k)), (k, 2, 2)) == piecewise_fold(
        Piecewise((x + y, Eq(i, 2)), (0, True)) +
        Piecewise((x + y, Eq(j, 2)), (0, True)))
    assert ds((x + y)*(Kd(i, k) + Kd(j, k)), (k, 3, 3)) == piecewise_fold(
        Piecewise((x + y, Eq(i, 3)), (0, True)) +
        Piecewise((x + y, Eq(j, 3)), (0, True)))
    assert ds((x + y)*(Kd(i, k) + Kd(j, k)), (k, 1, l)) == piecewise_fold(
        Piecewise((x + y, And(Integer(1) <= i, i <= l)), (0, True)) +
        Piecewise((x + y, And(Integer(1) <= j, j <= l)), (0, True)))
    assert ds((x + y)*(Kd(i, k) + Kd(j, k)), (k, l, 3)) == piecewise_fold(
        Piecewise((x + y, And(l <= i, i <= 3)), (0, True)) +
        Piecewise((x + y, And(l <= j, j <= 3)), (0, True)))
    assert ds((x + y)*(Kd(i, k) + Kd(j, k)), (k, l, m)) == piecewise_fold(
        Piecewise((x + y, And(l <= i, i <= m)), (0, True)) +
        Piecewise((x + y, And(l <= j, j <= m)), (0, True)))


def test_deltasummation_add_mul_x_y_mul_x_kd():
    assert ds(x*y + x*Kd(i, j), (j, 1, 3)) == \
        Piecewise((3*x*y + x, And(Integer(1) <= i, i <= 3)), (3*x*y, True))
    assert ds(x*y + x*Kd(i, j), (j, 1, 1)) == \
        Piecewise((x*y + x, Eq(i, 1)), (x*y, True))
    assert ds(x*y + x*Kd(i, j), (j, 2, 2)) == \
        Piecewise((x*y + x, Eq(i, 2)), (x*y, True))
    assert ds(x*y + x*Kd(i, j), (j, 3, 3)) == \
        Piecewise((x*y + x, Eq(i, 3)), (x*y, True))
    assert ds(x*y + x*Kd(i, j), (j, 1, k)) == \
        Piecewise((k*x*y + x, And(Integer(1) <= i, i <= k)), (k*x*y, True))
    assert ds(x*y + x*Kd(i, j), (j, k, 3)) == \
        Piecewise(((4 - k)*x*y + x, And(k <= i, i <= 3)), ((4 - k)*x*y, True))
    assert ds(x*y + x*Kd(i, j), (j, k, l)) == Piecewise(
        ((l - k + 1)*x*y + x, And(k <= i, i <= l)), ((l - k + 1)*x*y, True))


def test_deltasummation_mul_x_add_y_kd():
    assert ds(x*(y + Kd(i, j)), (j, 1, 3)) == \
        Piecewise((3*x*y + x, And(Integer(1) <= i, i <= 3)), (3*x*y, True))
    assert ds(x*(y + Kd(i, j)), (j, 1, 1)) == \
        Piecewise((x*y + x, Eq(i, 1)), (x*y, True))
    assert ds(x*(y + Kd(i, j)), (j, 2, 2)) == \
        Piecewise((x*y + x, Eq(i, 2)), (x*y, True))
    assert ds(x*(y + Kd(i, j)), (j, 3, 3)) == \
        Piecewise((x*y + x, Eq(i, 3)), (x*y, True))
    assert ds(x*(y + Kd(i, j)), (j, 1, k)) == \
        Piecewise((k*x*y + x, And(Integer(1) <= i, i <= k)), (k*x*y, True))
    assert ds(x*(y + Kd(i, j)), (j, k, 3)) == \
        Piecewise(((4 - k)*x*y + x, And(k <= i, i <= 3)), ((4 - k)*x*y, True))
    assert ds(x*(y + Kd(i, j)), (j, k, l)) == Piecewise(
        ((l - k + 1)*x*y + x, And(k <= i, i <= l)), ((l - k + 1)*x*y, True))


def test_deltasummation_mul_x_add_y_twokd():
    assert ds(x*(y + 2*Kd(i, j)), (j, 1, 3)) == \
        Piecewise((3*x*y + 2*x, And(Integer(1) <= i, i <= 3)), (3*x*y, True))
    assert ds(x*(y + 2*Kd(i, j)), (j, 1, 1)) == \
        Piecewise((x*y + 2*x, Eq(i, 1)), (x*y, True))
    assert ds(x*(y + 2*Kd(i, j)), (j, 2, 2)) == \
        Piecewise((x*y + 2*x, Eq(i, 2)), (x*y, True))
    assert ds(x*(y + 2*Kd(i, j)), (j, 3, 3)) == \
        Piecewise((x*y + 2*x, Eq(i, 3)), (x*y, True))
    assert ds(x*(y + 2*Kd(i, j)), (j, 1, k)) == \
        Piecewise((k*x*y + 2*x, And(Integer(1) <= i, i <= k)), (k*x*y, True))
    assert ds(x*(y + 2*Kd(i, j)), (j, k, 3)) == Piecewise(
        ((4 - k)*x*y + 2*x, And(k <= i, i <= 3)), ((4 - k)*x*y, True))
    assert ds(x*(y + 2*Kd(i, j)), (j, k, l)) == Piecewise(
        ((l - k + 1)*x*y + 2*x, And(k <= i, i <= l)), ((l - k + 1)*x*y, True))


def test_deltasummation_mul_add_x_y_add_y_kd():
    assert ds((x + y)*(y + Kd(i, j)), (j, 1, 3)) == Piecewise(
        (3*(x + y)*y + x + y, And(Integer(1) <= i, i <= 3)), (3*(x + y)*y, True))
    assert ds((x + y)*(y + Kd(i, j)), (j, 1, 1)) == \
        Piecewise(((x + y)*y + x + y, Eq(i, 1)), ((x + y)*y, True))
    assert ds((x + y)*(y + Kd(i, j)), (j, 2, 2)) == \
        Piecewise(((x + y)*y + x + y, Eq(i, 2)), ((x + y)*y, True))
    assert ds((x + y)*(y + Kd(i, j)), (j, 3, 3)) == \
        Piecewise(((x + y)*y + x + y, Eq(i, 3)), ((x + y)*y, True))
    assert ds((x + y)*(y + Kd(i, j)), (j, 1, k)) == Piecewise(
        (k*(x + y)*y + x + y, And(Integer(1) <= i, i <= k)), (k*(x + y)*y, True))
    assert ds((x + y)*(y + Kd(i, j)), (j, k, 3)) == Piecewise(
        ((4 - k)*(x + y)*y + x + y, And(k <= i, i <= 3)),
        ((4 - k)*(x + y)*y, True))
    assert ds((x + y)*(y + Kd(i, j)), (j, k, l)) == Piecewise(
        ((l - k + 1)*(x + y)*y + x + y, And(k <= i, i <= l)),
        ((l - k + 1)*(x + y)*y, True))


def test_deltasummation_mul_add_x_kd_add_y_kd():
    assert ds((x + Kd(i, k))*(y + Kd(i, j)), (j, 1, 3)) == piecewise_fold(
        Piecewise((Kd(i, k) + x, And(Integer(1) <= i, i <= 3)), (0, True)) +
        3*(Kd(i, k) + x)*y)
    assert ds((x + Kd(i, k))*(y + Kd(i, j)), (j, 1, 1)) == piecewise_fold(
        Piecewise((Kd(i, k) + x, Eq(i, 1)), (0, True)) +
        (Kd(i, k) + x)*y)
    assert ds((x + Kd(i, k))*(y + Kd(i, j)), (j, 2, 2)) == piecewise_fold(
        Piecewise((Kd(i, k) + x, Eq(i, 2)), (0, True)) +
        (Kd(i, k) + x)*y)
    assert ds((x + Kd(i, k))*(y + Kd(i, j)), (j, 3, 3)) == piecewise_fold(
        Piecewise((Kd(i, k) + x, Eq(i, 3)), (0, True)) +
        (Kd(i, k) + x)*y)
    assert ds((x + Kd(i, k))*(y + Kd(i, j)), (j, 1, k)) == piecewise_fold(
        Piecewise((Kd(i, k) + x, And(Integer(1) <= i, i <= k)), (0, True)) +
        k*(Kd(i, k) + x)*y)
    assert ds((x + Kd(i, k))*(y + Kd(i, j)), (j, k, 3)) == piecewise_fold(
        Piecewise((Kd(i, k) + x, And(k <= i, i <= 3)), (0, True)) +
        (4 - k)*(Kd(i, k) + x)*y)
    assert ds((x + Kd(i, k))*(y + Kd(i, j)), (j, k, l)) == piecewise_fold(
        Piecewise((Kd(i, k) + x, And(k <= i, i <= l)), (0, True)) +
        (l - k + 1)*(Kd(i, k) + x)*y)
