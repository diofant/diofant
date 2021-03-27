import pytest

from diofant import Function, Integer, Rational
from diofant.abc import h, x, y
from diofant.calculus.finite_diff import (apply_finite_diff, as_finite_diff,
                                          finite_diff_weights)


__all__ = ()


def test_apply_finite_diff():
    pytest.raises(ValueError, lambda: apply_finite_diff(1, [1, 2], [3], x))

    f = Function('f')
    assert (apply_finite_diff(1, [x-h, x+h], [f(x-h), f(x+h)], x) -
            (f(x+h)-f(x-h))/(2*h)).simplify() == 0

    assert (apply_finite_diff(1, [5, 6, 7], [f(5), f(6), f(7)], 5) -
            (-3*f(5)/2 + 2*f(6) - f(7)/2)).simplify() == 0


def test_finite_diff_weights():
    pytest.raises(ValueError, lambda: finite_diff_weights(-1, [5, 6], 5))
    pytest.raises(ValueError, lambda: finite_diff_weights(Rational(1, 3),
                                                          [5, 6], 5))

    d = finite_diff_weights(1, [5, 6, 7], 5)
    assert d[1][2] == [-Rational(3, 2), 2, -Rational(1, 2)]

    # Table 1, p. 702 in doi:10.1090/S0025-5718-1988-0935077-0
    # --------------------------------------------------------
    xl = [0, 1, -1, 2, -2, 3, -3, 4, -4]

    # d holds all coefficients
    d = finite_diff_weights(4, xl, Integer(0))

    # Zeroeth derivative
    for i in range(5):
        assert d[0][i] == [1] + [0]*8

    # First derivative
    assert d[1][0] == [0]*9
    assert d[1][2] == [0, Rational(1, 2), -Rational(1, 2)] + [0]*6
    assert d[1][4] == [0, Rational(2, 3), -Rational(2, 3), -Rational(1, 12), Rational(1, 12)] + [0]*4
    assert d[1][6] == [0, Rational(3, 4), -Rational(3, 4), -Rational(3, 20), Rational(3, 20),
                       Rational(1, 60), -Rational(1, 60)] + [0]*2
    assert d[1][8] == [0, Rational(4, 5), -Rational(4, 5), -Rational(1, 5), Rational(1, 5),
                       Rational(4, 105), -Rational(4, 105), -Rational(1, 280), Rational(1, 280)]

    # Second derivative
    for i in range(2):
        assert d[2][i] == [0]*9
    assert d[2][2] == [-2, 1, 1] + [0]*6
    assert d[2][4] == [-Rational(5, 2), Rational(4, 3), Rational(4, 3),
                       -Rational(1, 12), -Rational(1, 12)] + [0]*4
    assert d[2][6] == [-Rational(49, 18), Rational(3, 2), Rational(3, 2),
                       -Rational(3, 20), -Rational(3, 20),
                       Rational(1, 90), Rational(1, 90)] + [0]*2
    assert d[2][8] == [-Rational(205, 72), Rational(8, 5), Rational(8, 5), -Rational(1, 5), -Rational(1, 5),
                       Rational(8, 315), Rational(8, 315), -Rational(1, 560), -Rational(1, 560)]

    # Third derivative
    for i in range(3):
        assert d[3][i] == [0]*9
    assert d[3][4] == [0, -1, 1, Rational(1, 2), -Rational(1, 2)] + [0]*4
    assert d[3][6] == [0, -Rational(13, 8), Rational(13, 8), 1, -1,
                       -Rational(1, 8), Rational(1, 8)] + [0]*2
    assert d[3][8] == [0, -Rational(61, 30), Rational(61, 30), Rational(169, 120), -Rational(169, 120),
                       -Rational(3, 10), Rational(3, 10), Rational(7, 240), -Rational(7, 240)]

    # Fourth derivative
    for i in range(4):
        assert d[4][i] == [0]*9
    assert d[4][4] == [6, -4, -4, 1, 1] + [0]*4
    assert d[4][6] == [Rational(28, 3), -Rational(13, 2), -Rational(13, 2), 2, 2,
                       -Rational(1, 6), -Rational(1, 6)] + [0]*2
    assert d[4][8] == [Rational(91, 8), -Rational(122, 15), -Rational(122, 15), Rational(169, 60), Rational(169, 60),
                       -Rational(2, 5), -Rational(2, 5), Rational(7, 240), Rational(7, 240)]

    # Table 2, p. 703 in doi:10.1090/S0025-5718-1988-0935077-0
    # --------------------------------------------------------
    xl = [[Rational(j, 2) for j in list(range(-i*2+1, 0, 2))+list(range(1, i*2+1, 2))]
          for i in range(1, 5)]

    # d holds all coefficients
    d = [finite_diff_weights({0: 1, 1: 2, 2: 4, 3: 4}[i], xl[i], 0) for
         i in range(4)]

    # Zeroth derivative
    assert d[0][0][1] == [Rational(1, 2), Rational(1, 2)]
    assert d[1][0][3] == [-Rational(1, 16), Rational(9, 16), Rational(9, 16), -Rational(1, 16)]
    assert d[2][0][5] == [Rational(3, 256), -Rational(25, 256), Rational(75, 128), Rational(75, 128),
                          -Rational(25, 256), Rational(3, 256)]
    assert d[3][0][7] == [-Rational(5, 2048), Rational(49, 2048), -Rational(245, 2048), Rational(1225, 2048),
                          Rational(1225, 2048), -Rational(245, 2048), Rational(49, 2048), -Rational(5, 2048)]

    # First derivative
    assert d[0][1][1] == [-1, 1]
    assert d[1][1][3] == [Rational(1, 24), -Rational(9, 8), Rational(9, 8), -Rational(1, 24)]
    assert d[2][1][5] == [-Rational(3, 640), Rational(25, 384), -Rational(75, 64), Rational(75, 64),
                          -Rational(25, 384), Rational(3, 640)]
    assert d[3][1][7] == [Rational(5, 7168), -Rational(49, 5120),  Rational(245, 3072), Rational(-1225, 1024),
                          Rational(1225, 1024), -Rational(245, 3072), Rational(49, 5120), -Rational(5, 7168)]

    # Reasonably the rest of the table is also correct... (testing of that
    # deemed excessive at the moment)


def test_as_finite_diff():
    f = Function('f')

    # Central 1st derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x), [x-2, x-1, x, x+1, x+2]) -
            (Rational(1, 12)*(f(x-2)-f(x+2)) + Rational(2, 3)*(f(x+1)-f(x-1)))).simplify() == 0

    # Central 1st derivative "half-way"
    assert (as_finite_diff(f(x).diff(x)) -
            (f(x + Rational(1, 2))-f(x - Rational(1, 2)))).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), h) -
            (f(x + h/2)-f(x - h/2))/h).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x - 3*h, x-h, x+h, x + 3*h]) -
            (9/(8*2*h)*(f(x+h) - f(x-h)) +
             1/(24*2*h)*(f(x - 3*h) - f(x + 3*h)))).simplify() == 0

    # One sided 1st derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x), [0, 1, 2], 0) -
            (-Rational(3, 2)*f(0) + 2*f(1) - f(2)/2)).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x, x+h], x) -
            (f(x+h) - f(x))/h).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x-h, x, x+h], x-h) -
            (-3/(2*h)*f(x-h) + 2/h*f(x) -
             1/(2*h)*f(x+h))).simplify() == 0

    # One sided 1st derivative "half-way"
    assert (as_finite_diff(f(x).diff(x), [x-h, x+h, x + 3*h, x + 5*h, x + 7*h])
            - 1/(2*h)*(-11*f(x-h)/12 + 17*f(x+h)/24
                       + Rational(3, 8)*f(x + 3*h) - Rational(5, 24)*f(x + 5*h)
                       + Rational(1, 24)*f(x + 7*h))).simplify() == 0

    # Central 2nd derivative at gridpoint
    assert (as_finite_diff(f(x).diff((x, 2)), [x-h, x, x+h]) -
            h**-2 * (f(x-h) + f(x+h) - 2*f(x))).simplify() == 0

    assert (as_finite_diff(f(x).diff((x, 2)), [x - 2*h, x-h, x, x+h, x + 2*h]) -
            h**-2 * (-Rational(1, 12)*(f(x - 2*h) + f(x + 2*h)) +
                     Rational(4, 3)*(f(x+h) + f(x-h)) - Rational(5, 2)*f(x))).simplify() == 0

    # Central 2nd derivative "half-way"
    assert (as_finite_diff(f(x).diff((x, 2)), [x - 3*h, x-h, x+h, x + 3*h]) -
            (2*h)**-2 * ((f(x - 3*h) + f(x + 3*h))/2 -
                         (f(x+h) + f(x-h))/2)).simplify() == 0

    # One sided 2nd derivative at gridpoint
    assert (as_finite_diff(f(x).diff((x, 2)), [x, x+h, x + 2*h, x + 3*h]) -
            h**-2 * (2*f(x) - 5*f(x+h) +
                     4*f(x+2*h) - f(x+3*h))).simplify() == 0

    # One sided 2nd derivative at "half-way"
    assert (as_finite_diff(f(x).diff((x, 2)), [x-h, x+h, x + 3*h, x + 5*h]) -
            (2*h)**-2 * (Rational(3, 2)*f(x-h) - Rational(7, 2)*f(x+h) + Rational(5, 2)*f(x + 3*h) -
                         f(x + 5*h)/2)).simplify() == 0

    # Central 3rd derivative at gridpoint
    assert (as_finite_diff(f(x).diff((x, 3))) -
            (-f(x - Rational(3, 2)) + 3*f(x - Rational(1, 2)) -
             3*f(x + Rational(1, 2)) + f(x + Rational(3, 2)))).simplify() == 0

    assert (as_finite_diff(
        f(x).diff((x, 3)), [x - 3*h, x - 2*h, x-h, x, x+h, x + 2*h, x + 3*h]) -
        h**-3 * (Rational(1, 8)*(f(x - 3*h) - f(x + 3*h)) - f(x - 2*h) +
                 f(x + 2*h) + Rational(13, 8)*(f(x-h) - f(x+h)))).simplify() == 0

    # Central 3rd derivative at "half-way"
    assert (as_finite_diff(f(x).diff((x, 3)), [x - 3*h, x-h, x+h, x + 3*h]) -
            (2*h)**-3 * (f(x + 3*h)-f(x - 3*h) +
                         3*(f(x-h)-f(x+h)))).simplify() == 0

    # One sided 3rd derivative at gridpoint
    assert (as_finite_diff(f(x).diff((x, 3)), [x, x+h, x + 2*h, x + 3*h]) -
            h**-3 * (f(x + 3*h)-f(x) + 3*(f(x+h)-f(x + 2*h)))).simplify() == 0

    # One sided 3rd derivative at "half-way"
    assert (as_finite_diff(f(x).diff((x, 3)), [x-h, x+h, x + 3*h, x + 5*h]) -
            (2*h)**-3 * (f(x + 5*h)-f(x-h) +
                         3*(f(x+h)-f(x + 3*h)))).simplify() == 0

    assert as_finite_diff(f(x).diff((x, 2))) == -2*f(x) + f(x - 1) + f(x + 1)

    d2fdxdy = f(x, y).diff(x, y)
    assert as_finite_diff(d2fdxdy, wrt=x) == (-f(x - Rational(1, 2), y) +
                                              f(x + Rational(1, 2), y))
    pytest.raises(ValueError, lambda: as_finite_diff(d2fdxdy))
    pytest.raises(ValueError, lambda: as_finite_diff(f(x).diff((x, 2)),
                                                     [x, x + h]))
