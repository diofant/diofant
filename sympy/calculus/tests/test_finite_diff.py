from sympy import S, symbols, Function
from sympy.calculus.finite_diff import (apply_finite_diff,
                                        finite_diff_weights, as_finite_diff)


def test_apply_finite_diff():
    x, h = symbols('x h')
    f = Function('f')
    assert (apply_finite_diff(1, [x-h, x+h], [f(x-h), f(x+h)], x) -
            (f(x+h)-f(x-h))/(2*h)).simplify() == 0

    assert (apply_finite_diff(1, [5, 6, 7], [f(5), f(6), f(7)], 5) -
            (-Integer(3)/2*f(5) + 2*f(6) - Integer(1)/2*f(7))).simplify() == 0


def test_finite_diff_weights():

    d = finite_diff_weights(1, [5, 6, 7], 5)
    assert d[1][2] == [-Integer(3)/2, 2, -Integer(1)/2]

    # Table 1, p. 702 in doi:10.1090/S0025-5718-1988-0935077-0
    # --------------------------------------------------------
    xl = [0, 1, -1, 2, -2, 3, -3, 4, -4]

    # d holds all coefficients
    d = finite_diff_weights(4, xl, Integer(0))

    # Zeroeth derivative
    for i in range(5):
        assert d[0][i] == [Integer(1)] + [Integer(0)]*8

    # First derivative
    assert d[1][0] == [Integer(0)]*9
    assert d[1][2] == [Integer(0), Integer(1)/2, -Integer(1)/2] + [Integer(0)]*6
    assert d[1][4] == [Integer(0), Integer(2)/3, -Integer(2)/3, -Integer(1)/12, Integer(1)/12] + [Integer(0)]*4
    assert d[1][6] == [Integer(0), Integer(3)/4, -Integer(3)/4, -Integer(3)/20, Integer(3)/20,
                       Integer(1)/60, -Integer(1)/60] + [Integer(0)]*2
    assert d[1][8] == [Integer(0), Integer(4)/5, -Integer(4)/5, -Integer(1)/5, Integer(1)/5,
                       Integer(4)/105, -Integer(4)/105, -Integer(1)/280, Integer(1)/280]

    # Second derivative
    for i in range(2):
        assert d[2][i] == [Integer(0)]*9
    assert d[2][2] == [-Integer(2), Integer(1), Integer(1)] + [Integer(0)]*6
    assert d[2][4] == [-Integer(5)/2, Integer(4)/3, Integer(4)/3, -Integer(1)/12, -Integer(1)/12] + [Integer(0)]*4
    assert d[2][6] == [-Integer(49)/18, Integer(3)/2, Integer(3)/2, -Integer(3)/20, -Integer(3)/20,
                       Integer(1)/90, Integer(1)/90] + [Integer(0)]*2
    assert d[2][8] == [-Integer(205)/72, Integer(8)/5, Integer(8)/5, -Integer(1)/5, -Integer(1)/5,
                       Integer(8)/315, Integer(8)/315, -Integer(1)/560, -Integer(1)/560]

    # Third derivative
    for i in range(3):
        assert d[3][i] == [Integer(0)]*9
    assert d[3][4] == [Integer(0), -Integer(1), Integer(1), Integer(1)/2, -Integer(1)/2] + [Integer(0)]*4
    assert d[3][6] == [Integer(0), -Integer(13)/8, Integer(13)/8, Integer(1), -Integer(1),
                       -Integer(1)/8, Integer(1)/8] + [Integer(0)]*2
    assert d[3][8] == [Integer(0), -Integer(61)/30, Integer(61)/30, Integer(169)/120, -Integer(169)/120,
                       -Integer(3)/10, Integer(3)/10, Integer(7)/240, -Integer(7)/240]

    # Fourth derivative
    for i in range(4):
        assert d[4][i] == [Integer(0)]*9
    assert d[4][4] == [Integer(6), -Integer(4), -Integer(4), Integer(1), Integer(1)] + [Integer(0)]*4
    assert d[4][6] == [Integer(28)/3, -Integer(13)/2, -Integer(13)/2, Integer(2), Integer(2),
                       -Integer(1)/6, -Integer(1)/6] + [Integer(0)]*2
    assert d[4][8] == [Integer(91)/8, -Integer(122)/15, -Integer(122)/15, Integer(169)/60, Integer(169)/60,
                       -Integer(2)/5, -Integer(2)/5, Integer(7)/240, Integer(7)/240]

    # Table 2, p. 703 in doi:10.1090/S0025-5718-1988-0935077-0
    # --------------------------------------------------------
    xl = [[j/Integer(2) for j in list(range(-i*2+1, 0, 2))+list(range(1, i*2+1, 2))]
          for i in range(1, 5)]

    # d holds all coefficients
    d = [finite_diff_weights({0: 1, 1: 2, 2: 4, 3: 4}[i], xl[i], 0) for
         i in range(4)]

    # Zeroth derivative
    assert d[0][0][1] == [Integer(1)/2, Integer(1)/2]
    assert d[1][0][3] == [-Integer(1)/16, Integer(9)/16, Integer(9)/16, -Integer(1)/16]
    assert d[2][0][5] == [Integer(3)/256, -Integer(25)/256, Integer(75)/128, Integer(75)/128,
                          -Integer(25)/256, Integer(3)/256]
    assert d[3][0][7] == [-Integer(5)/2048, Integer(49)/2048, -Integer(245)/2048, Integer(1225)/2048,
                          Integer(1225)/2048, -Integer(245)/2048, Integer(49)/2048, -Integer(5)/2048]

    # First derivative
    assert d[0][1][1] == [-Integer(1), Integer(1)]
    assert d[1][1][3] == [Integer(1)/24, -Integer(9)/8, Integer(9)/8, -Integer(1)/24]
    assert d[2][1][5] == [-Integer(3)/640, Integer(25)/384, -Integer(75)/64, Integer(75)/64,
                          -Integer(25)/384, Integer(3)/640]
    assert d[3][1][7] == [Integer(5)/7168, -Integer(49)/5120,  Integer(245)/3072, Integer(-1225)/1024,
                          Integer(1225)/1024, -Integer(245)/3072, Integer(49)/5120, -Integer(5)/7168]

    # Reasonably the rest of the table is also correct... (testing of that
    # deemed excessive at the moment)


def test_as_finite_diff():
    x, h = symbols('x h')
    f = Function('f')

    # Central 1st derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x), [x-2, x-1, x, x+1, x+2]) -
            (Integer(1)/12*(f(x-2)-f(x+2)) + Integer(2)/3*(f(x+1)-f(x-1)))).simplify() == 0

    # Central 1st derivative "half-way"
    assert (as_finite_diff(f(x).diff(x)) -
            (f(x + Integer(1)/2)-f(x - Integer(1)/2))).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), h) -
            (f(x + h/Integer(2))-f(x - h/Integer(2)))/h).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x - 3*h, x-h, x+h, x + 3*h]) -
            (Integer(9)/(8*2*h)*(f(x+h) - f(x-h)) +
             Integer(1)/(24*2*h)*(f(x - 3*h) - f(x + 3*h)))).simplify() == 0

    # One sided 1st derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x), [0, 1, 2], 0) -
            (-Integer(3)/2*f(0) + 2*f(1) - f(2)/2)).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x, x+h], x) -
            (f(x+h) - f(x))/h).simplify() == 0
    assert (as_finite_diff(f(x).diff(x), [x-h, x, x+h], x-h) -
            (-Integer(3)/(2*h)*f(x-h) + 2/h*f(x) -
             Integer(1)/(2*h)*f(x+h))).simplify() == 0

    # One sided 1st derivative "half-way"
    assert (as_finite_diff(f(x).diff(x), [x-h, x+h, x + 3*h, x + 5*h, x + 7*h])
            - 1/(2*h)*(-Integer(11)/(12)*f(x-h) + Integer(17)/(24)*f(x+h)
                       + Integer(3)/8*f(x + 3*h) - Integer(5)/24*f(x + 5*h)
                       + Integer(1)/24*f(x + 7*h))).simplify() == 0

    # Central 2nd derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x, 2), [x-h, x, x+h]) -
            h**-2 * (f(x-h) + f(x+h) - 2*f(x))).simplify() == 0

    assert (as_finite_diff(f(x).diff(x, 2), [x - 2*h, x-h, x, x+h, x + 2*h]) -
            h**-2 * (-Integer(1)/12*(f(x - 2*h) + f(x + 2*h)) +
                     Integer(4)/3*(f(x+h) + f(x-h)) - Integer(5)/2*f(x))).simplify() == 0

    # Central 2nd derivative "half-way"
    assert (as_finite_diff(f(x).diff(x, 2), [x - 3*h, x-h, x+h, x + 3*h]) -
            (2*h)**-2 * (Integer(1)/2*(f(x - 3*h) + f(x + 3*h)) -
                         Integer(1)/2*(f(x+h) + f(x-h)))).simplify() == 0

    # One sided 2nd derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x, 2), [x, x+h, x + 2*h, x + 3*h]) -
            h**-2 * (2*f(x) - 5*f(x+h) +
                     4*f(x+2*h) - f(x+3*h))).simplify() == 0

    # One sided 2nd derivative at "half-way"
    assert (as_finite_diff(f(x).diff(x, 2), [x-h, x+h, x + 3*h, x + 5*h]) -
            (2*h)**-2 * (Integer(3)/2*f(x-h) - Integer(7)/2*f(x+h) + Integer(5)/2*f(x + 3*h) -
                         Integer(1)/2*f(x + 5*h))).simplify() == 0

    # Central 3rd derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x, 3)) -
            (-f(x - 3/Integer(2)) + 3*f(x - 1/Integer(2)) -
             3*f(x + 1/Integer(2)) + f(x + 3/Integer(2)))).simplify() == 0

    assert (as_finite_diff(
        f(x).diff(x, 3), [x - 3*h, x - 2*h, x-h, x, x+h, x + 2*h, x + 3*h]) -
        h**-3 * (Integer(1)/8*(f(x - 3*h) - f(x + 3*h)) - f(x - 2*h) +
                 f(x + 2*h) + Integer(13)/8*(f(x-h) - f(x+h)))).simplify() == 0

    # Central 3rd derivative at "half-way"
    assert (as_finite_diff(f(x).diff(x, 3), [x - 3*h, x-h, x+h, x + 3*h]) -
            (2*h)**-3 * (f(x + 3*h)-f(x - 3*h) +
                         3*(f(x-h)-f(x+h)))).simplify() == 0

    # One sided 3rd derivative at gridpoint
    assert (as_finite_diff(f(x).diff(x, 3), [x, x+h, x + 2*h, x + 3*h]) -
            h**-3 * (f(x + 3*h)-f(x) + 3*(f(x+h)-f(x + 2*h)))).simplify() == 0

    # One sided 3rd derivative at "half-way"
    assert (as_finite_diff(f(x).diff(x, 3), [x-h, x+h, x + 3*h, x + 5*h]) -
            (2*h)**-3 * (f(x + 5*h)-f(x-h) +
                         3*(f(x+h)-f(x + 3*h)))).simplify() == 0
