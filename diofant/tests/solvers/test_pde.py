import pytest

from diofant import (Eq, Function, I, Symbol, cos, diff, exp, log, sin,
                     symbols, tan)
from diofant.abc import a, b, c, t, x, y
from diofant.solvers.pde import (checkpdesol, classify_pde, pde_separate,
                                 pde_separate_add, pde_separate_mul, pdsolve)


__all__ = ()


def test_pde_separate():
    f, X, Y = map(Function, 'fXY')
    pytest.raises(ValueError,
                  lambda: pde_separate(f(x, y).diff(x), f(x, y),
                                       [X(x), Y(y)], strategy='spam'))

    eq = diff(f(x, y), x) - exp(f(x, y))*diff(f(x, y), y)
    assert pde_separate(eq, f(x, y),
                        [X(x), Y(y)],
                        strategy='add') == [exp(-X(x))*diff(X(x), x),
                                            exp(Y(y))*diff(Y(y), y)]


def test_pde_separate_add():
    x, y, z, t = symbols('x,y,z,t')
    F, T, X, Y, Z, u = map(Function, 'FTXYZu')

    eq = Eq(diff(u(x, t), x), diff(u(x, t), t)*exp(u(x, t)))
    res = pde_separate_add(eq, u(x, t), [X(x), T(t)])
    assert res == [diff(X(x), x)*exp(-X(x)), diff(T(t), t)*exp(T(t))]

    eq = Eq(u(x, y).diff(x), -exp(u(x, y).diff(y)))
    res = pde_separate_add(eq, u(x, y), [X(x), Y(y)])
    assert res == [diff(X(x), x), -exp(diff(Y(y), y))]


def test_pde_separate_mul():
    x, y, z, t = symbols('x,y,z,t')
    c = Symbol('C', extended_real=True)
    Phi = Function('Phi')
    F, R, T, X, Y, Z, u = map(Function, 'FRTXYZu')
    r, theta, z = symbols('r,theta,z')

    # Something simple :)
    eq = Eq(diff(F(x, y, z), x) + diff(F(x, y, z), y) + diff(F(x, y, z), z), 0)

    # Duplicate arguments in functions
    pytest.raises(ValueError, lambda: pde_separate_mul(eq, F(x, y, z), [X(x), u(z, z)]))
    # Wrong number of arguments
    pytest.raises(ValueError, lambda: pde_separate_mul(eq, F(x, y, z), [X(x), Y(y)]))
    # Wrong variables: [x, y] -> [x, z]
    pytest.raises(ValueError, lambda: pde_separate_mul(eq, F(x, y, z), [X(t), Y(x, y)]))

    assert pde_separate_mul(eq, F(x, y, z), [Y(y), u(x, z)]) == \
        [diff(Y(y), y)/Y(y), -diff(u(x, z), x)/u(x, z) - diff(u(x, z), z)/u(x, z)]
    assert pde_separate_mul(eq, F(x, y, z), [X(x), Y(y), Z(z)]) == \
        [diff(X(x), x)/X(x), -diff(Z(z), z)/Z(z) - diff(Y(y), y)/Y(y)]

    # wave equation
    wave = Eq(diff(u(x, t), t, t), c**2*diff(u(x, t), x, x))
    res = pde_separate_mul(wave, u(x, t), [X(x), T(t)])
    assert res == [diff(X(x), x, x)/X(x), diff(T(t), t, t)/(c**2*T(t))]

    # Laplace equation in cylindrical coords
    eq = Eq(1/r * diff(Phi(r, theta, z), r) + diff(Phi(r, theta, z), (r, 2)) +
            1/r**2 * diff(Phi(r, theta, z), (theta, 2)) + diff(Phi(r, theta, z), (z, 2)), 0)
    # Separate z
    res = pde_separate_mul(eq, Phi(r, theta, z), [Z(z), u(theta, r)])
    assert res == [diff(Z(z), z, z)/Z(z),
                   -diff(u(theta, r), r, r)/u(theta, r) -
                   diff(u(theta, r), r)/(r*u(theta, r)) -
                   diff(u(theta, r), theta, theta)/(r**2*u(theta, r))]
    # Lets use the result to create a new equation...
    eq = Eq(res[1], c)
    # ...and separate theta...
    res = pde_separate_mul(eq, u(theta, r), [T(theta), R(r)])
    assert res == [diff(T(theta), theta, theta)/T(theta),
                   -r*diff(R(r), r)/R(r) - r**2*diff(R(r), r, r)/R(r) - c*r**2]
    # ...or r...
    res = pde_separate_mul(eq, u(theta, r), [R(r), T(theta)])
    assert res == [r*diff(R(r), r)/R(r) + r**2*diff(R(r), r, r)/R(r) + c*r**2,
                   -diff(T(theta), theta, theta)/T(theta)]

    # Laplace eq in spherical coordinates
    r, phi, theta, C1 = symbols('r,phi,theta,C1')
    Xi = Function('Xi')
    R, Phi, Theta, u = map(Function, ['R', 'Phi', 'Theta', 'u'])
    eq = Eq(diff(Xi(r, phi, theta), (r, 2)) + 2/r * diff(Xi(r, phi, theta), r) +
            1/(r**2 * sin(phi)**2) * diff(Xi(r, phi, theta), (theta, 2)) +
            cos(phi)/(r**2 * sin(phi)) * diff(Xi(r, phi, theta), phi) +
            1/r**2 * diff(Xi(r, phi, theta), (phi, 2)), 0)
    res_theta = pde_separate(eq, Xi(r, phi, theta), [Theta(theta), u(r, phi)])
    eq_left = Eq(res_theta[1], -C1)
    res_theta = pde_separate(eq_left, u(r, phi), [Phi(phi), R(r)])
    assert (res_theta == [C1/sin(phi)**2 - diff(Phi(phi), phi, phi)/Phi(phi) -
                          diff(Phi(phi), phi)/(Phi(phi)*tan(phi)),
                          r**2*diff(R(r), r, r)/R(r) + 2*r*diff(R(r), r)/R(r)])

    # coverage tests
    assert pde_separate_mul(Eq(u(x, t).diff(x), u(x, t).diff(x, t)), u(x, t),
                            [X(x), T(t)]) == [-1, -T(t)/diff(T(t), t)]
    assert pde_separate(Eq((x + t)*u(x, t).diff(x), u(x, t).diff(t)),
                        u(x, t), [X(x), T(t)], strategy='mul') is None
    assert pde_separate(Eq(u(x, t).diff(x), u(x, t).diff(t) + t), u(x, t),
                        [X(x), T(t)], strategy='mul') is None
    assert pde_separate(Eq(u(x, t).diff(x), exp(u(x, t).diff(t))), u(x, t),
                        [X(x), T(t)], strategy='mul') is None


def test_pde_classify():
    # When more number of hints are added, add tests for classifying here.
    f, g = symbols('f g', cls=Function)
    u = f(x, y)
    eq1 = a*u + b*u.diff(x) + c*u.diff(y)
    eq2 = 3*u + 2*u.diff(x) + u.diff(y)
    eq3 = a*u + b*u.diff(x) + 2*u.diff(y)
    eq4 = x*u + u.diff(x) + 3*u.diff(y)
    eq5 = x**2*u + x*u.diff(x) + x*y*u.diff(y)
    eq6 = y*x**2*u + y*u.diff(x) + u.diff(y)
    for eq in [eq1, eq2, eq3]:
        assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    for eq in [eq4, eq5, eq6]:
        assert classify_pde(eq) == ('1st_linear_variable_coeff',)

    # coverage tests
    assert classify_pde(eq, u) == ('1st_linear_variable_coeff',)
    assert classify_pde(eq, g(x, y)) == ()
    assert classify_pde(eq, g(x, y), dict=True) == {'default': None, 'order': 0}
    assert classify_pde(u.diff(x) + I*u.diff(y) + y) == ()
    assert classify_pde(u.diff(x, y) + x) == ()
    assert classify_pde(x*u.diff(x) + u*u.diff(y) + y) == ()
    assert classify_pde(x*u.diff(x) + u*u.diff(y) + y,
                        dict=True) == {'default': None, 'ordered_hints': (),
                                       'order': 1}
    assert classify_pde(2*u.diff(x, y, y)) == ()

    eq = Eq(1 + (2*(u.diff(x)/u)) + (3*(u.diff(y)/u)), 0)
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)


def test_checkpdesol():
    f, F = map(Function, ['f', 'F'])
    eq1 = a*f(x, y) + b*f(x, y).diff(x) + c*f(x, y).diff(y)
    eq2 = 3*f(x, y) + 2*f(x, y).diff(x) + f(x, y).diff(y)
    eq3 = a*f(x, y) + b*f(x, y).diff(x) + 2*f(x, y).diff(y)
    for eq in [eq1, eq2, eq3]:
        assert checkpdesol(eq, pdsolve(eq))[0]
    eq4 = x*f(x, y) + f(x, y).diff(x) + 3*f(x, y).diff(y)
    eq5 = 2*f(x, y) + 1*f(x, y).diff(x) + 3*f(x, y).diff(y)
    eq6 = f(x, y) + 1*f(x, y).diff(x) + 3*f(x, y).diff(y)
    assert checkpdesol(eq4, [pdsolve(eq5), pdsolve(eq6)]) == [
        (False, (x - 2)*F(3*x - y)*exp(-x/5 - 3*y/5)),
        (False, (x - 1)*F(3*x - y)*exp(-x/10 - 3*y/10))]
    for eq in [eq4, eq5, eq6]:
        assert checkpdesol(eq, pdsolve(eq))[0]
    sol = pdsolve(eq4)
    sol4 = Eq(sol.lhs - sol.rhs, 0)
    pytest.raises(NotImplementedError,
                  lambda: checkpdesol(eq4, sol4, solve_for_func=False))
    assert checkpdesol(f(x, y).diff(x) - x*f(x, y).diff(y) + 1,
                       F(x**2/2 + y) - x)[0]


def test_solvefun():
    f, F, G, H = map(Function, ['f', 'F', 'G', 'H'])
    eq1 = f(x, y) + f(x, y).diff(x) + f(x, y).diff(y)
    assert pdsolve(eq1) == Eq(f(x, y), F(x - y)*exp(-x/2 - y/2))
    assert pdsolve(eq1, solvefun=G) == Eq(f(x, y), G(x - y)*exp(-x/2 - y/2))
    assert pdsolve(eq1, solvefun=H) == Eq(f(x, y), H(x - y)*exp(-x/2 - y/2))


def test_pde_1st_linear_constant_coeff_homogeneous():
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = 2*u + u.diff(x) + u.diff(y)
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(x - y)*exp(-x - y))
    assert checkpdesol(eq, sol)[0]

    eq = 4 + (3*u.diff(x)/u) + (2*u.diff(y)/u)
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(2*x - 3*y)*exp(-12*x/13 - 8*y/13))
    assert checkpdesol(eq, sol)[0]

    eq = u + (6*u.diff(x)) + (7*u.diff(y))
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(7*x - 6*y)*exp(-6*x/85 - 7*y/85))
    assert checkpdesol(eq, sol)[0]

    eq = a*u + b*u.diff(x) + c*u.diff(y)
    sol = pdsolve(eq)
    assert checkpdesol(eq, sol)[0]


def test_pde_1st_linear_constant_coeff():
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = -2*u.diff(x) + 4*u.diff(y) + 5*u - exp(x + 3*y)
    sol = pdsolve(eq)
    assert sol == Eq(f(x, y),
                     (F(4*x + 2*y) + exp(x/2 + 4*y)/15)*exp(x/2 - y))
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
                                '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = (u.diff(x)/u) + (u.diff(y)/u) + 1 - (exp(x + y)/u)
    sol = pdsolve(eq)
    assert sol == Eq(f(x, y), F(x - y)*exp(-x/2 - y/2) + exp(x + y)/3)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
                                '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = 2*u + -u.diff(x) + 3*u.diff(y) + sin(x)
    sol = pdsolve(eq)
    assert sol == Eq(f(x, y),
                     F(3*x + y)*exp(x/5 - 3*y/5) - 2*sin(x)/5 - cos(x)/5)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
                                '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = u + u.diff(x) + u.diff(y) + x*y
    sol = pdsolve(eq)
    assert sol == Eq(f(x, y),
                     -x*y + x + y + F(x - y)*exp(-x/2 - y/2) - 2)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
                                '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = u + u.diff(x) + u.diff(y) + log(x)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
                                '1st_linear_constant_coeff_Integral')


def test_pdsolve_all():
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = u + u.diff(x) + u.diff(y) + x**2*y
    sol = pdsolve(eq, hint='all')
    keys = ['1st_linear_constant_coeff',
            '1st_linear_constant_coeff_Integral', 'default', 'order']
    assert sorted(sol) == keys
    assert sol['order'] == 1
    assert sol['default'] == '1st_linear_constant_coeff'
    assert sol['1st_linear_constant_coeff'] == Eq(f(x, y),
                                                  -x**2*y + x**2 + 2*x*y - 4*x - 2*y + F(x - y)*exp(-x/2 - y/2) + 6)


def test_pdsolve_variable_coeff():
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = x*(u.diff(x)) - y*(u.diff(y)) + y**2*u - y**2
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, F(x*y)*exp(y**2/2) + 1)
    assert checkpdesol(eq, sol)[0]

    eq = x**2*u + x*u.diff(x) + x*y*u.diff(y)
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, F(y*exp(-x))*exp(-x**2/2))
    assert checkpdesol(eq, sol)[0]

    eq = y*x**2*u + y*u.diff(x) + u.diff(y)
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, F(-2*x + y**2)*exp(-x**3/3))
    assert checkpdesol(eq, sol)[0]

    eq = exp(x)**2*(u.diff(x)) + y
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, y*exp(-2*x)/2 + F(y))
    assert checkpdesol(eq, sol)[0]

    eq = exp(2*x)*(u.diff(y)) + y*u - u
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, exp((-y**2 + 2*y + 2*F(x))*exp(-2*x)/2))

    eq = u.diff(x) + y*u - 3
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, 3/y + exp(-x*y)*exp(y*F(y))/y)

    eq = u.diff(x) + x*u.diff(y) - 3
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, 3*x + F(-x**2/2 + y))

    eq = x*u.diff(y) - y
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert pdsolve(eq) == Eq(u, F(x) + y**2/x/2)


def test_sympyissue_11726():
    f, X, T = symbols('f X T', cls=Function)

    u = f(x, t)
    eq = u.diff((x, 2)) - u.diff((t, 2))
    res = pde_separate(eq, u, [T(x), X(t)])

    assert res == [diff(T(x), x, x)/T(x), diff(X(t), t, t)/X(t)]
