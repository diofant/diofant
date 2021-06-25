import pytest

from diofant import (Derivative, Function, Integer, Matrix, Subs, Symbol,
                     atan2, cos, exp, simplify, sin, sqrt, sstr, symbols,
                     trigsimp)
from diofant.abc import t, x, y
from diofant.diffgeom import (BaseCovarDerivativeOp, Commutator, CoordSystem,
                              CovarDerivativeOp, Differential, LieDerivative,
                              Manifold, Patch, TensorProduct, WedgeProduct,
                              contravariant_order, covariant_order,
                              intcurve_diffequ, intcurve_series,
                              metric_to_Christoffel_1st,
                              metric_to_Christoffel_2nd,
                              metric_to_Ricci_components,
                              metric_to_Riemann_components, twoform_to_matrix)
from diofant.diffgeom.rn import R2, R2_p, R2_r, R3_c, R3_r, R3_s


__all__ = ()

TP = TensorProduct


def test_R2():
    x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0', extended_real=True)
    point_r = R2_r.point([x0, y0])
    point_p = R2_p.point([r0, theta0])

    # r**2 = x**2 + y**2
    assert (R2.r**2 - R2.x**2 - R2.y**2).rcall(point_r) == 0
    assert trigsimp( (R2.r**2 - R2.x**2 - R2.y**2).rcall(point_p) ) == 0
    assert trigsimp(R2.e_r(R2.x**2 + R2.y**2).rcall(point_p).doit()) == 2*r0

    # polar->rect->polar == Id
    a, b = symbols('a b', positive=True)
    m = Matrix([[a], [b]])
    # TODO assert m == R2_r.coord_tuple_transform_to(R2_p, R2_p.coord_tuple_transform_to(R2_r, [a, b])).applyfunc(simplify)
    assert m == R2_p.coord_tuple_transform_to(
        R2_r, R2_r.coord_tuple_transform_to(R2_p, m)).applyfunc(simplify)


def test_R3():
    a, b, c = symbols('a b c', positive=True)
    m = Matrix([[a], [b], [c]])
    assert m == R3_c.coord_tuple_transform_to(
        R3_r, R3_r.coord_tuple_transform_to(R3_c, m)).applyfunc(simplify)
    # TODO assert m == R3_r.coord_tuple_transform_to(R3_c, R3_c.coord_tuple_transform_to(R3_r, m)).applyfunc(simplify)
    assert m == R3_s.coord_tuple_transform_to(
        R3_r, R3_r.coord_tuple_transform_to(R3_s, m)).applyfunc(simplify)
    # TODO assert m == R3_r.coord_tuple_transform_to(R3_s, R3_s.coord_tuple_transform_to(R3_r, m)).applyfunc(simplify)
    assert m == R3_s.coord_tuple_transform_to(
        R3_c, R3_c.coord_tuple_transform_to(R3_s, m)).applyfunc(simplify)
    # TODO assert m == R3_c.coord_tuple_transform_to(R3_s, R3_s.coord_tuple_transform_to(R3_c, m)).applyfunc(simplify)


def test_point():
    p = R2_r.point([x, y])
    # TODO assert p.free_symbols() == {x, y}
    assert p.coords(R2_r) == p.coords() == Matrix([x, y])
    assert p.coords(R2_p) == Matrix([sqrt(x**2 + y**2), atan2(y, x)])


def test_coords():
    r, theta = symbols('r, theta')
    m = Manifold('M', 2)
    patch = Patch('P', m)
    rect = CoordSystem('rect', patch)
    polar = CoordSystem('polar', patch)
    polar.connect_to(rect, [r, theta], [r*cos(theta), r*sin(theta)])
    polar.coord_tuple_transform_to(rect, [0, 2]) == Matrix([[0], [0]])


def test_commutator():
    assert Commutator(R2.e_x, R2.e_y) == 0
    assert Commutator(R2.x*R2.e_x, R2.x*R2.e_x) == 0
    assert Commutator(R2.x*R2.e_x, R2.x*R2.e_y) == R2.x*R2.e_y
    c = Commutator(R2.e_x, R2.e_r)
    assert c(R2.x) == R2.y*(R2.x**2 + R2.y**2)**(-1)*sin(R2.theta)


def test_differential():
    xdy = R2.x*R2.dy
    dxdy = Differential(xdy)
    assert xdy.rcall(None) == xdy
    assert dxdy(R2.e_x, R2.e_y) == 1
    assert dxdy(R2.e_x, R2.x*R2.e_y) == R2.x
    assert Differential(dxdy) == 0


def test_products():
    assert TensorProduct(
        R2.dx, R2.dy)(R2.e_x, R2.e_y) == R2.dx(R2.e_x)*R2.dy(R2.e_y) == 1
    assert WedgeProduct(R2.dx, R2.dy)(R2.e_x, R2.e_y) == 1
    assert TensorProduct(R2.dx, R2.dy)(None, R2.e_y) == R2.dx
    assert TensorProduct(R2.dx, R2.dy)(R2.e_x, None) == R2.dy
    assert TensorProduct(R2.dx, R2.dy)(R2.e_x) == R2.dy
    assert TensorProduct(R2.x, R2.dx) == R2.x*R2.dx


def test_lie_derivative():
    assert LieDerivative(R2.e_x, R2.y) == R2.e_x(R2.y) == 0
    assert LieDerivative(R2.e_x, R2.x) == R2.e_x(R2.x) == 1
    assert LieDerivative(R2.e_x, R2.e_x) == Commutator(R2.e_x, R2.e_x) == 0
    assert LieDerivative(R2.e_x, R2.e_r) == Commutator(R2.e_x, R2.e_r)
    assert LieDerivative(R2.e_x + R2.e_y, R2.x) == 1
    assert LieDerivative(
        R2.e_x, TensorProduct(R2.dx, R2.dy))(R2.e_x, R2.e_y) == 0


def test_covar_deriv():
    ch = metric_to_Christoffel_2nd(TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy))
    cvd = BaseCovarDerivativeOp(R2_r, 0, ch)
    assert cvd(R2.x) == 1
    assert cvd(R2.x*R2.e_x) == R2.e_x
    cvd = CovarDerivativeOp(R2.x*R2.e_x, ch)
    assert cvd(R2.x) == R2.x
    assert cvd(R2.x*R2.e_x) == R2.x*R2.e_x


def test_intcurve_diffequ():
    start_point = R2_r.point([1, 0])
    vector_field = -R2.y*R2.e_x + R2.x*R2.e_y
    equations, init_cond = intcurve_diffequ(vector_field, t, start_point)
    assert sstr(equations) == '[f_1(t) + Derivative(f_0(t), t), -f_0(t) + Derivative(f_1(t), t)]'
    assert sstr(init_cond) == '[f_0(0) - 1, f_1(0)]'
    equations, init_cond = intcurve_diffequ(vector_field, t, start_point, R2_p)
    assert sstr(
        equations) == '[Derivative(f_0(t), t), Derivative(f_1(t), t) - 1]'
    assert sstr(init_cond) == '[f_0(0) - 1, f_1(0)]'

    start_point = R2_r.point([x, y])
    vector_field = R2_r.e_x
    assert intcurve_series(vector_field, t, start_point,
                           n=3) == Matrix([[t + x], [y]])


def test_helpers_and_coordinate_dependent():
    one_form = R2.dr + R2.dx
    two_form = Differential(R2.x*R2.dr + R2.r*R2.dx)
    three_form = Differential(
        R2.y*two_form) + Differential(R2.x*Differential(R2.r*R2.dr))
    metric = TensorProduct(R2.dx, R2.dx) + TensorProduct(R2.dy, R2.dy)
    metric_ambig = TensorProduct(R2.dx, R2.dx) + TensorProduct(R2.dr, R2.dr)
    misform_a = TensorProduct(R2.dr, R2.dr) + R2.dr
    misform_b = R2.dr**4
    misform_c = R2.dx*R2.dy
    twoform_not_sym = TensorProduct(R2.dx, R2.dx) + TensorProduct(R2.dx, R2.dy)
    twoform_not_TP = WedgeProduct(R2.dx, R2.dy)

    assert covariant_order(one_form) == 1
    assert covariant_order(two_form) == 2
    assert covariant_order(three_form) == 3
    assert covariant_order(two_form + metric) == 2
    assert covariant_order(two_form + metric_ambig) == 2
    assert covariant_order(two_form + twoform_not_sym) == 2
    assert covariant_order(two_form + twoform_not_TP) == 2

    pytest.raises(ValueError, lambda: covariant_order(misform_a))
    pytest.raises(ValueError, lambda: covariant_order(misform_b))
    pytest.raises(ValueError, lambda: covariant_order(misform_c))

    assert twoform_to_matrix(metric) == Matrix([[1, 0], [0, 1]])
    assert twoform_to_matrix(twoform_not_sym) == Matrix([[1, 0], [1, 0]])
    assert twoform_to_matrix(twoform_not_TP) == Matrix([[0, -1], [1, 0]])

    pytest.raises(ValueError, lambda: twoform_to_matrix(one_form))
    pytest.raises(ValueError, lambda: twoform_to_matrix(three_form))
    pytest.raises(ValueError, lambda: twoform_to_matrix(metric_ambig))

    pytest.raises(ValueError, lambda: metric_to_Christoffel_1st(twoform_not_sym))
    pytest.raises(ValueError, lambda: metric_to_Christoffel_2nd(twoform_not_sym))
    pytest.raises(ValueError, lambda: metric_to_Riemann_components(twoform_not_sym))
    pytest.raises(ValueError, lambda: metric_to_Ricci_components(twoform_not_sym))


def test_correct_arguments():
    pytest.raises(ValueError, lambda: R2.e_x(R2.e_x))
    pytest.raises(ValueError, lambda: R2.e_x(R2.dx))

    pytest.raises(ValueError, lambda: Commutator(R2.e_x, R2.x))
    pytest.raises(ValueError, lambda: Commutator(R2.dx, R2.e_x))

    pytest.raises(ValueError, lambda: Differential(Differential(R2.e_x)))

    pytest.raises(ValueError, lambda: R2.dx(R2.x))

    pytest.raises(ValueError, lambda: TensorProduct(R2.e_x, R2.dx))

    pytest.raises(ValueError, lambda: LieDerivative(R2.dx, R2.dx))
    pytest.raises(ValueError, lambda: LieDerivative(R2.x, R2.dx))

    pytest.raises(ValueError, lambda: CovarDerivativeOp(R2.dx, []))
    pytest.raises(ValueError, lambda: CovarDerivativeOp(R2.x, []))

    a = Symbol('a')
    pytest.raises(ValueError, lambda: intcurve_series(R2.dx, a, R2_r.point([1, 2])))
    pytest.raises(ValueError, lambda: intcurve_series(R2.x, a, R2_r.point([1, 2])))

    pytest.raises(ValueError, lambda: intcurve_diffequ(R2.dx, a, R2_r.point([1, 2])))
    pytest.raises(ValueError, lambda: intcurve_diffequ(R2.x, a, R2_r.point([1, 2])))

    pytest.raises(ValueError, lambda: contravariant_order(R2.e_x + R2.dx))
    pytest.raises(ValueError, lambda: contravariant_order(R2.dx**2))
    pytest.raises(ValueError, lambda: covariant_order(R2.e_x + R2.dx))

    pytest.raises(ValueError, lambda: contravariant_order(R2.e_x*R2.e_y))
    pytest.raises(ValueError, lambda: covariant_order(R2.dx*R2.dy))

    assert covariant_order(Integer(0), True) == -1
    assert contravariant_order(Integer(0), True) == -1


def test_simplify():
    x, y = R2_r.coord_functions()
    dx, dy = R2_r.base_oneforms()
    ex, ey = R2_r.base_vectors()
    assert simplify(x) == x
    assert simplify(x*y) == x*y
    assert simplify(dx*dy) == dx*dy
    assert simplify(ex*ey) == ex*ey
    assert ((1-x)*dx)/(1-x)**2 == dx/(1-x)


def test_schwarzschild():
    m = Manifold('Schwarzschild', 4)
    p = Patch('origin', m)
    cs = CoordSystem('spherical', p, ['t', 'r', 'theta', 'phi'])
    t, r, theta, phi = cs.coord_functions()
    dt, dr, dtheta, dphi = cs.base_oneforms()
    f, g = symbols('f g', cls=Function)
    metric = (exp(2*f(r))*TP(dt, dt) - exp(2*g(r))*TP(dr, dr) -
              r**2*TP(dtheta, dtheta) - r**2*sin(theta)**2*TP(dphi, dphi))
    ricci = metric_to_Ricci_components(metric)
    assert all(ricci[i, j] == 0 for i in range(4) for j in range(4) if i != j)
    R = Symbol('R')
    eq1 = simplify((ricci[0, 0]/exp(2*f(r) - 2*g(r)) +
                    ricci[1, 1])*r/2).subs({r: R}).doit()
    assert eq1 == f(R).diff(R) + g(R).diff(R)
    eq2 = simplify(ricci[1, 1].replace(g, lambda x: -f(x)).replace(r, R).doit())
    assert eq2 == -2*f(R).diff(R)**2 - f(R).diff((R, 2)) - 2*f(R).diff(R)/R


def test_sympyissue_11799():
    n = 2
    M = Manifold('M', n)
    P = Patch('P', M)

    coord = CoordSystem('coord', P, [f'x{i}' for i in range(n)])
    x = coord.coord_functions()
    dx = coord.base_oneforms()

    f = Function('f')
    g = [[f(x[0], x[1])**2, 0], [0, f(x[0], x[1])**2]]
    metric = sum(g[i][j]*TP(dx[i], dx[j]) for i in range(n) for j in range(n))

    R = metric_to_Riemann_components(metric)
    d = Symbol('d')

    assert (R[0, 1, 0, 1] ==
            -Subs(Derivative(f(d, x[1]), d, d), (d, x[0]))/f(x[0], x[1]) -
            Subs(Derivative(f(x[0], d), d, d), (d, x[1]))/f(x[0], x[1]) +
            Subs(Derivative(f(d, x[1]), d), (d, x[0]))**2/f(x[0], x[1])**2 +
            Subs(Derivative(f(x[0], d), d), (d, x[1]))**2/f(x[0], x[1])**2)
