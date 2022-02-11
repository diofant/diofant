import itertools

import pytest

from diofant import (Add, Eq, Function, Idx, ImmutableDenseMatrix,
                     ImmutableSparseMatrix, IndexedBase, Matrix, MatrixSymbol,
                     MutableDenseMatrix, MutableSparseMatrix, O, Piecewise,
                     Pow, Rational, RootOf, Subs, Symbol, Tuple, cos, cse, exp,
                     meijerg, sin, sqrt, symbols, sympify, true)
from diofant.abc import a, b, w, x, y, z
from diofant.simplify import cse_main, cse_opts
from diofant.simplify.cse_opts import sub_post, sub_pre


__all__ = ()

x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12 = symbols('x:13')


def test_numbered_symbols():
    ns = cse_main.numbered_symbols(prefix='y')
    assert list(itertools.islice(
        ns, 0, 10)) == [Symbol(f'y{i}') for i in range(10)]
    ns = cse_main.numbered_symbols(prefix='y')
    assert list(itertools.islice(
        ns, 10, 20)) == [Symbol(f'y{i}') for i in range(10, 20)]
    ns = cse_main.numbered_symbols()
    assert list(itertools.islice(
        ns, 0, 10)) == [Symbol(f'x{i}') for i in range(10)]

# Dummy "optimization" functions for testing.


def opt1(expr):
    return expr + y


def opt2(expr):
    return expr*z


def test_preprocess_for_cse():
    assert cse_main.preprocess_for_cse(x, [(opt1, None)]) == x + y
    assert cse_main.preprocess_for_cse(x, [(None, opt1)]) == x
    assert cse_main.preprocess_for_cse(x, [(None, None)]) == x
    assert cse_main.preprocess_for_cse(x, [(opt1, opt2)]) == x + y
    assert cse_main.preprocess_for_cse(
        x, [(opt1, None), (opt2, None)]) == (x + y)*z


def test_postprocess_for_cse():
    assert cse_main.postprocess_for_cse(x, [(opt1, None)]) == x
    assert cse_main.postprocess_for_cse(x, [(None, opt1)]) == x + y
    assert cse_main.postprocess_for_cse(x, [(None, None)]) == x
    assert cse_main.postprocess_for_cse(x, [(opt1, opt2)]) == x*z
    # Note the reverse order of application.
    assert cse_main.postprocess_for_cse(
        x, [(None, opt1), (None, opt2)]) == x*z + y


def test_cse_single():
    # Simple substitution.
    e = Add(Pow(x + y, 2), sqrt(x + y))
    substs, reduced = cse([e])
    assert substs == [(x0, x + y)]
    assert reduced == [sqrt(x0) + x0**2]
    assert cse([e], order='none') == cse([e])


def test_cse_single2():
    # Simple substitution, test for being able to pass the expression directly
    e = Add(Pow(x + y, 2), sqrt(x + y))
    substs, reduced = cse(e)
    assert substs == [(x0, x + y)]
    assert reduced == [sqrt(x0) + x0**2]
    substs, reduced = cse(Matrix([[1]]))
    assert isinstance(reduced[0], Matrix)


def test_cse_not_possible():
    # No substitution possible.
    e = Add(x, y)
    substs, reduced = cse([e])
    assert not substs
    assert reduced == [x + y]
    # issue sympy/sympy#6329
    eq = (meijerg((1, 2), (y, 4), (5,), [], x) +
          meijerg((1, 3), (y, 4), (5,), [], x))
    assert cse(eq) == ([], [eq])


def test_nested_substitution():
    # Substitution within a substitution.
    e = Add(Pow(w*x + y, 2), sqrt(w*x + y))
    substs, reduced = cse([e])
    assert substs == [(x0, w*x + y)]
    assert reduced == [sqrt(x0) + x0**2]


def test_subtraction_opt():
    # Make sure subtraction is optimized.
    e = (x - y)*(z - y) + exp((x - y)*(z - y))
    substs, reduced = cse(
        [e], optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)])
    assert substs == [(x0, (x - y)*(y - z))]
    assert reduced == [-x0 + exp(-x0)]
    e = -(x - y)*(z - y) + exp(-(x - y)*(z - y))
    substs, reduced = cse(
        [e], optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)])
    assert substs == [(x0, (x - y)*(y - z))]
    assert reduced == [x0 + exp(x0)]
    # issue sympy/sympy#4077
    n = -1 + 1/x
    e = n/x/(-n)**2 - 1/n/x
    assert cse(e, optimizations=[(cse_opts.sub_pre, cse_opts.sub_post)]) == \
        ([], [0])


def test_multiple_expressions():
    e1 = (x + y)*z
    e2 = (x + y)*w
    substs, reduced = cse([e1, e2])
    assert substs == [(x0, x + y)]
    assert reduced == [x0*z, x0*w]
    l = [w*x*y + z, w*y]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    assert substs == rsubsts
    assert reduced == [z + x*x0, x0]
    l = [w*x*y, w*x*y + z, w*y]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    assert substs == rsubsts
    assert reduced == [x1, x1 + z, x0]
    l = [(x - z)*(y - z), x - z, y - z]
    substs, reduced = cse(l)
    rsubsts, _ = cse(reversed(l))
    assert substs == [(x0, -z), (x1, x + x0), (x2, x0 + y)]
    assert rsubsts == [(x0, -z), (x1, x0 + y), (x2, x + x0)]
    assert reduced == [x1*x2, x1, x2]
    l = [w*y + w + x + y + z, w*x*y]
    assert cse(l) == ([(x0, w*y)], [w + x + x0 + y + z, x*x0])
    assert cse([x + y, x + y + z]) == ([(x0, x + y)], [x0, z + x0])
    assert cse([x + y, x + z]) == ([], [x + y, x + z])
    assert cse([x*y, z + x*y, x*y*z + 3]) == \
        ([(x0, x*y)], [x0, z + x0, 3 + x0*z])


def test_non_commutative_cse():
    A, B, C = symbols('A B C', commutative=False)
    l = [A*B*C, A*C]
    assert cse(l) == ([], l)


@pytest.mark.xfail
def test_non_commutative_cse_mul():
    x0 = symbols('x0', commutative=False)
    A, B, C = symbols('A B C', commutative=False)
    l = [A*B*C, A*B]
    assert cse(l) == ([(x0, A*B)], [x0*C, x0])


# Test if CSE of non-commutative Mul terms is disabled
def test_bypass_non_commutatives():
    A, B, C = symbols('A B C', commutative=False)
    l = [A*B*C, A*C]
    assert cse(l) == ([], l)
    l = [A*B*C, A*B]
    assert cse(l) == ([], l)
    l = [B*C, A*B*C]
    assert cse(l) == ([], l)


def test_non_commutative_order():
    A, B, C = symbols('A B C', commutative=False)
    x0 = symbols('x0', commutative=False)
    l = [B+C, A*(B+C)]
    assert cse(l) == ([(x0, B+C)], [x0, A*x0])
    l = [(A - B)**2 + A - B]
    assert cse(l) == ([(x0, A - B)], [x0**2 + x0])


@pytest.mark.xfail
def test_powers():
    assert cse(x*y**2 + x*y) == ([(x0, x*y)], [x0*y + x0])


def test_basic_optimization():
    # issue sympy/sympy#4498
    assert cse(w/(x - y) + z/(y - x), optimizations='basic') == \
        ([], [(w - z)/(x - y)])


def test_sympyissue_4020():
    assert cse(x**5 + x**4 + x**3 + x**2, optimizations='basic') \
        == ([(x0, x**2)], [x0*(x**3 + x + x0 + 1)])


def test_sympyissue_4203():
    assert cse(sin(x**x)/x**x) == ([(x0, x**x)], [sin(x0)/x0])


def test_sympyissue_6263():
    e = Eq(x*(-x + 1) + x*(x - 1), 0)
    assert cse(e, optimizations='basic') == ([], [True])


def test_dont_cse_tuples():
    f = Function('f')
    g = Function('g')

    name_val, (expr,) = cse(Subs(f(x, y), (x, 0), (y, 1)) +
                            Subs(g(x, y), (x, 0), (y, 1)))

    assert not name_val
    assert expr == (Subs(f(x, y), (x, 0), (y, 1))
                    + Subs(g(x, y), (x, 0), (y, 1)))

    name_val, (expr,) = cse(Subs(f(x, y), (x, 0), (y, x + y)) +
                            Subs(g(x, y), (x, 0), (y, x + y)))

    assert name_val == [(x0, x + y)]
    assert expr == Subs(f(x, y), (x, 0), (y, x0)) + Subs(g(x, y), (x, 0), (y, x0))


def test_pow_invpow():
    assert cse(1/x**2 + x**2) == \
        ([(x0, x**2)], [x0 + 1/x0])
    assert cse(x**2 + (1 + 1/x**2)/x**2) == \
        ([(x0, x**2), (x1, 1/x0)], [x0 + x1*(x1 + 1)])
    assert cse(1/x**2 + (1 + 1/x**2)*x**2) == \
        ([(x0, x**2), (x1, 1/x0)], [x0*(x1 + 1) + x1])
    assert cse(cos(1/x**2) + sin(1/x**2)) == \
        ([(x0, x**(-2))], [sin(x0) + cos(x0)])
    assert cse(cos(x**2) + sin(x**2)) == \
        ([(x0, x**2)], [sin(x0) + cos(x0)])
    assert cse(y/(2 + x**2) + z/x**2/y) == \
        ([(x0, x**2)], [y/(x0 + 2) + z/(x0*y)])
    assert cse(exp(x**2) + x**2*cos(1/x**2)) == \
        ([(x0, x**2)], [x0*cos(1/x0) + exp(x0)])
    assert cse((1 + 1/x**2)/x**2) == \
        ([(x0, x**(-2))], [x0*(x0 + 1)])
    assert cse(x**(2*y) + x**(-2*y)) == \
        ([(x0, x**(2*y))], [x0 + 1/x0])


def test_postprocess():
    eq = (x + 1 + exp((x + 1)/(y + 1)) + cos(y + 1))
    assert cse([eq, Eq(x, z + 1), z - 2, (z + 1)*(x + 1)],
               postprocess=cse_main.cse_separate) == \
        [[(x1, y + 1), (x2, z + 1), (x, x2), (x0, x + 1)],
         [x0 + exp(x0/x1) + cos(x1), z - 2, x0*x2]]


def test_sympyissue_4499():
    # previously, this gave 16 constants
    B = Function('B')
    G = Function('G')
    t = Tuple(*
              (a, a + Rational(1, 2), 2*a, b, 2*a - b + 1, (sqrt(z)/2)**(-2*a + 1)*B(2*a -
                                                                                     b, sqrt(z))*B(b - 1, sqrt(z))*G(b)*G(2*a - b + 1),
               sqrt(z)*(sqrt(z)/2)**(-2*a + 1)*B(b, sqrt(z))*B(2*a - b,
                                                               sqrt(z))*G(b)*G(2*a - b + 1), sqrt(z)*(sqrt(z)/2)**(-2*a + 1)*B(b - 1,
                                                                                                                               sqrt(z))*B(2*a - b + 1, sqrt(z))*G(b)*G(2*a - b + 1),
               (sqrt(z)/2)**(-2*a + 1)*B(b, sqrt(z))*B(2*a - b + 1,
                                                       sqrt(z))*G(b)*G(2*a - b + 1), 1, 0, Rational(1, 2), z/2, -b + 1, -2*a + b,
                  -2*a))
    c = cse(t)
    ans = (
        [(x0, 2*a), (x1, -b), (x2, x1 + 1), (x3, x0 + x2), (x4, sqrt(z)), (x5,
                                                                           B(x0 + x1, x4)), (x6, G(b)), (x7, G(x3)), (x8, -x0), (x9,
                                                                                                                                 (x4/2)**(x8 + 1)), (x10, x6*x7*x9*B(b - 1, x4)), (x11, x6*x7*x9*B(b,
                                                                                                                                                                                                   x4)), (x12, B(x3, x4))], [(a, a + Rational(1, 2), x0, b, x3, x10*x5,
                                                                                                                                                                                                                              x11*x4*x5, x10*x12*x4, x11*x12, 1, 0, Rational(1, 2), z/2, x2, b + x8, x8)])
    assert ans == c


def test_sympyissue_6169():
    r = RootOf(x**6 - 4*x**5 - 2, 1)
    assert cse(r) == ([], [r])
    # and a check that the right thing is done with the new
    # mechanism
    assert sub_post(sub_pre((-x - y)*z - x - y)) == -z*(x + y) - x - y


def test_cse_Indexed():
    len_y = 5
    y = IndexedBase('y', shape=(len_y,))
    x = IndexedBase('x', shape=(len_y,))
    i = Idx('i', len_y-1)

    expr1 = (y[i+1]-y[i])/(x[i+1]-x[i])
    expr2 = 1/(x[i+1]-x[i])
    replacements, _ = cse([expr1, expr2])
    assert len(replacements) > 0


@pytest.mark.xfail
def test_cse_MatrixSymbol():
    A = MatrixSymbol('A', 3, 3)
    y = MatrixSymbol('y', 3, 1)

    expr1 = (A.T*A).inverse() * A * y
    expr2 = (A.T*A) * A * y
    replacements, _ = cse([expr1, expr2])
    assert len(replacements) > 0


def test_Piecewise():
    f = Piecewise((-z + x*y, Eq(y, 0)), (-z - x*y, True))
    ans = cse(f)
    actual_ans = ([(x0, -z), (x1, x*y)], [Piecewise((x0+x1, Eq(y, 0)), (x0 - x1, True))])
    assert ans == actual_ans


def test_ignore_order_terms():
    eq = exp(x).series(x, 0, 3) + sin(y + x**3) - 1
    assert cse(eq) == ([], [sin(x**3 + y) + x + x**2/2 + O(x**3)])


def test_name_conflict():
    z1 = x0 + y
    z2 = x2 + x3
    l = [cos(z1) + z1, cos(z2) + z2, x0 + x2]
    substs, reduced = cse(l)
    assert [e.subs(reversed(substs)) for e in reduced] == l


def test_name_conflict_cust_symbols():
    z1 = x0 + y
    z2 = x2 + x3
    l = [cos(z1) + z1, cos(z2) + z2, x0 + x2]
    substs, reduced = cse(l, symbols('x:10'))
    assert [e.subs(reversed(substs)) for e in reduced] == l


def test_symbols_exhausted_error():
    l = cos(x+y)+x+y+cos(w+y)+sin(w+y)
    sym = [x, y, z]
    with pytest.raises(ValueError):
        cse(l, symbols=sym)


def test_sympyissue_7840():
    # daveknippers' example
    C393 = sympify(
        'Piecewise((C391 - 1.65, C390 < 0.5), (Piecewise((C391 - 1.65, \
        C391 > 2.35), (C392, True)), True))'
    )
    C391 = sympify(
        'Piecewise((2.05*C390**(-1.03), C390 < 0.5), (2.5*C390**(-0.625), True))'
    )
    C393 = C393.subs({'C391': C391})
    # simple substitution
    sub = {}
    sub['C390'] = 0.703451854
    sub['C392'] = 1.01417794
    ss_answer = C393.subs(sub)
    # cse
    substitutions, new_eqn = cse(C393)
    for pair in substitutions:
        sub[pair[0].name] = pair[1].subs(sub)
    cse_answer = new_eqn[0].subs(sub)
    # both methods should be the same
    assert ss_answer == cse_answer

    # GitRay's example
    expr = Piecewise((Symbol('ON'), Eq(Symbol('mode'), Symbol('ON'))),
                     (Piecewise((Piecewise((Symbol('OFF'), Symbol('x') < Symbol('threshold')),
                                           (Symbol('ON'), true)), Eq(Symbol('mode'), Symbol('AUTO'))),
                                (Symbol('OFF'), true)), true))
    substitutions, new_eqn = cse(expr)
    # this Piecewise should be exactly the same
    assert new_eqn[0] == expr
    # there should not be any replacements
    assert len(substitutions) < 1


def test_matrices():
    # issue sympy/sympy#8891
    for cls in (MutableDenseMatrix, MutableSparseMatrix,
                ImmutableDenseMatrix, ImmutableSparseMatrix):
        m = cls(2, 2, [x + y, 0, 0, 0])
        res = cse([x + y, m])
        ans = ([(x0, x + y)], [x0, cls([[x0, 0], [0, 0]])])
        assert res == ans
        assert isinstance(res[1][-1], cls)


def test_cse_ignore():
    exprs = [exp(y)*(3*y + 3*sqrt(x+1)), exp(y)*(5*y + 5*sqrt(x+1))]
    subst1, _ = cse(exprs)
    assert any(y in sub.free_symbols for _, sub in subst1), 'cse failed to identify any term with y'

    subst2, _ = cse(exprs, ignore=(y,))  # y is not allowed in substitutions
    assert not any(y in sub.free_symbols for _, sub in subst2), 'Sub-expressions containing y must be ignored'
    assert any(sub - sqrt(x + 1) == 0 for _, sub in subst2), 'cse failed to identify sqrt(x + 1) as sub-expression'
