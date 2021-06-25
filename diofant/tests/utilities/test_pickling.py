import copy
import pickle
import warnings

import pytest

from diofant import (QQ, ZZ, Abs, Add, Atom, Basic, Catalan, CoercionFailed,
                     Derivative, DiracDelta, DomainError, Dummy, E, Eijk,
                     Equality, EulerGamma, EvaluationFailed, ExpressionDomain,
                     ExtraneousFactors, FlagError, Float, FractionField,
                     Function, FunctionClass, GeneratorsError,
                     GeneratorsNeeded, GoldenRatio, GreaterThan, GroebnerBasis,
                     Heaviside, HeuristicGCDFailed, HomomorphismFailed, I,
                     Integer, Integral, Interval, IsomorphismFailed, Lambda,
                     LambertW, LessThan, Limit, Matrix, Monomial, Mul,
                     MultivariatePolynomialError, NotAlgebraic, NotInvertible,
                     NotReversible, OptionError, Options, Order, Piecewise,
                     Poly, PolynomialError, PolynomialRing, Pow, Product,
                     PurePoly, PythonRational, Rational, RefinementFailed,
                     Relational, RootOf, RootSum, Sieve, SparseMatrix,
                     StrictGreaterThan, StrictLessThan, Sum, Symbol,
                     Unequality, UnificationFailed, UnivariatePolynomialError,
                     Wild, WildFunction, acos, acosh, acot, acoth, arg, asin,
                     asinh, assoc_legendre, atan, atan2, atanh, bell,
                     bernoulli, binomial, ceiling, chebyshevt, chebyshevt_root,
                     chebyshevu, chebyshevu_root, conjugate, cos, cosh, cot,
                     coth, dirichlet_eta, erf, exp, factorial, ff, fibonacci,
                     floor, gamma, harmonic, hermite, im, legendre, ln, log,
                     loggamma, lowergamma, lucas, nan, oo, pi, polygamma, re,
                     rf, sign, sin, sinh, sqrt, tan, tanh, uppergamma,
                     vectorize, zeta, zoo)
from diofant.abc import x, y, z
from diofant.core.compatibility import HAS_GMPY
from diofant.core.logic import Logic
from diofant.core.singleton import S, SingletonRegistry
from diofant.domains import AlgebraicField, ComplexField, RealField
from diofant.domains.finitefield import GMPYFiniteField, PythonFiniteField
from diofant.domains.integerring import GMPYIntegerRing, PythonIntegerRing
from diofant.domains.rationalfield import (GMPYRationalField,
                                           PythonRationalField)
from diofant.geometry import (Circle, Ellipse, Line, Point, Polygon, Ray,
                              RegularPolygon, Segment, Triangle)
from diofant.geometry.entity import GeometryEntity
from diofant.geometry.line import LinearEntity
from diofant.plotting.plot import Plot
from diofant.polys.orderings import (GradedLexOrder, InverseOrder, LexOrder,
                                     ProductOrder, ReversedGradedLexOrder)
from diofant.printing.latex import LatexPrinter
from diofant.printing.mathml import MathMLPrinter
from diofant.printing.pretty.pretty import PrettyPrinter
from diofant.printing.pretty.stringpict import prettyForm, stringPict
from diofant.printing.printer import Printer
from diofant.printing.python import PythonPrinter
from diofant.utilities.exceptions import DiofantDeprecationWarning


__all__ = ()

# XXX we need this to initialize singleton registry properly
half = Rational(1, 2)
ø = S.EmptySet
ℕ = S.Naturals0
Id = Lambda(x, x)

excluded_attrs = {'_assumptions', '_hash', '__dict__'}


def check(a, exclude=[], check_attr=True):
    """Check that pickling and copying round-trips."""
    # Python 2.6+ warns about BasicException.message, for example.
    warnings.filterwarnings('ignore', category=DeprecationWarning)

    protocols = list(range(5)) + [copy.copy, copy.deepcopy]
    for protocol in protocols:
        if protocol in exclude:
            continue

        if callable(protocol):
            if isinstance(a, type):
                # Classes can't be copied, but that's okay.
                return
            b = protocol(a)
        else:
            b = pickle.loads(pickle.dumps(a, protocol))

        d1 = dir(a)
        d2 = dir(b)
        assert set(d1) == set(d2)

        if not check_attr:
            continue

        def c(a, b, d):
            for i in d:
                if i in excluded_attrs:
                    continue
                try:
                    if not hasattr(a, i):
                        continue
                except NotImplementedError:
                    continue
                attr = getattr(a, i)
                if not hasattr(attr, '__call__'):
                    assert hasattr(b, i), i
                    assert getattr(b, i) == attr, f'{getattr(b, i)} != {attr}'
        c(a, b, d1)
        c(b, a, d2)

    # reset filters
    warnings.simplefilter('default', category=DeprecationWarning)
    warnings.simplefilter('error', category=DiofantDeprecationWarning)

# ================= core =========================


def test_core_basic():
    for c in (Atom, Atom(), Basic, Basic(),
              SingletonRegistry, SingletonRegistry()):
        check(c)


def test_core_symbol():
    # make the Symbol a unique name that doesn't class with any other
    # testing variable in this file since after this test the symbol
    # having the same name will be cached as noncommutative
    for c in (Dummy, Dummy('x', commutative=False), Symbol,
              Symbol('_sympyissue_6229', commutative=False),
              Wild, Wild('x')):
        check(c)


def test_core_numbers():
    for c in (Integer(2), Rational(2, 3), Float('1.2')):
        check(c)


def test_sympyissue_19988():
    c = pi.evalf(100)

    check(c)


def test_core_relational():
    for c in (Equality, Equality(x, y), GreaterThan, GreaterThan(x, y),
              LessThan, LessThan(x, y), Relational, Relational(x, y),
              StrictGreaterThan, StrictGreaterThan(x, y), StrictLessThan,
              StrictLessThan(x, y), Unequality, Unequality(x, y)):
        check(c)


def test_core_add():
    for c in (Add, Add(x, 4)):
        check(c)


def test_core_mul():
    for c in (Mul, Mul(x, 4)):
        check(c)


def test_core_power():
    for c in (Pow, Pow(x, 4)):
        check(c)


def test_core_function():
    for f in (Derivative, Derivative(x), Function, FunctionClass, Lambda,
              WildFunction):
        check(f)


@pytest.mark.xfail
def test_core_dynamicfunctions():
    f = Function('f')
    check(f)


def test_core_interval():
    for c in (Interval, Interval(0, 2)):
        check(c)


def test_core_multidimensional():
    for c in (vectorize, vectorize(0)):
        check(c)


def test_Singletons():
    protocols = [0, 1, 2, 3]
    copiers = [copy.copy, copy.deepcopy]
    copiers += [lambda x: pickle.loads(pickle.dumps(x, p)) for p in protocols]

    for obj in (Integer(-1), Integer(0), Integer(1), Rational(1, 2), pi,
                E, I, oo, -oo, zoo, nan, GoldenRatio, EulerGamma,
                Catalan, S.EmptySet, S.IdentityFunction):
        for func in copiers:
            assert func(obj) is obj


def test_functions():
    one_var = (acosh, ln, Heaviside, factorial, bernoulli, coth, tanh,
               sign, arg, asin, DiracDelta, re, Abs, sinh, cos, cot, acos,
               acot, gamma, bell, harmonic, LambertW, zeta, log, factorial,
               asinh, acoth, cosh, dirichlet_eta, loggamma, erf, ceiling,
               im, fibonacci, conjugate, tan, floor, atanh, sin, atan,
               lucas, exp)
    two_var = (rf, ff, lowergamma, chebyshevu, chebyshevt, binomial,
               atan2, polygamma, hermite, legendre, uppergamma)
    others = (chebyshevt_root, chebyshevu_root, Eijk(x, y, z),
              Piecewise( (0, x < -1), (x**2, x <= 1), (x**3, True)),
              assoc_legendre)
    for cls in one_var:
        check(cls)
        c = cls(x)
        check(c)
    for cls in two_var:
        check(cls)
        c = cls(x, y)
        check(c)
    for cls in others:
        check(cls)


def test_geometry():
    p1 = Point(1, 2)
    p2 = Point(2, 3)
    p3 = Point(0, 0)
    p4 = Point(0, 1)
    for c in (GeometryEntity, GeometryEntity(), Point, p1, Circle,
              Circle(p1, 2), Ellipse, Ellipse(p1, 3, 4), Line, Line(p1, p2),
              LinearEntity, LinearEntity(p1, p2), Ray, Ray(p1, p2), Segment,
              Segment(p1, p2), Polygon, Polygon(p1, p2, p3, p4),
              RegularPolygon, RegularPolygon(p1, 4, 5), Triangle,
              Triangle(p1, p2, p3)):
        check(c, check_attr=False)


def test_integrals():
    for c in (Integral, Integral(x)):
        check(c)


def test_logic():
    for c in (Logic, Logic(1)):
        check(c)


def test_matrices():
    for c in (Matrix, Matrix([1, 2, 3]), SparseMatrix,
              SparseMatrix([[1, 2], [3, 4]])):
        check(c)


def test_ntheory():
    for c in (Sieve, Sieve()):
        check(c)


def test_plotting():
    for c in (Plot, Plot(1, visible=False)):
        check(c)


# ================= polys =======================

def test_pickling_polys_polytools():
    for c in (Poly, Poly(x, x), PurePoly, PurePoly(x)):
        check(c)

    for c in (GroebnerBasis, GroebnerBasis([x**2 - 1], x)):
        check(c)


def test_pickling_polys_rings():
    # NOTE: can't use protocols < 2 because we have to execute __new__ to
    # make sure caching of rings works properly.

    ring = PolynomialRing(ZZ, 'x,y,z')

    for c in (PolynomialRing, ring):
        check(c, exclude=[0, 1])

    for c in (ring.one, ring.x):
        check(c, exclude=[0, 1])


def test_pickling_polys_fields():
    # NOTE: can't use protocols < 2 because we have to execute __new__ to
    # make sure caching of fields works properly.

    field = FractionField(ZZ, 'x,y,z')

    for c in (FractionField, field):
        check(c, exclude=[0, 1])

    for c in (field.one, field.x):
        check(c, exclude=[0, 1])


def test_pickling_polys_elements():
    for c in (PythonRational, PythonRational(1, 7)):
        check(c)

    gf17 = PythonFiniteField(17)
    gf64 = PythonFiniteField(64)

    for c in (gf17(5), gf64(12)):
        check(c, exclude=[0, 1])

    A = AlgebraicField(QQ, sqrt(2))

    for c in (A.one, A.unit, A([2, 1])):
        check(c, exclude=[0, 1])

    R = RealField(100)

    for c in (R.zero, R.one, R(1.2345)):
        check(c, exclude=[0, 1])

    C = ComplexField(100)

    for c in (C.zero, C.one, C(1.2345)):
        check(c, exclude=[0, 1])


def test_pickling_polys_domains():
    for c in (PythonIntegerRing, PythonIntegerRing()):
        check(c)

    for c in (PythonRationalField, PythonRationalField()):
        check(c)

    for c in (PythonFiniteField, PythonFiniteField(7), PythonFiniteField(64)):
        check(c, exclude=[0, 1])

    if HAS_GMPY:
        for c in (GMPYIntegerRing, GMPYIntegerRing()):
            check(c)

        for c in (GMPYRationalField, GMPYRationalField()):
            check(c)

        for c in (GMPYFiniteField, GMPYFiniteField(7), GMPYFiniteField(64)):
            check(c, exclude=[0, 1])

    for c in (RealField, RealField(100)):
        check(c, exclude=[0, 1])

    for c in (ComplexField, ComplexField(100)):
        check(c, exclude=[0, 1])

    for c in (AlgebraicField, AlgebraicField(QQ, sqrt(2))):
        check(c, exclude=[0, 1])

    EX = ExpressionDomain()
    for c in (ExpressionDomain, EX, EX(sin(x))):
        check(c)


def first_two(m):
    return m[:2]


def after_second(m):
    return m[2:]


def test_pickling_polys_orderings():
    for c in (LexOrder, LexOrder()):
        check(c)

    for c in (GradedLexOrder, GradedLexOrder()):
        check(c)

    for c in (ReversedGradedLexOrder, ReversedGradedLexOrder()):
        check(c)

    for c in (ProductOrder, ProductOrder((LexOrder(), first_two),
                                         (GradedLexOrder(), after_second))):
        check(c)

    for c in (InverseOrder, InverseOrder(LexOrder())):
        check(c)


def test_pickling_polys_monomials():
    for c in (Monomial, Monomial((1, 2, 3), (x, y, z))):
        check(c)


def test_pickling_polys_errors():
    # TODO: TypeError: __init__() takes at least 3 arguments (1 given)
    # for c in (ExactQuotientFailed, ExactQuotientFailed(x, 3*x, ZZ)):
    #    check(c)

    # TODO: TypeError: can't pickle instancemethod objects
    # for c in (OperationNotSupported, OperationNotSupported(Poly(x), Poly.gcd)):
    #    check(c)

    for c in (HeuristicGCDFailed, HeuristicGCDFailed(),
              HomomorphismFailed, HomomorphismFailed(),
              IsomorphismFailed, IsomorphismFailed(),
              ExtraneousFactors, ExtraneousFactors(),
              EvaluationFailed, EvaluationFailed(),
              RefinementFailed, RefinementFailed(),
              CoercionFailed, CoercionFailed(),
              NotInvertible, NotInvertible(),
              NotReversible, NotReversible(),
              NotAlgebraic, NotAlgebraic(),
              DomainError, DomainError(),
              PolynomialError, PolynomialError(),
              UnificationFailed, UnificationFailed(),
              GeneratorsError, GeneratorsError(),
              GeneratorsNeeded, GeneratorsNeeded()):
        check(c)

    # TODO: PicklingError: Can't pickle <function <lambda> at 0x38578c0>: it's not found as __main__.<lambda>
    # for c in (ComputationFailed, ComputationFailed(lambda t: t, 3, None)):
    #    check(c)

    for c in (UnivariatePolynomialError, UnivariatePolynomialError(),
              MultivariatePolynomialError, MultivariatePolynomialError()):
        check(c)

    # TODO: TypeError: __init__() takes at least 3 arguments (1 given)
    # for c in (PolificationFailed, PolificationFailed({}, x, x, False)):
    #    check(c)

    for c in (OptionError, OptionError(), FlagError, FlagError()):
        check(c)


@pytest.mark.xfail
def test_pickling_polys_options():
    for c in (Options, Options((), {'domain': 'ZZ', 'polys': False})):
        check(c)


def test_pickling_polys_rootoftools():
    f = x**3 + x + 3

    for c in (RootOf, RootOf(f, 0)):
        check(c)

    for c in (RootSum, RootSum(f, Lambda(x, exp(x)))):
        check(c)


def test_printing():
    for c in (LatexPrinter, LatexPrinter(), MathMLPrinter,
              PrettyPrinter, prettyForm, stringPict, stringPict('a'),
              Printer, Printer(), PythonPrinter, PythonPrinter()):
        check(c)


def test_series():
    for c in (Limit, Limit(y, x, 1), Order, Order(y)):
        check(c)


def test_concrete():
    for c in (Product, Product(x, (x, 2, 4)), Sum, Sum(x, (x, 2, 4))):
        check(c)


def test_sympyissue_7457():
    pickle.loads(pickle.dumps(Point(1.1, 2.1).evalf()))  # not raises

    a = Float('1.2')
    b = pickle.loads(pickle.dumps(a))
    b.evalf(strict=False)  # not raises
    assert a == b


def test_sympyissue_21121():
    x1 = Symbol('x', real=True)
    x2 = Symbol('x', real=True)

    assert hash(x1) == hash(x2)

    x2_pickled = pickle.dumps(x2)
    x2_unpickled = pickle.loads(x2_pickled)

    assert hash(x1) == hash(x2_unpickled)
