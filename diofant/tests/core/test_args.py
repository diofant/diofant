"""Test whether all elements of cls.args are instances of Basic."""

# NOTE: keep tests sorted by (module, class name) key.

import inspect
import io
import os
import re
import warnings

from diofant import (ITE, Add, Adjoint, And, Atom, AtomicExpr, Basic,
                     BlockDiagMatrix, BlockMatrix, Complement, Contains,
                     Derivative, Determinant, DiagonalMatrix, DiagonalOf, Dict,
                     Dummy, Equality, Equivalent, Expr, FiniteSet, Float,
                     FunctionMatrix, GrayCode, GreaterThan, GroebnerBasis,
                     HadamardProduct, Identity, Idx, ImageSet, ImmutableMatrix,
                     ImmutableSparseMatrix, Implies, Indexed, IndexedBase,
                     Integer, IntegerPartition, Integral, Intersection,
                     Interval, Inverse, Lambda, LessThan, Limit, MatAdd,
                     MatMul, MatPow, MatrixSlice, MatrixSymbol, Mod, Mul, Nand,
                     Nor, Not, Number, Or, Order, Partition, Permutation,
                     PermutationGroup, Poly, Polyhedron, Pow, Product,
                     ProductSet, Prufer, PurePoly, Range, Rational, RootOf,
                     RootSum, S, Set, StrictGreaterThan, StrictLessThan, Subs,
                     Subset, Sum, Symbol, SymmetricDifference, Trace,
                     Transpose, Tuple, Unequality, Union, Wild, WildFunction,
                     Xor, ZeroMatrix, divisor_sigma, false, mobius, oo, sin,
                     symbols, totient, true)
from diofant.abc import a, b, c, w, x, y, z
from diofant.concrete.expr_with_intlimits import ExprWithIntLimits
from diofant.concrete.expr_with_limits import AddWithLimits, ExprWithLimits
from diofant.core.function import Application, AppliedUndef
from diofant.core.numbers import (Catalan, ComplexInfinity, EulerGamma, Exp1,
                                  GoldenRatio, Half, ImaginaryUnit, Infinity,
                                  NaN, NegativeInfinity, NegativeOne,
                                  NumberSymbol, One, Pi, Zero)
from diofant.core.trace import Tr
from diofant.diffgeom import (BaseCovarDerivativeOp, BaseScalarField,
                              BaseVectorField, Commutator, CoordSystem,
                              CovarDerivativeOp, Differential, LieDerivative,
                              Manifold, Patch)
from diofant.diffgeom import Point as DiffgeomPoint
from diofant.diffgeom import TensorProduct, WedgeProduct
from diofant.functions import (Chi, Ci, DiracDelta, Ei, FallingFactorial,
                               Heaviside, KroneckerDelta, LambertW, LeviCivita,
                               Li, Max, Min, Piecewise, RisingFactorial, Shi,
                               Si, Ynm, Znm, acos, acosh, acot, acoth, acsc,
                               adjoint, airyai, airyaiprime, airybi,
                               airybiprime, arg, asec, asin, asinh,
                               assoc_laguerre, assoc_legendre, atan, atan2,
                               atanh, bell, bernoulli, besseli, besselj,
                               besselk, bessely, beta, binomial, catalan,
                               ceiling, chebyshevt, chebyshevt_root,
                               chebyshevu, chebyshevu_root, conjugate, cos,
                               cosh, cot, coth, csc, csch, dirichlet_eta,
                               elliptic_e, elliptic_f, elliptic_k, elliptic_pi,
                               erf, erf2, erf2inv, erfc, erfcinv, erfi, erfinv,
                               euler, exp, exp_polar, expint, factorial,
                               factorial2, fibonacci, floor, fresnelc,
                               fresnels, gamma, gegenbauer, genocchi, hankel1,
                               hankel2, harmonic, hermite, hyper, im, jacobi,
                               jn, laguerre, legendre, lerchphi, li, log,
                               loggamma, lowergamma, lucas, meijerg,
                               periodic_argument, polar_lift, polygamma,
                               polylog, principal_branch)
from diofant.functions import re as _re
from diofant.functions import (sec, sech, sign, sinh, subfactorial, tan, tanh,
                               transpose, uppergamma, yn, zeta)
from diofant.functions.elementary.miscellaneous import IdentityFunction
from diofant.functions.elementary.piecewise import ExprCondPair
from diofant.functions.special.error_functions import _erfs
from diofant.functions.special.hyper import (HyperRep_asin1, HyperRep_asin2,
                                             HyperRep_atanh, HyperRep_cosasin,
                                             HyperRep_log1, HyperRep_log2,
                                             HyperRep_power1, HyperRep_power2,
                                             HyperRep_sinasin, HyperRep_sqrts1,
                                             HyperRep_sqrts2)
from diofant.geometry import (Circle, Curve, Ellipse, Line, Point, Polygon,
                              Ray, RegularPolygon, Segment, Triangle)
from diofant.geometry.entity import GeometryEntity
from diofant.integrals.risch import NonElementaryIntegral
from diofant.integrals.transforms import (CosineTransform, FourierTransform,
                                          HankelTransform,
                                          InverseCosineTransform,
                                          InverseFourierTransform,
                                          InverseHankelTransform,
                                          InverseLaplaceTransform,
                                          InverseMellinTransform,
                                          InverseSineTransform,
                                          LaplaceTransform, MellinTransform,
                                          SineTransform)
from diofant.logic.boolalg import BooleanFunction
from diofant.matrices.expressions.fourier import DFT, IDFT
from diofant.matrices.expressions.matexpr import MatrixElement
from diofant.printing.codeprinter import Assignment
from diofant.sets.fancysets import (ExtendedReals, Integers, Naturals,
                                    Naturals0, Rationals, Reals)
from diofant.sets.sets import EmptySet, UniversalSet
from diofant.simplify.hyperexpand import G_Function, Hyper_Function
from diofant.stats.crv import (ConditionalContinuousDomain,
                               ContinuousDistributionHandmade,
                               ContinuousDomain, ContinuousPSpace,
                               ProductContinuousDomain,
                               ProductContinuousPSpace, SingleContinuousDomain,
                               SingleContinuousPSpace)
from diofant.stats.crv_types import (ArcsinDistribution, BeniniDistribution,
                                     BetaDistribution, BetaPrimeDistribution,
                                     CauchyDistribution, ChiDistribution,
                                     ChiNoncentralDistribution,
                                     ChiSquaredDistribution, DagumDistribution,
                                     ExponentialDistribution,
                                     FDistributionDistribution,
                                     FisherZDistribution, FrechetDistribution,
                                     GammaDistribution,
                                     GammaInverseDistribution,
                                     KumaraswamyDistribution,
                                     LaplaceDistribution, LogisticDistribution,
                                     LogNormalDistribution,
                                     MaxwellDistribution, NakagamiDistribution,
                                     Normal, NormalDistribution,
                                     ParetoDistribution,
                                     QuadraticUDistribution,
                                     RaisedCosineDistribution,
                                     RayleighDistribution,
                                     StudentTDistribution,
                                     TriangularDistribution,
                                     UniformDistribution,
                                     UniformSumDistribution,
                                     VonMisesDistribution, WeibullDistribution,
                                     WignerSemicircleDistribution)
from diofant.stats.drv import SingleDiscreteDomain, SingleDiscretePSpace
from diofant.stats.drv_types import GeometricDistribution, PoissonDistribution
from diofant.stats.frv import (ConditionalFiniteDomain, FiniteDomain,
                               FinitePSpace, ProductFiniteDomain,
                               ProductFinitePSpace, SingleFiniteDomain,
                               SingleFinitePSpace)
from diofant.stats.frv_types import (BernoulliDistribution,
                                     BinomialDistribution, DieDistribution,
                                     DiscreteUniformDistribution,
                                     FiniteDistributionHandmade,
                                     HypergeometricDistribution,
                                     RademacherDistribution)
from diofant.stats.rv import (ConditionalDomain, Density, ProductDomain,
                              ProductPSpace, PSpace, RandomDomain,
                              RandomSymbol, SingleDomain)
from diofant.tensor import ImmutableDenseNDimArray, ImmutableSparseNDimArray
from diofant.tensor.tensor import (TensAdd, TensorHead, TensorIndex,
                                   TensorIndexType, TensorSymmetry, TensorType,
                                   get_symmetric_group_sgs, tensor_indices)
from diofant.utilities.exceptions import DiofantDeprecationWarning
from diofant.vector import (AxisOrienter, BaseDyadic, BaseScalar, BaseVector,
                            BodyOrienter, CoordSysCartesian, Del, DyadicAdd,
                            DyadicMul, DyadicZero)
from diofant.vector import Point as VPoint
from diofant.vector import (QuaternionOrienter, SpaceOrienter, VectorAdd,
                            VectorMul, VectorZero)


__all__ = ()


def test_all_classes_are_tested():
    this = os.path.split(__file__)[0]
    path = os.path.join(this, os.pardir, os.pardir)
    diofant_path = os.path.abspath(path)
    prefix = os.path.split(diofant_path)[0] + os.sep

    re_cls = re.compile(r'^class ([A-Za-z][A-Za-z0-9_]*)\s*\(', re.MULTILINE)

    modules = {}

    for root, dirs, files in os.walk(diofant_path):
        module = root.replace(prefix, '').replace(os.sep, '.')

        for file in files:
            if file.startswith(('_', 'test_', 'bench_')):
                continue
            if not file.endswith('.py'):
                continue

            with io.open(os.path.join(root, file), 'r', encoding='utf-8') as f:
                text = f.read()

            submodule = module + '.' + file[:-3]
            names = re_cls.findall(text)

            if not names:
                continue

            try:
                mod = __import__(submodule, fromlist=names)
            except ImportError:
                continue

            def is_Basic(name):
                cls = getattr(mod, name)
                return issubclass(cls, Basic)

            names = list(filter(is_Basic, names))

            if names:
                modules[submodule] = names

    ns = globals()
    failed = []

    for module, names in modules.items():
        mod = module.replace('.', '__')

        for name in names:
            test = 'test_' + mod + '__' + name

            if test not in ns:
                failed.append(module + '.' + name)

    # reset all DiofantDeprecationWarning into errors
    warnings.simplefilter('error', category=DiofantDeprecationWarning)

    assert not failed, f"Missing classes: {', '.join(failed)}.  Please add tests for these to diofant/core/tests/test_args.py."


def _test_args(obj):
    all_basic = all(isinstance(arg, Basic) for arg in obj.args)

    # Ideally obj.func(*obj.args) would always recreate the object, but for
    # now, we only require it for objects with non-empty .args
    recreatable = not obj.args or obj.func(*obj.args) == obj

    res = all_basic and recreatable

    if hasattr(obj, 'doit'):
        doit = obj.doit
        if inspect.ismethod(doit):
            spec = inspect.getfullargspec(doit)
            res &= (len(spec.args) == 1 and spec.varargs is None
                    and spec.varkw is not None)
        else:
            res &= False

    return res


def test_diofant__combinatorics__graycode__GrayCode():
    # an integer is given and returned from GrayCode as the arg
    assert _test_args(GrayCode(3, start='100'))
    assert _test_args(GrayCode(3, rank=1))


def test_diofant__combinatorics__subsets__Subset():
    assert _test_args(Subset([0, 1], [0, 1, 2, 3]))
    assert _test_args(Subset(['c', 'd'], ['a', 'b', 'c', 'd']))


def test_diofant__combinatorics__permutations__Permutation():
    assert _test_args(Permutation([0, 1, 2, 3]))


def test_diofant__combinatorics__perm_groups__PermutationGroup():
    assert _test_args(PermutationGroup([Permutation([0, 1])]))


def test_diofant__combinatorics__polyhedron__Polyhedron():
    pgroup = [Permutation([[0, 1, 2], [3]]),
              Permutation([[0, 1, 3], [2]]),
              Permutation([[0, 2, 3], [1]]),
              Permutation([[1, 2, 3], [0]]),
              Permutation([[0, 1], [2, 3]]),
              Permutation([[0, 2], [1, 3]]),
              Permutation([[0, 3], [1, 2]]),
              Permutation([[0, 1, 2, 3]])]
    corners = [w, x, y, z]
    faces = [(w, x, y), (w, y, z), (w, z, x), (x, y, z)]
    assert _test_args(Polyhedron(corners, faces, pgroup))


def test_diofant__combinatorics__prufer__Prufer():
    assert _test_args(Prufer([[0, 1], [0, 2], [0, 3]], 4))


def test_diofant__combinatorics__partitions__Partition():
    assert _test_args(Partition([1]))


def test_diofant__combinatorics__partitions__IntegerPartition():
    assert _test_args(IntegerPartition([1]))


def test_diofant__concrete__products__Product():
    assert _test_args(Product(x, (x, 0, 10)))
    assert _test_args(Product(x, (x, 0, y), (y, 0, 10)))


def test_diofant__concrete__expr_with_limits__ExprWithLimits():
    assert _test_args(ExprWithLimits(x, (x, 0, 10)))
    assert _test_args(ExprWithLimits(x*y, (x, 0, 10.), (y, 1., 3)))


def test_diofant__concrete__expr_with_limits__AddWithLimits():
    assert _test_args(AddWithLimits(x, (x, 0, 10)))
    assert _test_args(AddWithLimits(x*y, (x, 0, 10), (y, 1, 3)))


def test_diofant__concrete__expr_with_intlimits__ExprWithIntLimits():
    assert _test_args(ExprWithIntLimits(x, (x, 0, 10)))
    assert _test_args(ExprWithIntLimits(x*y, (x, 0, 10), (y, 1, 3)))


def test_diofant__concrete__summations__Sum():
    assert _test_args(Sum(x, (x, 0, 10)))
    assert _test_args(Sum(x, (x, 0, y), (y, 0, 10)))


def test_diofant__core__add__Add():
    assert _test_args(Add(x, y, z, 2))


def test_diofant__core__basic__Atom():
    assert _test_args(Atom())


def test_diofant__core__basic__Basic():
    assert _test_args(Basic())


def test_diofant__core__containers__Dict():
    assert _test_args(Dict({x: y, y: z}))


def test_diofant__core__containers__Tuple():
    assert _test_args(Tuple(x, y, z, 2))


def test_diofant__core__expr__AtomicExpr():
    assert _test_args(AtomicExpr())


def test_diofant__core__expr__Expr():
    assert _test_args(Expr())


def test_diofant__core__function__Application():
    assert _test_args(Application(1, 2, 3))


def test_diofant__core__function__AppliedUndef():
    assert _test_args(AppliedUndef(1, 2, 3))


def test_diofant__core__function__Derivative():
    assert _test_args(Derivative(2, x, (y, 3)))


def test_diofant__core__function__Function():
    pass


def test_diofant__core__function__Lambda():
    assert _test_args(Lambda((x, y), x + y + z))


def test_diofant__core__function__Subs():
    assert _test_args(Subs(x + y, (x, 2)))
    assert _test_args(Subs(x + y, (x, 2), (y, 1)))


def test_diofant__core__function__WildFunction():
    assert _test_args(WildFunction('f'))


def test_diofant__core__mod__Mod():
    assert _test_args(Mod(x, 2))


def test_diofant__core__mul__Mul():
    assert _test_args(Mul(2, x, y, z))


def test_diofant__core__numbers__Catalan():
    assert _test_args(Catalan())


def test_diofant__core__numbers__ComplexInfinity():
    assert _test_args(ComplexInfinity())


def test_diofant__core__numbers__EulerGamma():
    assert _test_args(EulerGamma())


def test_diofant__core__numbers__Exp1():
    assert _test_args(Exp1())


def test_diofant__core__numbers__Float():
    assert _test_args(Float(1.23))


def test_diofant__core__numbers__GoldenRatio():
    assert _test_args(GoldenRatio())


def test_diofant__core__numbers__Half():
    assert _test_args(Half())


def test_diofant__core__numbers__ImaginaryUnit():
    assert _test_args(ImaginaryUnit())


def test_diofant__core__numbers__Infinity():
    assert _test_args(Infinity())


def test_diofant__core__numbers__Integer():
    assert _test_args(Integer(7))


def test_diofant__core__numbers__IntegerConstant():
    pass


def test_diofant__core__numbers__NaN():
    assert _test_args(NaN())


def test_diofant__core__numbers__NegativeInfinity():
    assert _test_args(NegativeInfinity())


def test_diofant__core__numbers__NegativeOne():
    assert _test_args(NegativeOne())


def test_diofant__core__numbers__Number():
    assert _test_args(Number(1, 7))


def test_diofant__core__numbers__NumberSymbol():
    assert _test_args(NumberSymbol())


def test_diofant__core__numbers__One():
    assert _test_args(One())


def test_diofant__core__numbers__Pi():
    assert _test_args(Pi())


def test_diofant__core__numbers__Rational():
    assert _test_args(Rational(1, 7))


def test_diofant__core__numbers__RationalConstant():
    pass


def test_diofant__core__numbers__Zero():
    assert _test_args(Zero())


def test_diofant__core__operations__AssocOp():
    pass


def test_diofant__core__operations__LatticeOp():
    pass


def test_diofant__core__power__Pow():
    assert _test_args(Pow(x, 2))


def test_diofant__core__relational__Equality():
    assert _test_args(Equality(x, 2))


def test_diofant__core__relational__GreaterThan():
    assert _test_args(GreaterThan(x, 2))


def test_diofant__core__relational__LessThan():
    assert _test_args(LessThan(x, 2))


def test_diofant__core__relational__Relational():
    pass


def test_diofant__core__relational__StrictGreaterThan():
    assert _test_args(StrictGreaterThan(x, 2))


def test_diofant__core__relational__StrictLessThan():
    assert _test_args(StrictLessThan(x, 2))


def test_diofant__core__relational__Unequality():
    assert _test_args(Unequality(x, 2))


def test_diofant__sets__sets__EmptySet():
    assert _test_args(EmptySet())


def test_diofant__sets__sets__UniversalSet():
    assert _test_args(UniversalSet())


def test_diofant__sets__sets__FiniteSet():
    assert _test_args(FiniteSet(x, y, z))


def test_diofant__sets__sets__Interval():
    assert _test_args(Interval(0, 1))


def test_diofant__sets__sets__ProductSet():
    assert _test_args(ProductSet(Interval(0, 1), Interval(0, 1)))


def test_diofant__sets__sets__Set():
    assert _test_args(Set())


def test_diofant__sets__sets__Intersection():
    assert _test_args(Intersection(Interval(0, 3), Interval(2, x)))


def test_diofant__sets__sets__Union():
    assert _test_args(Union(Interval(0, 1), Interval(2, 3)))


def test_diofant__sets__sets__Complement():
    assert _test_args(Complement(Interval(0, 2), Interval(0, 1)))


def test_diofant__sets__sets__SymmetricDifference():
    assert _test_args(SymmetricDifference(FiniteSet(1, 2, 3),
                                          FiniteSet(2, 3, 4)))


def test_diofant__core__trace__Tr():
    a, b = symbols('a b')
    assert _test_args(Tr(a + b))


def test_diofant__sets__fancysets__Naturals():
    assert _test_args(Naturals())


def test_diofant__sets__fancysets__Naturals0():
    assert _test_args(Naturals0())


def test_diofant__sets__fancysets__Integers():
    assert _test_args(Integers())


def test_diofant__sets__fancysets__Rationals():
    assert _test_args(Rationals())


def test_diofant__sets__fancysets__Reals():
    assert _test_args(Reals())


def test_diofant__sets__fancysets__ExtendedReals():
    assert _test_args(ExtendedReals())


def test_diofant__sets__fancysets__ImageSet():
    x = Symbol('x')
    assert _test_args(ImageSet(Lambda(x, x**2), S.Naturals))


def test_diofant__sets__fancysets__Range():
    assert _test_args(Range(1, 5, 1))


def test_diofant__sets__contains__Contains():
    assert _test_args(Contains(x, Range(0, 10, 2)))


# STATS


nd = NormalDistribution(0, 1)
die = DieDistribution(6)


def test_diofant__stats__crv__ContinuousDomain():
    assert _test_args(ContinuousDomain({x}, Interval(-oo, oo)))


def test_diofant__stats__crv__SingleContinuousDomain():
    assert _test_args(SingleContinuousDomain(x, Interval(-oo, oo)))


def test_diofant__stats__crv__ProductContinuousDomain():
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    E = SingleContinuousDomain(y, Interval(0, oo))
    assert _test_args(ProductContinuousDomain(D, E))


def test_diofant__stats__crv__ConditionalContinuousDomain():
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    assert _test_args(ConditionalContinuousDomain(D, x > 0))


def test_diofant__stats__crv__ContinuousPSpace():
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    assert _test_args(ContinuousPSpace(D, nd))


def test_diofant__stats__crv__SingleContinuousPSpace():
    assert _test_args(SingleContinuousPSpace(x, nd))


def test_diofant__stats__crv__ProductContinuousPSpace():
    A = SingleContinuousPSpace(x, nd)
    B = SingleContinuousPSpace(y, nd)
    assert _test_args(ProductContinuousPSpace(A, B))


def test_diofant__stats__crv__SingleContinuousDistribution():
    pass


def test_diofant__stats__drv__SingleDiscreteDomain():
    assert _test_args(SingleDiscreteDomain(x, S.Naturals))


def test_diofant__stats__drv__SingleDiscretePSpace():
    assert _test_args(SingleDiscretePSpace(x, PoissonDistribution(1)))


def test_diofant__stats__drv__SingleDiscreteDistribution():
    pass


def test_diofant__stats__rv__RandomDomain():
    assert _test_args(RandomDomain(FiniteSet(x), FiniteSet(1, 2, 3)))


def test_diofant__stats__rv__SingleDomain():
    assert _test_args(SingleDomain(x, FiniteSet(1, 2, 3)))


def test_diofant__stats__rv__ConditionalDomain():
    D = RandomDomain(FiniteSet(x), FiniteSet(1, 2))
    assert _test_args(ConditionalDomain(D, x > 1))


def test_diofant__stats__rv__PSpace():
    D = RandomDomain(FiniteSet(x), FiniteSet(1, 2, 3, 4, 5, 6))
    assert _test_args(PSpace(D, die))


def test_diofant__stats__rv__SinglePSpace():
    pass


def test_diofant__stats__rv__RandomSymbol():
    A = SingleContinuousPSpace(x, nd)
    assert _test_args(RandomSymbol(A, x))


def test_diofant__stats__rv__ProductPSpace():
    A = SingleContinuousPSpace(x, nd)
    B = SingleContinuousPSpace(y, nd)
    assert _test_args(ProductPSpace(A, B))


def test_diofant__stats__rv__ProductDomain():
    D = SingleDomain(x, Interval(-oo, oo))
    E = SingleDomain(y, Interval(0, oo))
    assert _test_args(ProductDomain(D, E))


def test_diofant__stats__frv_types__DiscreteUniformDistribution():
    assert _test_args(DiscreteUniformDistribution(Tuple(*range(6))))


def test_diofant__stats__frv_types__DieDistribution():
    assert _test_args(DieDistribution(6))


def test_diofant__stats__frv_types__BernoulliDistribution():
    assert _test_args(BernoulliDistribution(Rational(1, 2), 0, 1))


def test_diofant__stats__frv_types__BinomialDistribution():
    assert _test_args(BinomialDistribution(5, Rational(1, 2), 1, 0))


def test_diofant__stats__frv_types__HypergeometricDistribution():
    assert _test_args(HypergeometricDistribution(10, 5, 3))


def test_diofant__stats__frv_types__RademacherDistribution():
    assert _test_args(RademacherDistribution())


def test_diofant__stats__frv__FiniteDomain():
    assert _test_args(FiniteDomain({(x, 1), (x, 2)}))  # x can be 1 or 2


def test_diofant__stats__frv__SingleFiniteDomain():
    assert _test_args(SingleFiniteDomain(x, {1, 2}))  # x can be 1 or 2


def test_diofant__stats__frv__ProductFiniteDomain():
    xd = SingleFiniteDomain(x, {1, 2})
    yd = SingleFiniteDomain(y, {1, 2})
    assert _test_args(ProductFiniteDomain(xd, yd))


def test_diofant__stats__frv__ConditionalFiniteDomain():
    xd = SingleFiniteDomain(x, {1, 2})
    assert _test_args(ConditionalFiniteDomain(xd, x > 1))


def test_diofant__stats__frv__FinitePSpace():
    xd = SingleFiniteDomain(x, {1, 2, 3, 4, 5, 6})
    xd = SingleFiniteDomain(x, {1, 2})
    assert _test_args(FinitePSpace(xd, {(x, 1): Rational(1, 2), (x, 2): Rational(1, 2)}))


def test_diofant__stats__frv__SingleFinitePSpace():
    assert _test_args(SingleFinitePSpace(Symbol('x'), die))


def test_diofant__stats__frv__ProductFinitePSpace():
    xp = SingleFinitePSpace(Symbol('x'), die)
    yp = SingleFinitePSpace(Symbol('y'), die)
    assert _test_args(ProductFinitePSpace(xp, yp))


def test_diofant__stats__frv__SingleFiniteDistribution():
    pass


def test_diofant__stats__crv__ContinuousDistribution():
    pass


def test_diofant__stats__frv_types__FiniteDistributionHandmade():
    assert _test_args(FiniteDistributionHandmade({1: 1}))


def test_diofant__stats__crv__ContinuousDistributionHandmade():
    assert _test_args(ContinuousDistributionHandmade(Symbol('x'),
                                                     Interval(0, 2)))


def test_diofant__stats__rv__Density():
    assert _test_args(Density(Normal('x', 0, 1)))


def test_diofant__stats__crv_types__ArcsinDistribution():
    assert _test_args(ArcsinDistribution(0, 1))


def test_diofant__stats__crv_types__BeniniDistribution():
    assert _test_args(BeniniDistribution(1, 1, 1))


def test_diofant__stats__crv_types__BetaDistribution():
    assert _test_args(BetaDistribution(1, 1))


def test_diofant__stats__crv_types__BetaPrimeDistribution():
    assert _test_args(BetaPrimeDistribution(1, 1))


def test_diofant__stats__crv_types__CauchyDistribution():
    assert _test_args(CauchyDistribution(0, 1))


def test_diofant__stats__crv_types__ChiDistribution():
    assert _test_args(ChiDistribution(1))


def test_diofant__stats__crv_types__ChiNoncentralDistribution():
    assert _test_args(ChiNoncentralDistribution(1, 1))


def test_diofant__stats__crv_types__ChiSquaredDistribution():
    assert _test_args(ChiSquaredDistribution(1))


def test_diofant__stats__crv_types__DagumDistribution():
    assert _test_args(DagumDistribution(1, 1, 1))


def test_diofant__stats__crv_types__ExponentialDistribution():
    assert _test_args(ExponentialDistribution(1))


def test_diofant__stats__crv_types__FDistributionDistribution():
    assert _test_args(FDistributionDistribution(1, 1))


def test_diofant__stats__crv_types__FisherZDistribution():
    assert _test_args(FisherZDistribution(1, 1))


def test_diofant__stats__crv_types__FrechetDistribution():
    assert _test_args(FrechetDistribution(1, 1, 1))


def test_diofant__stats__crv_types__GammaInverseDistribution():
    assert _test_args(GammaInverseDistribution(1, 1))


def test_diofant__stats__crv_types__GammaDistribution():
    assert _test_args(GammaDistribution(1, 1))


def test_diofant__stats__crv_types__KumaraswamyDistribution():
    assert _test_args(KumaraswamyDistribution(1, 1))


def test_diofant__stats__crv_types__LaplaceDistribution():
    assert _test_args(LaplaceDistribution(0, 1))


def test_diofant__stats__crv_types__LogisticDistribution():
    assert _test_args(LogisticDistribution(0, 1))


def test_diofant__stats__crv_types__LogNormalDistribution():
    assert _test_args(LogNormalDistribution(0, 1))


def test_diofant__stats__crv_types__MaxwellDistribution():
    assert _test_args(MaxwellDistribution(1))


def test_diofant__stats__crv_types__NakagamiDistribution():
    assert _test_args(NakagamiDistribution(1, 1))


def test_diofant__stats__crv_types__NormalDistribution():
    assert _test_args(NormalDistribution(0, 1))


def test_diofant__stats__crv_types__ParetoDistribution():
    assert _test_args(ParetoDistribution(1, 1))


def test_diofant__stats__crv_types__QuadraticUDistribution():
    assert _test_args(QuadraticUDistribution(1, 2))


def test_diofant__stats__crv_types__RaisedCosineDistribution():
    assert _test_args(RaisedCosineDistribution(1, 1))


def test_diofant__stats__crv_types__RayleighDistribution():
    assert _test_args(RayleighDistribution(1))


def test_diofant__stats__crv_types__StudentTDistribution():
    assert _test_args(StudentTDistribution(1))


def test_diofant__stats__crv_types__TriangularDistribution():
    assert _test_args(TriangularDistribution(-1, 0, 1))


def test_diofant__stats__crv_types__UniformDistribution():
    assert _test_args(UniformDistribution(0, 1))


def test_diofant__stats__crv_types__UniformSumDistribution():
    assert _test_args(UniformSumDistribution(1))


def test_diofant__stats__crv_types__VonMisesDistribution():
    assert _test_args(VonMisesDistribution(1, 1))


def test_diofant__stats__crv_types__WeibullDistribution():
    assert _test_args(WeibullDistribution(1, 1))


def test_diofant__stats__crv_types__WignerSemicircleDistribution():
    assert _test_args(WignerSemicircleDistribution(1))


def test_diofant__stats__drv_types__PoissonDistribution():
    assert _test_args(PoissonDistribution(1))


def test_diofant__stats__drv_types__GeometricDistribution():
    assert _test_args(GeometricDistribution(.5))


def test_diofant__core__symbol__BaseSymbol():
    pass


def test_diofant__core__symbol__Dummy():
    assert _test_args(Dummy('t'))


def test_diofant__core__symbol__Symbol():
    assert _test_args(Symbol('t'))


def test_diofant__core__symbol__Wild():
    assert _test_args(Wild('x', exclude=[x]))


def test_diofant__functions__combinatorial__factorials__CombinatorialFunction():
    pass


def test_diofant__functions__combinatorial__factorials__FallingFactorial():
    assert _test_args(FallingFactorial(2, x))


def test_diofant__functions__combinatorial__factorials__RisingFactorial():
    assert _test_args(RisingFactorial(2, x))


def test_diofant__functions__combinatorial__factorials__binomial():
    assert _test_args(binomial(2, x))


def test_diofant__functions__combinatorial__factorials__subfactorial():
    assert _test_args(subfactorial(1))


def test_diofant__functions__combinatorial__factorials__factorial():
    assert _test_args(factorial(x))


def test_diofant__functions__combinatorial__factorials__factorial2():
    assert _test_args(factorial2(x))


def test_diofant__functions__combinatorial__numbers__bell():
    assert _test_args(bell(x, y))


def test_diofant__functions__combinatorial__numbers__bernoulli():
    assert _test_args(bernoulli(x))


def test_diofant__functions__combinatorial__numbers__catalan():
    assert _test_args(catalan(x))


def test_diofant__functions__combinatorial__numbers__genocchi():
    assert _test_args(genocchi(x))


def test_diofant__functions__combinatorial__numbers__euler():
    assert _test_args(euler(x))


def test_diofant__functions__combinatorial__numbers__fibonacci():
    assert _test_args(fibonacci(x))


def test_diofant__functions__combinatorial__numbers__harmonic():
    assert _test_args(harmonic(x, 2))


def test_diofant__functions__combinatorial__numbers__lucas():
    assert _test_args(lucas(x))


def test_diofant__functions__elementary__complexes__Abs():
    assert _test_args(abs(x))


def test_diofant__functions__elementary__complexes__adjoint():
    assert _test_args(adjoint(x))


def test_diofant__functions__elementary__complexes__arg():
    assert _test_args(arg(x))


def test_diofant__functions__elementary__complexes__conjugate():
    assert _test_args(conjugate(x))


def test_diofant__functions__elementary__complexes__im():
    assert _test_args(im(x))


def test_diofant__functions__elementary__complexes__re():
    assert _test_args(_re(x))


def test_diofant__functions__elementary__complexes__sign():
    assert _test_args(sign(x))


def test_diofant__functions__elementary__complexes__polar_lift():
    assert _test_args(polar_lift(x))


def test_diofant__functions__elementary__complexes__periodic_argument():
    assert _test_args(periodic_argument(x, y))


def test_diofant__functions__elementary__complexes__principal_branch():
    assert _test_args(principal_branch(x, y))


def test_diofant__functions__elementary__complexes__transpose():
    assert _test_args(transpose(x))


def test_diofant__functions__elementary__exponential__LambertW():
    assert _test_args(LambertW(2))


def test_diofant__functions__elementary__exponential__exp():
    assert _test_args(exp(2))


def test_diofant__functions__elementary__exponential__exp_polar():
    assert _test_args(exp_polar(2))


def test_diofant__functions__elementary__exponential__log():
    assert _test_args(log(2))


def test_diofant__functions__elementary__hyperbolic__HyperbolicFunction():
    pass


def test_diofant__functions__elementary__hyperbolic__ReciprocalHyperbolicFunction():
    pass


def test_diofant__functions__elementary__hyperbolic__acosh():
    assert _test_args(acosh(2))


def test_diofant__functions__elementary__hyperbolic__acoth():
    assert _test_args(acoth(2))


def test_diofant__functions__elementary__hyperbolic__asinh():
    assert _test_args(asinh(2))


def test_diofant__functions__elementary__hyperbolic__atanh():
    assert _test_args(atanh(2))


def test_diofant__functions__elementary__hyperbolic__cosh():
    assert _test_args(cosh(2))


def test_diofant__functions__elementary__hyperbolic__coth():
    assert _test_args(coth(2))


def test_diofant__functions__elementary__hyperbolic__csch():
    assert _test_args(csch(2))


def test_diofant__functions__elementary__hyperbolic__sech():
    assert _test_args(sech(2))


def test_diofant__functions__elementary__hyperbolic__sinh():
    assert _test_args(sinh(2))


def test_diofant__functions__elementary__hyperbolic__tanh():
    assert _test_args(tanh(2))


def test_diofant__functions__elementary__integers__RoundFunction():
    pass


def test_diofant__functions__elementary__integers__ceiling():
    assert _test_args(ceiling(x))


def test_diofant__functions__elementary__integers__floor():
    assert _test_args(floor(x))


def test_diofant__functions__elementary__miscellaneous__IdentityFunction():
    assert _test_args(IdentityFunction())


def test_diofant__functions__elementary__miscellaneous__Max():
    assert _test_args(Max(x, 2))


def test_diofant__functions__elementary__miscellaneous__Min():
    assert _test_args(Min(x, 2))


def test_diofant__functions__elementary__miscellaneous__MinMaxBase():
    pass


def test_diofant__functions__elementary__piecewise__ExprCondPair():
    assert _test_args(ExprCondPair(1, True))


def test_diofant__functions__elementary__piecewise__Piecewise():
    assert _test_args(Piecewise((1, x >= 0), (0, True)))


def test_diofant__functions__elementary__trigonometric__TrigonometricFunction():
    pass


def test_diofant__functions__elementary__trigonometric__ReciprocalTrigonometricFunction():
    pass


def test_diofant__functions__elementary__trigonometric__InverseTrigonometricFunction():
    pass


def test_diofant__functions__elementary__trigonometric__acos():
    assert _test_args(acos(2))


def test_diofant__functions__elementary__trigonometric__acot():
    assert _test_args(acot(2))


def test_diofant__functions__elementary__trigonometric__asin():
    assert _test_args(asin(2))


def test_diofant__functions__elementary__trigonometric__asec():
    assert _test_args(asec(2))


def test_diofant__functions__elementary__trigonometric__acsc():
    assert _test_args(acsc(2))


def test_diofant__functions__elementary__trigonometric__atan():
    assert _test_args(atan(2))


def test_diofant__functions__elementary__trigonometric__atan2():
    assert _test_args(atan2(2, 3))


def test_diofant__functions__elementary__trigonometric__cos():
    assert _test_args(cos(2))


def test_diofant__functions__elementary__trigonometric__csc():
    assert _test_args(csc(2))


def test_diofant__functions__elementary__trigonometric__cot():
    assert _test_args(cot(2))


def test_diofant__functions__elementary__trigonometric__sin():
    assert _test_args(sin(2))


def test_diofant__functions__elementary__trigonometric__sec():
    assert _test_args(sec(2))


def test_diofant__functions__elementary__trigonometric__tan():
    assert _test_args(tan(2))


def test_diofant__functions__special__bessel__BesselBase():
    pass


def test_diofant__functions__special__bessel__SphericalBesselBase():
    pass


def test_diofant__functions__special__bessel__besseli():
    assert _test_args(besseli(x, 1))


def test_diofant__functions__special__bessel__besselj():
    assert _test_args(besselj(x, 1))


def test_diofant__functions__special__bessel__besselk():
    assert _test_args(besselk(x, 1))


def test_diofant__functions__special__bessel__bessely():
    assert _test_args(bessely(x, 1))


def test_diofant__functions__special__bessel__hankel1():
    assert _test_args(hankel1(x, 1))


def test_diofant__functions__special__bessel__hankel2():
    assert _test_args(hankel2(x, 1))


def test_diofant__functions__special__bessel__jn():
    assert _test_args(jn(0, x))


def test_diofant__functions__special__bessel__yn():
    assert _test_args(yn(0, x))


def test_diofant__functions__special__bessel__AiryBase():
    pass


def test_diofant__functions__special__bessel__airyai():
    assert _test_args(airyai(2))


def test_diofant__functions__special__bessel__airybi():
    assert _test_args(airybi(2))


def test_diofant__functions__special__bessel__airyaiprime():
    assert _test_args(airyaiprime(2))


def test_diofant__functions__special__bessel__airybiprime():
    assert _test_args(airybiprime(2))


def test_diofant__functions__special__elliptic_integrals__elliptic_k():
    assert _test_args(elliptic_k(x))


def test_diofant__functions__special__elliptic_integrals__elliptic_f():
    assert _test_args(elliptic_f(x, y))


def test_diofant__functions__special__elliptic_integrals__elliptic_e():
    assert _test_args(elliptic_e(x))
    assert _test_args(elliptic_e(x, y))


def test_diofant__functions__special__elliptic_integrals__elliptic_pi():
    assert _test_args(elliptic_pi(x, y))
    assert _test_args(elliptic_pi(x, y, z))


def test_diofant__functions__special__delta_functions__DiracDelta():
    assert _test_args(DiracDelta(x, 1))


def test_diofant__functions__special__delta_functions__Heaviside():
    assert _test_args(Heaviside(x))


def test_diofant__functions__special__error_functions__erf():
    assert _test_args(erf(2))


def test_diofant__functions__special__error_functions__erfc():
    assert _test_args(erfc(2))


def test_diofant__functions__special__error_functions__erfi():
    assert _test_args(erfi(2))


def test_diofant__functions__special__error_functions__erf2():
    assert _test_args(erf2(2, 3))


def test_diofant__functions__special__error_functions__erfinv():
    assert _test_args(erfinv(2))


def test_diofant__functions__special__error_functions__erfcinv():
    assert _test_args(erfcinv(2))


def test_diofant__functions__special__error_functions__erf2inv():
    assert _test_args(erf2inv(2, 3))


def test_diofant__functions__special__error_functions__FresnelIntegral():
    pass


def test_diofant__functions__special__error_functions__fresnels():
    assert _test_args(fresnels(2))


def test_diofant__functions__special__error_functions__fresnelc():
    assert _test_args(fresnelc(2))


def test_diofant__functions__special__error_functions__erfs():
    assert _test_args(_erfs(2))


def test_diofant__functions__special__error_functions__Ei():
    assert _test_args(Ei(2))


def test_diofant__functions__special__error_functions__li():
    assert _test_args(li(2))


def test_diofant__functions__special__error_functions__Li():
    assert _test_args(Li(2))


def test_diofant__functions__special__error_functions__TrigonometricIntegral():
    pass


def test_diofant__functions__special__error_functions__Si():
    assert _test_args(Si(2))


def test_diofant__functions__special__error_functions__Ci():
    assert _test_args(Ci(2))


def test_diofant__functions__special__error_functions__Shi():
    assert _test_args(Shi(2))


def test_diofant__functions__special__error_functions__Chi():
    assert _test_args(Chi(2))


def test_diofant__functions__special__error_functions__expint():
    assert _test_args(expint(y, x))


def test_diofant__functions__special__gamma_functions__gamma():
    assert _test_args(gamma(x))


def test_diofant__functions__special__gamma_functions__loggamma():
    assert _test_args(loggamma(2))


def test_diofant__functions__special__gamma_functions__lowergamma():
    assert _test_args(lowergamma(x, 2))


def test_diofant__functions__special__gamma_functions__polygamma():
    assert _test_args(polygamma(x, 2))


def test_diofant__functions__special__gamma_functions__uppergamma():
    assert _test_args(uppergamma(x, 2))


def test_diofant__functions__special__beta_functions__beta():
    assert _test_args(beta(x, x))


def test_diofant__functions__special__hyper__TupleParametersBase():
    pass


def test_diofant__functions__special__hyper__TupleArg():
    pass


def test_diofant__functions__special__hyper__hyper():
    assert _test_args(hyper([1, 2, 3], [4, 5], x))


def test_diofant__functions__special__hyper__meijerg():
    assert _test_args(meijerg([1, 2, 3], [4, 5], [6], [], x))


def test_diofant__functions__special__hyper__HyperRep():
    pass


def test_diofant__functions__special__hyper__HyperRep_power1():
    assert _test_args(HyperRep_power1(x, y))


def test_diofant__functions__special__hyper__HyperRep_power2():
    assert _test_args(HyperRep_power2(x, y))


def test_diofant__functions__special__hyper__HyperRep_log1():
    assert _test_args(HyperRep_log1(x))


def test_diofant__functions__special__hyper__HyperRep_atanh():
    assert _test_args(HyperRep_atanh(x))


def test_diofant__functions__special__hyper__HyperRep_asin1():
    assert _test_args(HyperRep_asin1(x))


def test_diofant__functions__special__hyper__HyperRep_asin2():
    assert _test_args(HyperRep_asin2(x))


def test_diofant__functions__special__hyper__HyperRep_sqrts1():
    assert _test_args(HyperRep_sqrts1(x, y))


def test_diofant__functions__special__hyper__HyperRep_sqrts2():
    assert _test_args(HyperRep_sqrts2(x, y))


def test_diofant__functions__special__hyper__HyperRep_log2():
    assert _test_args(HyperRep_log2(x))


def test_diofant__functions__special__hyper__HyperRep_cosasin():
    assert _test_args(HyperRep_cosasin(x, y))


def test_diofant__functions__special__hyper__HyperRep_sinasin():
    assert _test_args(HyperRep_sinasin(x, y))


def test_diofant__functions__special__polynomials__OrthogonalPolynomial():
    pass


def test_diofant__functions__special__polynomials__jacobi():
    assert _test_args(jacobi(x, 2, 2, 2))


def test_diofant__functions__special__polynomials__gegenbauer():
    assert _test_args(gegenbauer(x, 2, 2))


def test_diofant__functions__special__polynomials__chebyshevt():
    assert _test_args(chebyshevt(x, 2))


def test_diofant__functions__special__polynomials__chebyshevt_root():
    assert _test_args(chebyshevt_root(3, 2))


def test_diofant__functions__special__polynomials__chebyshevu():
    assert _test_args(chebyshevu(x, 2))


def test_diofant__functions__special__polynomials__chebyshevu_root():
    assert _test_args(chebyshevu_root(3, 2))


def test_diofant__functions__special__polynomials__hermite():
    assert _test_args(hermite(x, 2))


def test_diofant__functions__special__polynomials__legendre():
    assert _test_args(legendre(x, 2))


def test_diofant__functions__special__polynomials__assoc_legendre():
    assert _test_args(assoc_legendre(x, 0, y))


def test_diofant__functions__special__polynomials__laguerre():
    assert _test_args(laguerre(x, 2))


def test_diofant__functions__special__polynomials__assoc_laguerre():
    assert _test_args(assoc_laguerre(x, 0, y))


def test_diofant__functions__special__spherical_harmonics__Ynm():
    assert _test_args(Ynm(1, 1, x, y))


def test_diofant__functions__special__spherical_harmonics__Znm():
    assert _test_args(Znm(1, 1, x, y))


def test_diofant__functions__special__tensor_functions__LeviCivita():
    assert _test_args(LeviCivita(x, y, 2))


def test_diofant__functions__special__tensor_functions__KroneckerDelta():
    assert _test_args(KroneckerDelta(x, y))


def test_diofant__functions__special__zeta_functions__dirichlet_eta():
    assert _test_args(dirichlet_eta(x))


def test_diofant__functions__special__zeta_functions__zeta():
    assert _test_args(zeta(101))


def test_diofant__functions__special__zeta_functions__lerchphi():
    assert _test_args(lerchphi(x, y, z))


def test_diofant__functions__special__zeta_functions__polylog():
    assert _test_args(polylog(x, y))


def test_diofant__integrals__integrals__Integral():
    assert _test_args(Integral(2, (x, 0, 1)))


def test_diofant__integrals__risch__NonElementaryIntegral():
    assert _test_args(NonElementaryIntegral(exp(-x**2), x))


def test_diofant__integrals__transforms__IntegralTransform():
    pass


def test_diofant__integrals__transforms__MellinTransform():
    assert _test_args(MellinTransform(2, x, y))


def test_diofant__integrals__transforms__InverseMellinTransform():
    assert _test_args(InverseMellinTransform(2, x, y, 0, 1))


def test_diofant__integrals__transforms__LaplaceTransform():
    assert _test_args(LaplaceTransform(2, x, y))


def test_diofant__integrals__transforms__InverseLaplaceTransform():
    assert _test_args(InverseLaplaceTransform(2, x, y, 0))


def test_diofant__integrals__transforms__FourierTypeTransform():
    pass


def test_diofant__integrals__transforms__InverseFourierTransform():
    assert _test_args(InverseFourierTransform(2, x, y))


def test_diofant__integrals__transforms__FourierTransform():
    assert _test_args(FourierTransform(2, x, y))


def test_diofant__integrals__transforms__SineCosineTypeTransform():
    pass


def test_diofant__integrals__transforms__InverseSineTransform():
    assert _test_args(InverseSineTransform(2, x, y))


def test_diofant__integrals__transforms__SineTransform():
    assert _test_args(SineTransform(2, x, y))


def test_diofant__integrals__transforms__InverseCosineTransform():
    assert _test_args(InverseCosineTransform(2, x, y))


def test_diofant__integrals__transforms__CosineTransform():
    assert _test_args(CosineTransform(2, x, y))


def test_diofant__integrals__transforms__HankelTypeTransform():
    pass


def test_diofant__integrals__transforms__InverseHankelTransform():
    assert _test_args(InverseHankelTransform(2, x, y, 0))


def test_diofant__integrals__transforms__HankelTransform():
    assert _test_args(HankelTransform(2, x, y, 0))


def test_diofant__logic__boolalg__And():
    assert _test_args(And(x, y, 2))


def test_diofant__logic__boolalg__Boolean():
    pass


def test_diofant__logic__boolalg__BooleanFunction():
    assert _test_args(BooleanFunction(1, 2, 3))


def test_diofant__logic__boolalg__BooleanAtom():
    pass


def test_diofant__logic__boolalg__BooleanTrue():
    assert _test_args(true)


def test_diofant__logic__boolalg__BooleanFalse():
    assert _test_args(false)


def test_diofant__logic__boolalg__Equivalent():
    assert _test_args(Equivalent(x, 2))


def test_diofant__logic__boolalg__ITE():
    assert _test_args(ITE(x, y, 2))


def test_diofant__logic__boolalg__Implies():
    assert _test_args(Implies(x, y))


def test_diofant__logic__boolalg__Nand():
    assert _test_args(Nand(x, y, 2))


def test_diofant__logic__boolalg__Nor():
    assert _test_args(Nor(x, y))


def test_diofant__logic__boolalg__Not():
    assert _test_args(Not(x))


def test_diofant__logic__boolalg__Or():
    assert _test_args(Or(x, y))


def test_diofant__logic__boolalg__Xor():
    assert _test_args(Xor(x, y, 2))


def test_diofant__matrices__expressions__matexpr__MatrixBase():
    pass


def test_diofant__matrices__immutable__ImmutableMatrix():
    m = ImmutableMatrix([[1, 2], [3, 4]])
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))
    m = ImmutableMatrix(1, 1, [1])
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))
    m = ImmutableMatrix(2, 2, lambda i, j: 1)
    assert m[0, 0] is Integer(1)
    m = ImmutableMatrix(2, 2, lambda i, j: 1/(1 + i) + 1/(1 + j))
    assert m[1, 1] is Integer(1)  # true div. will give 1.0 if i,j not sympified
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))


def test_diofant__matrices__immutable__ImmutableSparseMatrix():
    m = ImmutableSparseMatrix([[1, 2], [3, 4]])
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))
    m = ImmutableSparseMatrix(1, 1, {(0, 0): 1})
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))
    m = ImmutableSparseMatrix(1, 1, [1])
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))
    m = ImmutableSparseMatrix(2, 2, lambda i, j: 1)
    assert m[0, 0] is Integer(1)
    m = ImmutableSparseMatrix(2, 2, lambda i, j: 1/(1 + i) + 1/(1 + j))
    assert m[1, 1] is Integer(1)  # true div. will give 1.0 if i,j not sympified
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))


def test_diofant__matrices__expressions__slice__MatrixSlice():
    X = MatrixSymbol('X', 4, 4)
    assert _test_args(MatrixSlice(X, (0, 2), (0, 2)))


def test_diofant__matrices__expressions__blockmatrix__BlockDiagMatrix():
    X = MatrixSymbol('X', x, x)
    Y = MatrixSymbol('Y', y, y)
    assert _test_args(BlockDiagMatrix(X, Y))


def test_diofant__matrices__expressions__blockmatrix__BlockMatrix():
    X = MatrixSymbol('X', x, x)
    Y = MatrixSymbol('Y', y, y)
    Z = MatrixSymbol('Z', x, y)
    O = ZeroMatrix(y, x)
    assert _test_args(BlockMatrix([[X, Z], [O, Y]]))


def test_diofant__matrices__expressions__inverse__Inverse():
    assert _test_args(Inverse(MatrixSymbol('A', 3, 3)))


def test_diofant__matrices__expressions__matadd__MatAdd():
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', x, y)
    assert _test_args(MatAdd(X, Y))


def test_diofant__matrices__expressions__matexpr__Identity():
    assert _test_args(Identity(3))


def test_diofant__matrices__expressions__matexpr__MatrixExpr():
    pass


def test_diofant__matrices__expressions__matexpr__MatrixElement():
    assert _test_args(MatrixElement(MatrixSymbol('A', 3, 5), Integer(2), Integer(3)))


def test_diofant__matrices__expressions__matexpr__MatrixSymbol():
    assert _test_args(MatrixSymbol('A', 3, 5))


def test_diofant__matrices__expressions__matexpr__ZeroMatrix():
    assert _test_args(ZeroMatrix(3, 5))


def test_diofant__matrices__expressions__matmul__MatMul():
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', y, x)
    assert _test_args(MatMul(X, Y))


def test_diofant__matrices__expressions__diagonal__DiagonalMatrix():
    x = MatrixSymbol('x', 10, 1)
    assert _test_args(DiagonalMatrix(x))


def test_diofant__matrices__expressions__diagonal__DiagonalOf():
    X = MatrixSymbol('x', 10, 10)
    assert _test_args(DiagonalOf(X))


def test_diofant__matrices__expressions__hadamard__HadamardProduct():
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', x, y)
    assert _test_args(HadamardProduct(X, Y))


def test_diofant__matrices__expressions__matpow__MatPow():
    X = MatrixSymbol('X', x, x)
    assert _test_args(MatPow(X, 2))


def test_diofant__matrices__expressions__transpose__Transpose():
    assert _test_args(Transpose(MatrixSymbol('A', 3, 5)))


def test_diofant__matrices__expressions__adjoint__Adjoint():
    assert _test_args(Adjoint(MatrixSymbol('A', 3, 5)))


def test_diofant__matrices__expressions__trace__Trace():
    assert _test_args(Trace(MatrixSymbol('A', 3, 3)))


def test_diofant__matrices__expressions__determinant__Determinant():
    assert _test_args(Determinant(MatrixSymbol('A', 3, 3)))


def test_diofant__matrices__expressions__funcmatrix__FunctionMatrix():
    i, j = symbols('i,j')
    assert _test_args(FunctionMatrix(3, 3, Lambda((i, j), i - j) ))


def test_diofant__matrices__expressions__fourier__DFT():
    assert _test_args(DFT(Integer(2)))


def test_diofant__matrices__expressions__fourier__IDFT():
    assert _test_args(IDFT(Integer(2)))


def test_diofant__polys__polytools__GroebnerBasis():
    assert _test_args(GroebnerBasis([x, y, z], x, y, z))


def test_diofant__polys__polytools__Poly():
    assert _test_args(Poly(2, x, y))


def test_diofant__polys__polytools__PurePoly():
    assert _test_args(PurePoly(2, x, y))


def test_diofant__polys__rootoftools__RootOf():
    assert _test_args(RootOf(x**3 + x + 1, 0))


def test_diofant__polys__rootoftools__RootSum():
    assert _test_args(RootSum(x**3 + x + 1, sin))


def test_diofant__series__limits__Limit():
    assert _test_args(Limit(x, x, 0, dir='-'))


def test_diofant__series__order__Order():
    assert _test_args(Order(1, x, y))


def test_diofant__simplify__hyperexpand__Hyper_Function():
    assert _test_args(Hyper_Function([2], [1]))


def test_diofant__simplify__hyperexpand__G_Function():
    assert _test_args(G_Function([2], [1], [], []))


def test_diofant__tensor__array__ndim_array__ImmutableNDimArray():
    pass


def test_diofant__tensor__array__dense_ndim_array__ImmutableDenseNDimArray():
    densarr = ImmutableDenseNDimArray(range(10, 34), (2, 3, 4))
    assert _test_args(densarr)


def test_diofant__tensor__array__sparse_ndim_array__ImmutableSparseNDimArray():
    sparr = ImmutableSparseNDimArray(range(10, 34), (2, 3, 4))
    assert _test_args(sparr)


def test_diofant__tensor__indexed__Idx():
    assert _test_args(Idx('test'))
    assert _test_args(Idx(1, (0, 10)))


def test_diofant__tensor__indexed__Indexed():
    assert _test_args(Indexed('A', Idx('i'), Idx('j')))


def test_diofant__tensor__indexed__IndexedBase():
    assert _test_args(IndexedBase('A', shape=(x, y)))
    assert _test_args(IndexedBase('A', 1))
    assert _test_args(IndexedBase('A')[0, 1])


def test_diofant__tensor__tensor__TensorIndexType():
    assert _test_args(TensorIndexType('Lorentz', metric=False))


def test_diofant__tensor__tensor__TensorSymmetry():
    assert _test_args(TensorSymmetry(get_symmetric_group_sgs(2)))


def test_diofant__tensor__tensor__TensorType():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    assert _test_args(TensorType([Lorentz], sym))


def test_diofant__tensor__tensor__TensorHead():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    assert _test_args(TensorHead('p', S1, 0))


def test_diofant__tensor__tensor__TensorIndex():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    assert _test_args(TensorIndex('i', Lorentz))


def test_diofant__tensor__tensor__TensExpr():
    pass


def test_diofant__tensor__tensor__TensAdd():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p, q = S1('p,q')
    t1 = p(a)
    t2 = q(a)
    assert _test_args(TensAdd(t1, t2))


def test_diofant__tensor__tensor__Tensor():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p = S1('p')
    assert _test_args(p(a))


def test_diofant__tensor__tensor__TensMul():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p = S1('p')
    q = S1('q')
    assert _test_args(3*p(a)*q(b))


def test_as_coeff_add():
    assert (7, (3*x, 4*x**2)) == (7 + 3*x + 4*x**2).as_coeff_add()


def test_diofant__geometry__curve__Curve():
    assert _test_args(Curve((x, 1), (x, 0, 1)))


def test_diofant__geometry__point__Point():
    assert _test_args(Point(0, 1))


def test_diofant__geometry__ellipse__Ellipse():
    assert _test_args(Ellipse((0, 1), 2, 3))


def test_diofant__geometry__ellipse__Circle():
    assert _test_args(Circle((0, 1), 2))


def test_diofant__geometry__line__LinearEntity():
    pass


def test_diofant__geometry__line__Line():
    assert _test_args(Line((0, 1), (2, 3)))


def test_diofant__geometry__line__Ray():
    assert _test_args(Ray((0, 1), (2, 3)))


def test_diofant__geometry__line__Segment():
    assert _test_args(Segment((0, 1), (2, 3)))


def test_diofant__geometry__polygon__Polygon():
    assert _test_args(Polygon((0, 1), (2, 3), (4, 5), (6, 7)))


def test_diofant__geometry__polygon__RegularPolygon():
    assert _test_args(RegularPolygon((0, 1), 2, 3, 4))


def test_diofant__geometry__polygon__Triangle():
    assert _test_args(Triangle((0, 1), (2, 3), (4, 5)))


def test_diofant__geometry__entity__GeometryEntity():
    assert _test_args(GeometryEntity(Point(1, 0), 1, [1, 2]))


def test_diofant__geometry__entity__GeometrySet():
    pass


def test_diofant__diffgeom__diffgeom__Manifold():
    assert _test_args(Manifold('name', 3))


def test_diofant__diffgeom__diffgeom__Patch():
    assert _test_args(Patch('name', Manifold('name', 3)))


def test_diofant__diffgeom__diffgeom__CoordSystem():
    assert _test_args(CoordSystem('name', Patch('name', Manifold('name', 3))))


def test_diofant__diffgeom__diffgeom__Point():
    assert _test_args(DiffgeomPoint(
        CoordSystem('name', Patch('name', Manifold('name', 3))), [x, y]))


def test_diofant__diffgeom__diffgeom__BaseScalarField():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseScalarField(cs, 0))


def test_diofant__diffgeom__diffgeom__BaseVectorField():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseVectorField(cs, 0))


def test_diofant__diffgeom__diffgeom__Differential():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(Differential(BaseScalarField(cs, 0)))


def test_diofant__diffgeom__diffgeom__Commutator():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    cs1 = CoordSystem('name1', Patch('name', Manifold('name', 3)))
    v = BaseVectorField(cs, 0)
    v1 = BaseVectorField(cs1, 0)
    assert _test_args(Commutator(v, v1))


def test_diofant__diffgeom__diffgeom__TensorProduct():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    assert _test_args(TensorProduct(d, d))


def test_diofant__diffgeom__diffgeom__WedgeProduct():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    d1 = Differential(BaseScalarField(cs, 1))
    assert _test_args(WedgeProduct(d, d1))


def test_diofant__diffgeom__diffgeom__LieDerivative():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    v = BaseVectorField(cs, 0)
    assert _test_args(LieDerivative(v, d))


def test_diofant__diffgeom__diffgeom__BaseCovarDerivativeOp():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseCovarDerivativeOp(cs, 0, [[[0, ]*3, ]*3, ]*3))


def test_diofant__diffgeom__diffgeom__CovarDerivativeOp():
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    v = BaseVectorField(cs, 0)
    _test_args(CovarDerivativeOp(v, [[[0, ]*3, ]*3, ]*3))


def test_diofant__ntheory__factor___totient():
    k = symbols('k', integer=True)
    t = totient(k)
    assert _test_args(t)


def test_diofant__ntheory__factor___divisor_sigma():
    k = symbols('k', integer=True)
    n = symbols('n', integer=True)
    t = divisor_sigma(n, k)
    assert _test_args(t)


def test_diofant__ntheory__residue_ntheory__mobius():
    assert _test_args(mobius(2))


def test_diofant__printing__codeprinter__Assignment():
    assert _test_args(Assignment(x, y))


def test_diofant__vector__coordsysrect__CoordSysCartesian():
    assert _test_args(CoordSysCartesian('C'))


def test_diofant__vector__point__Point():
    assert _test_args(VPoint('P'))


def test_diofant__vector__basisdependent__BasisDependent():
    # from diofant.vector.basisdependent import BasisDependent
    pass
    # These classes have been created to maintain an OOP hierarchy
    # for Vectors and Dyadics. Are NOT meant to be initialized


def test_diofant__vector__basisdependent__BasisDependentMul():
    # from diofant.vector.basisdependent import BasisDependentMul
    pass
    # These classes have been created to maintain an OOP hierarchy
    # for Vectors and Dyadics. Are NOT meant to be initialized


def test_diofant__vector__basisdependent__BasisDependentAdd():
    # from diofant.vector.basisdependent import BasisDependentAdd
    pass
    # These classes have been created to maintain an OOP hierarchy
    # for Vectors and Dyadics. Are NOT meant to be initialized


def test_diofant__vector__basisdependent__BasisDependentZero():
    # from diofant.vector.basisdependent import BasisDependentZero
    pass
    # These classes have been created to maintain an OOP hierarchy
    # for Vectors and Dyadics. Are NOT meant to be initialized


def test_diofant__vector__vector__BaseVector():
    C = CoordSysCartesian('C')
    assert _test_args(BaseVector('Ci', 0, C, ' ', ' '))


def test_diofant__vector__vector__VectorAdd():
    C = CoordSysCartesian('C')
    v1 = a*C.i + b*C.j + c*C.k
    v2 = x*C.i + y*C.j + z*C.k
    assert _test_args(VectorAdd(v1, v2))
    assert _test_args(VectorMul(x, v1))


def test_diofant__vector__vector__VectorMul():
    C = CoordSysCartesian('C')
    assert _test_args(VectorMul(a, C.i))


def test_diofant__vector__vector__VectorZero():
    assert _test_args(VectorZero())


def test_diofant__vector__vector__Vector():
    # from diofant.vector.vector import Vector
    # Vector is never to be initialized using args
    pass


def test_diofant__vector__dyadic__Dyadic():
    # from diofant.vector.dyadic import Dyadic
    # Dyadic is never to be initialized using args
    pass


def test_diofant__vector__dyadic__BaseDyadic():
    C = CoordSysCartesian('C')
    assert _test_args(BaseDyadic(C.i, C.j))


def test_diofant__vector__dyadic__DyadicMul():
    C = CoordSysCartesian('C')
    assert _test_args(DyadicMul(3, BaseDyadic(C.i, C.j)))


def test_diofant__vector__dyadic__DyadicAdd():
    C = CoordSysCartesian('C')
    assert _test_args(2 * DyadicAdd(BaseDyadic(C.i, C.i),
                                    BaseDyadic(C.i, C.j)))


def test_diofant__vector__dyadic__DyadicZero():
    assert _test_args(DyadicZero())


def test_diofant__vector__deloperator__Del():
    C = CoordSysCartesian('C')
    assert _test_args(Del(C))


def test_diofant__vector__orienters__Orienter():
    # from diofant.vector.orienters import Orienter
    pass
    # Not to be initialized


def test_diofant__vector__orienters__ThreeAngleOrienter():
    # from diofant.vector.orienters import ThreeAngleOrienter
    pass
    # Not to be initialized


def test_diofant__vector__orienters__AxisOrienter():
    C = CoordSysCartesian('C')
    assert _test_args(AxisOrienter(x, C.i))


def test_diofant__vector__orienters__BodyOrienter():
    assert _test_args(BodyOrienter(x, y, z, '123'))


def test_diofant__vector__orienters__SpaceOrienter():
    assert _test_args(SpaceOrienter(x, y, z, '123'))


def test_diofant__vector__orienters__QuaternionOrienter():
    a, b, c, d = symbols('a b c d')
    assert _test_args(QuaternionOrienter(a, b, c, d))


def test_diofant__vector__scalar__BaseScalar():
    C = CoordSysCartesian('C')
    assert _test_args(BaseScalar('Cx', 0, C, ' ', ' '))
