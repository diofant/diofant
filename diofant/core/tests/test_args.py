"""Test whether all elements of cls.args are instances of Basic. """

# NOTE: keep tests sorted by (module, class name) key.

import os
import re
import warnings
import io
import inspect

import pytest

from diofant import Basic, S, symbols, sqrt, sin, oo, Interval, exp, Integer
from diofant.utilities.exceptions import DiofantDeprecationWarning

from diofant.abc import x, y, z


def test_all_classes_are_tested():
    this = os.path.split(__file__)[0]
    path = os.path.join(this, os.pardir, os.pardir)
    diofant_path = os.path.abspath(path)
    prefix = os.path.split(diofant_path)[0] + os.sep

    re_cls = re.compile(r"^class ([A-Za-z][A-Za-z0-9_]*)\s*\(", re.MULTILINE)

    modules = {}

    for root, dirs, files in os.walk(diofant_path):
        module = root.replace(prefix, "").replace(os.sep, ".")

        for file in files:
            if file.startswith(("_", "test_", "bench_")):
                continue
            if not file.endswith(".py"):
                continue

            with io.open(os.path.join(root, file), "r", encoding='utf-8') as f:
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
    warnings.simplefilter("error", category=DiofantDeprecationWarning)

    assert not failed, "Missing classes: %s.  Please add tests for these to diofant/core/tests/test_args.py." % ", ".join(failed)


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


@pytest.mark.xfail
def test_diofant__combinatorics__graycode__GrayCode():
    from diofant.combinatorics.graycode import GrayCode
    # an integer is given and returned from GrayCode as the arg
    assert _test_args(GrayCode(3, start='100'))
    assert _test_args(GrayCode(3, rank=1))


def test_diofant__combinatorics__subsets__Subset():
    from diofant.combinatorics.subsets import Subset
    assert _test_args(Subset([0, 1], [0, 1, 2, 3]))
    assert _test_args(Subset(['c', 'd'], ['a', 'b', 'c', 'd']))


@pytest.mark.xfail
def test_diofant__combinatorics__permutations__Permutation():
    from diofant.combinatorics.permutations import Permutation
    assert _test_args(Permutation([0, 1, 2, 3]))


def test_diofant__combinatorics__perm_groups__PermutationGroup():
    from diofant.combinatorics.permutations import Permutation
    from diofant.combinatorics.perm_groups import PermutationGroup
    assert _test_args(PermutationGroup([Permutation([0, 1])]))


def test_diofant__combinatorics__polyhedron__Polyhedron():
    from diofant.combinatorics.permutations import Permutation
    from diofant.combinatorics.polyhedron import Polyhedron
    from diofant.abc import w, x, y, z
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


@pytest.mark.xfail
def test_diofant__combinatorics__prufer__Prufer():
    from diofant.combinatorics.prufer import Prufer
    assert _test_args(Prufer([[0, 1], [0, 2], [0, 3]], 4))


def test_diofant__combinatorics__partitions__Partition():
    from diofant.combinatorics.partitions import Partition
    assert _test_args(Partition([1]))


@pytest.mark.xfail
def test_diofant__combinatorics__partitions__IntegerPartition():
    from diofant.combinatorics.partitions import IntegerPartition
    assert _test_args(IntegerPartition([1]))


def test_diofant__concrete__products__Product():
    from diofant.concrete.products import Product
    assert _test_args(Product(x, (x, 0, 10)))
    assert _test_args(Product(x, (x, 0, y), (y, 0, 10)))


def test_diofant__concrete__expr_with_limits__ExprWithLimits():
    from diofant.concrete.expr_with_limits import ExprWithLimits
    assert _test_args(ExprWithLimits(x, (x, 0, 10)))
    assert _test_args(ExprWithLimits(x*y, (x, 0, 10.), (y, 1., 3)))


def test_diofant__concrete__expr_with_limits__AddWithLimits():
    from diofant.concrete.expr_with_limits import AddWithLimits
    assert _test_args(AddWithLimits(x, (x, 0, 10)))
    assert _test_args(AddWithLimits(x*y, (x, 0, 10), (y, 1, 3)))


def test_diofant__concrete__expr_with_intlimits__ExprWithIntLimits():
    from diofant.concrete.expr_with_intlimits import ExprWithIntLimits
    assert _test_args(ExprWithIntLimits(x, (x, 0, 10)))
    assert _test_args(ExprWithIntLimits(x*y, (x, 0, 10), (y, 1, 3)))


def test_diofant__concrete__summations__Sum():
    from diofant.concrete.summations import Sum
    assert _test_args(Sum(x, (x, 0, 10)))
    assert _test_args(Sum(x, (x, 0, y), (y, 0, 10)))


def test_diofant__core__add__Add():
    from diofant.core.add import Add
    assert _test_args(Add(x, y, z, 2))


def test_diofant__core__basic__Atom():
    from diofant.core.basic import Atom
    assert _test_args(Atom())


def test_diofant__core__basic__Basic():
    from diofant.core.basic import Basic
    assert _test_args(Basic())


def test_diofant__core__containers__Dict():
    from diofant.core.containers import Dict
    assert _test_args(Dict({x: y, y: z}))


def test_diofant__core__containers__Tuple():
    from diofant.core.containers import Tuple
    assert _test_args(Tuple(x, y, z, 2))


def test_diofant__core__expr__AtomicExpr():
    from diofant.core.expr import AtomicExpr
    assert _test_args(AtomicExpr())


def test_diofant__core__expr__Expr():
    from diofant.core.expr import Expr
    assert _test_args(Expr())


def test_diofant__core__function__Application():
    from diofant.core.function import Application
    assert _test_args(Application(1, 2, 3))


def test_diofant__core__function__AppliedUndef():
    from diofant.core.function import AppliedUndef
    assert _test_args(AppliedUndef(1, 2, 3))


def test_diofant__core__function__Derivative():
    from diofant.core.function import Derivative
    assert _test_args(Derivative(2, x, y, 3))


def test_diofant__core__function__Function():
    pass


def test_diofant__core__function__Lambda():
    from diofant.core.function import Lambda
    assert _test_args(Lambda((x, y), x + y + z))


def test_diofant__core__function__Subs():
    from diofant.core.function import Subs
    assert _test_args(Subs(x + y, x, 2))


def test_diofant__core__function__WildFunction():
    from diofant.core.function import WildFunction
    assert _test_args(WildFunction('f'))


def test_diofant__core__mod__Mod():
    from diofant.core.mod import Mod
    assert _test_args(Mod(x, 2))


def test_diofant__core__mul__Mul():
    from diofant.core.mul import Mul
    assert _test_args(Mul(2, x, y, z))


def test_diofant__core__numbers__Catalan():
    from diofant.core.numbers import Catalan
    assert _test_args(Catalan())


def test_diofant__core__numbers__ComplexInfinity():
    from diofant.core.numbers import ComplexInfinity
    assert _test_args(ComplexInfinity())


def test_diofant__core__numbers__EulerGamma():
    from diofant.core.numbers import EulerGamma
    assert _test_args(EulerGamma())


def test_diofant__core__numbers__Exp1():
    from diofant.core.numbers import Exp1
    assert _test_args(Exp1())


def test_diofant__core__numbers__Float():
    from diofant.core.numbers import Float
    assert _test_args(Float(1.23))


def test_diofant__core__numbers__GoldenRatio():
    from diofant.core.numbers import GoldenRatio
    assert _test_args(GoldenRatio())


def test_diofant__core__numbers__Half():
    from diofant.core.numbers import Half
    assert _test_args(Half())


def test_diofant__core__numbers__ImaginaryUnit():
    from diofant.core.numbers import ImaginaryUnit
    assert _test_args(ImaginaryUnit())


def test_diofant__core__numbers__Infinity():
    from diofant.core.numbers import Infinity
    assert _test_args(Infinity())


def test_diofant__core__numbers__Integer():
    from diofant.core.numbers import Integer
    assert _test_args(Integer(7))


def test_diofant__core__numbers__IntegerConstant():
    pass


def test_diofant__core__numbers__NaN():
    from diofant.core.numbers import NaN
    assert _test_args(NaN())


def test_diofant__core__numbers__NegativeInfinity():
    from diofant.core.numbers import NegativeInfinity
    assert _test_args(NegativeInfinity())


def test_diofant__core__numbers__NegativeOne():
    from diofant.core.numbers import NegativeOne
    assert _test_args(NegativeOne())


def test_diofant__core__numbers__Number():
    from diofant.core.numbers import Number
    assert _test_args(Number(1, 7))


def test_diofant__core__numbers__NumberSymbol():
    from diofant.core.numbers import NumberSymbol
    assert _test_args(NumberSymbol())


def test_diofant__core__numbers__One():
    from diofant.core.numbers import One
    assert _test_args(One())


def test_diofant__core__numbers__Pi():
    from diofant.core.numbers import Pi
    assert _test_args(Pi())


def test_diofant__core__numbers__Rational():
    from diofant.core.numbers import Rational
    assert _test_args(Rational(1, 7))


def test_diofant__core__numbers__RationalConstant():
    pass


def test_diofant__core__numbers__Zero():
    from diofant.core.numbers import Zero
    assert _test_args(Zero())


def test_diofant__core__operations__AssocOp():
    pass


def test_diofant__core__operations__LatticeOp():
    pass


def test_diofant__core__power__Pow():
    from diofant.core.power import Pow
    assert _test_args(Pow(x, 2))


def test_diofant__core__relational__Equality():
    from diofant.core.relational import Equality
    assert _test_args(Equality(x, 2))


def test_diofant__core__relational__GreaterThan():
    from diofant.core.relational import GreaterThan
    assert _test_args(GreaterThan(x, 2))


def test_diofant__core__relational__LessThan():
    from diofant.core.relational import LessThan
    assert _test_args(LessThan(x, 2))


def test_diofant__core__relational__Relational():
    pass


def test_diofant__core__relational__StrictGreaterThan():
    from diofant.core.relational import StrictGreaterThan
    assert _test_args(StrictGreaterThan(x, 2))


def test_diofant__core__relational__StrictLessThan():
    from diofant.core.relational import StrictLessThan
    assert _test_args(StrictLessThan(x, 2))


def test_diofant__core__relational__Unequality():
    from diofant.core.relational import Unequality
    assert _test_args(Unequality(x, 2))


def test_diofant__sets__sets__EmptySet():
    from diofant.sets.sets import EmptySet
    assert _test_args(EmptySet())


def test_diofant__sets__sets__UniversalSet():
    from diofant.sets.sets import UniversalSet
    assert _test_args(UniversalSet())


def test_diofant__sets__sets__FiniteSet():
    from diofant.sets.sets import FiniteSet
    assert _test_args(FiniteSet(x, y, z))


def test_diofant__sets__sets__Interval():
    from diofant.sets.sets import Interval
    assert _test_args(Interval(0, 1))


def test_diofant__sets__sets__ProductSet():
    from diofant.sets.sets import ProductSet, Interval
    assert _test_args(ProductSet(Interval(0, 1), Interval(0, 1)))


def test_diofant__sets__sets__Set():
    from diofant.sets.sets import Set
    assert _test_args(Set())


def test_diofant__sets__sets__Intersection():
    from diofant.sets.sets import Intersection, Interval
    assert _test_args(Intersection(Interval(0, 3), Interval(2, x)))


def test_diofant__sets__sets__Union():
    from diofant.sets.sets import Union, Interval
    assert _test_args(Union(Interval(0, 1), Interval(2, 3)))


def test_diofant__sets__sets__Complement():
    from diofant.sets.sets import Complement
    assert _test_args(Complement(Interval(0, 2), Interval(0, 1)))


def test_diofant__sets__sets__SymmetricDifference():
    from diofant.sets.sets import FiniteSet, SymmetricDifference
    assert _test_args(SymmetricDifference(FiniteSet(1, 2, 3),
           FiniteSet(2, 3, 4)))


def test_diofant__core__trace__Tr():
    from diofant.core.trace import Tr
    a, b = symbols('a b')
    assert _test_args(Tr(a + b))


def test_diofant__sets__fancysets__Naturals():
    from diofant.sets.fancysets import Naturals
    assert _test_args(Naturals())


def test_diofant__sets__fancysets__Naturals0():
    from diofant.sets.fancysets import Naturals0
    assert _test_args(Naturals0())


def test_diofant__sets__fancysets__Integers():
    from diofant.sets.fancysets import Integers
    assert _test_args(Integers())


def test_diofant__sets__fancysets__Rationals():
    from diofant.sets.fancysets import Rationals
    assert _test_args(Rationals())


def test_diofant__sets__fancysets__Reals():
    from diofant.sets.fancysets import Reals
    assert _test_args(Reals())


def test_diofant__sets__fancysets__ImageSet():
    from diofant.sets.fancysets import ImageSet
    from diofant import S, Lambda, Symbol
    x = Symbol('x')
    assert _test_args(ImageSet(Lambda(x, x**2), S.Naturals))


def test_diofant__sets__fancysets__Range():
    from diofant.sets.fancysets import Range
    assert _test_args(Range(1, 5, 1))


def test_diofant__sets__contains__Contains():
    from diofant.sets.fancysets import Range
    from diofant.sets.contains import Contains
    assert _test_args(Contains(x, Range(0, 10, 2)))


# STATS


from diofant.stats.crv_types import NormalDistribution
nd = NormalDistribution(0, 1)
from diofant.stats.frv_types import DieDistribution
die = DieDistribution(6)


def test_diofant__stats__crv__ContinuousDomain():
    from diofant.stats.crv import ContinuousDomain
    assert _test_args(ContinuousDomain({x}, Interval(-oo, oo)))


def test_diofant__stats__crv__SingleContinuousDomain():
    from diofant.stats.crv import SingleContinuousDomain
    assert _test_args(SingleContinuousDomain(x, Interval(-oo, oo)))


def test_diofant__stats__crv__ProductContinuousDomain():
    from diofant.stats.crv import SingleContinuousDomain, ProductContinuousDomain
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    E = SingleContinuousDomain(y, Interval(0, oo))
    assert _test_args(ProductContinuousDomain(D, E))


def test_diofant__stats__crv__ConditionalContinuousDomain():
    from diofant.stats.crv import (SingleContinuousDomain,
            ConditionalContinuousDomain)
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    assert _test_args(ConditionalContinuousDomain(D, x > 0))


def test_diofant__stats__crv__ContinuousPSpace():
    from diofant.stats.crv import ContinuousPSpace, SingleContinuousDomain
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    assert _test_args(ContinuousPSpace(D, nd))


def test_diofant__stats__crv__SingleContinuousPSpace():
    from diofant.stats.crv import SingleContinuousPSpace
    assert _test_args(SingleContinuousPSpace(x, nd))


def test_diofant__stats__crv__ProductContinuousPSpace():
    from diofant.stats.crv import ProductContinuousPSpace, SingleContinuousPSpace
    A = SingleContinuousPSpace(x, nd)
    B = SingleContinuousPSpace(y, nd)
    assert _test_args(ProductContinuousPSpace(A, B))


def test_diofant__stats__crv__SingleContinuousDistribution():
    pass


def test_diofant__stats__drv__SingleDiscreteDomain():
    from diofant.stats.drv import SingleDiscreteDomain
    assert _test_args(SingleDiscreteDomain(x, S.Naturals))


def test_diofant__stats__drv__SingleDiscretePSpace():
    from diofant.stats.drv import SingleDiscretePSpace
    from diofant.stats.drv_types import PoissonDistribution
    assert _test_args(SingleDiscretePSpace(x, PoissonDistribution(1)))


def test_diofant__stats__drv__SingleDiscreteDistribution():
    pass


def test_diofant__stats__rv__RandomDomain():
    from diofant.stats.rv import RandomDomain
    from diofant.sets.sets import FiniteSet
    assert _test_args(RandomDomain(FiniteSet(x), FiniteSet(1, 2, 3)))


def test_diofant__stats__rv__SingleDomain():
    from diofant.stats.rv import SingleDomain
    from diofant.sets.sets import FiniteSet
    assert _test_args(SingleDomain(x, FiniteSet(1, 2, 3)))


def test_diofant__stats__rv__ConditionalDomain():
    from diofant.stats.rv import ConditionalDomain, RandomDomain
    from diofant.sets.sets import FiniteSet
    D = RandomDomain(FiniteSet(x), FiniteSet(1, 2))
    assert _test_args(ConditionalDomain(D, x > 1))


def test_diofant__stats__rv__PSpace():
    from diofant.stats.rv import PSpace, RandomDomain
    from diofant import FiniteSet
    D = RandomDomain(FiniteSet(x), FiniteSet(1, 2, 3, 4, 5, 6))
    assert _test_args(PSpace(D, die))


def test_diofant__stats__rv__SinglePSpace():
    pass


def test_diofant__stats__rv__RandomSymbol():
    from diofant.stats.rv import RandomSymbol
    from diofant.stats.crv import SingleContinuousPSpace
    A = SingleContinuousPSpace(x, nd)
    assert _test_args(RandomSymbol(A, x))


def test_diofant__stats__rv__ProductPSpace():
    from diofant.stats.rv import ProductPSpace
    from diofant.stats.crv import SingleContinuousPSpace
    A = SingleContinuousPSpace(x, nd)
    B = SingleContinuousPSpace(y, nd)
    assert _test_args(ProductPSpace(A, B))


def test_diofant__stats__rv__ProductDomain():
    from diofant.stats.rv import ProductDomain, SingleDomain
    D = SingleDomain(x, Interval(-oo, oo))
    E = SingleDomain(y, Interval(0, oo))
    assert _test_args(ProductDomain(D, E))


def test_diofant__stats__frv_types__DiscreteUniformDistribution():
    from diofant.stats.frv_types import DiscreteUniformDistribution
    from diofant.core.containers import Tuple
    assert _test_args(DiscreteUniformDistribution(Tuple(*range(6))))


def test_diofant__stats__frv_types__DieDistribution():
    from diofant.stats.frv_types import DieDistribution
    assert _test_args(DieDistribution(6))


def test_diofant__stats__frv_types__BernoulliDistribution():
    from diofant.stats.frv_types import BernoulliDistribution
    assert _test_args(BernoulliDistribution(S.Half, 0, 1))


def test_diofant__stats__frv_types__BinomialDistribution():
    from diofant.stats.frv_types import BinomialDistribution
    assert _test_args(BinomialDistribution(5, S.Half, 1, 0))


def test_diofant__stats__frv_types__HypergeometricDistribution():
    from diofant.stats.frv_types import HypergeometricDistribution
    assert _test_args(HypergeometricDistribution(10, 5, 3))


def test_diofant__stats__frv_types__RademacherDistribution():
    from diofant.stats.frv_types import RademacherDistribution
    assert _test_args(RademacherDistribution())


def test_diofant__stats__frv__FiniteDomain():
    from diofant.stats.frv import FiniteDomain
    assert _test_args(FiniteDomain({(x, 1), (x, 2)}))  # x can be 1 or 2


def test_diofant__stats__frv__SingleFiniteDomain():
    from diofant.stats.frv import SingleFiniteDomain
    assert _test_args(SingleFiniteDomain(x, {1, 2}))  # x can be 1 or 2


def test_diofant__stats__frv__ProductFiniteDomain():
    from diofant.stats.frv import SingleFiniteDomain, ProductFiniteDomain
    xd = SingleFiniteDomain(x, {1, 2})
    yd = SingleFiniteDomain(y, {1, 2})
    assert _test_args(ProductFiniteDomain(xd, yd))


def test_diofant__stats__frv__ConditionalFiniteDomain():
    from diofant.stats.frv import SingleFiniteDomain, ConditionalFiniteDomain
    xd = SingleFiniteDomain(x, {1, 2})
    assert _test_args(ConditionalFiniteDomain(xd, x > 1))


def test_diofant__stats__frv__FinitePSpace():
    from diofant.stats.frv import FinitePSpace, SingleFiniteDomain
    xd = SingleFiniteDomain(x, {1, 2, 3, 4, 5, 6})
    p = 1.0/6
    xd = SingleFiniteDomain(x, {1, 2})
    assert _test_args(FinitePSpace(xd, {(x, 1): S.Half, (x, 2): S.Half}))


def test_diofant__stats__frv__SingleFinitePSpace():
    from diofant.stats.frv import SingleFinitePSpace
    from diofant import Symbol

    assert _test_args(SingleFinitePSpace(Symbol('x'), die))


def test_diofant__stats__frv__ProductFinitePSpace():
    from diofant.stats.frv import SingleFinitePSpace, ProductFinitePSpace
    from diofant import Symbol
    xp = SingleFinitePSpace(Symbol('x'), die)
    yp = SingleFinitePSpace(Symbol('y'), die)
    assert _test_args(ProductFinitePSpace(xp, yp))


def test_diofant__stats__frv__SingleFiniteDistribution():
    pass


def test_diofant__stats__crv__ContinuousDistribution():
    pass


def test_diofant__stats__frv_types__FiniteDistributionHandmade():
    from diofant.stats.frv_types import FiniteDistributionHandmade
    assert _test_args(FiniteDistributionHandmade({1: 1}))


def test_diofant__stats__crv__ContinuousDistributionHandmade():
    from diofant.stats.crv import ContinuousDistributionHandmade
    from diofant import Symbol, Interval
    assert _test_args(ContinuousDistributionHandmade(Symbol('x'),
                                                     Interval(0, 2)))


def test_diofant__stats__rv__Density():
    from diofant.stats.rv import Density
    from diofant.stats.crv_types import Normal
    assert _test_args(Density(Normal('x', 0, 1)))


def test_diofant__stats__crv_types__ArcsinDistribution():
    from diofant.stats.crv_types import ArcsinDistribution
    assert _test_args(ArcsinDistribution(0, 1))


def test_diofant__stats__crv_types__BeniniDistribution():
    from diofant.stats.crv_types import BeniniDistribution
    assert _test_args(BeniniDistribution(1, 1, 1))


def test_diofant__stats__crv_types__BetaDistribution():
    from diofant.stats.crv_types import BetaDistribution
    assert _test_args(BetaDistribution(1, 1))


def test_diofant__stats__crv_types__BetaPrimeDistribution():
    from diofant.stats.crv_types import BetaPrimeDistribution
    assert _test_args(BetaPrimeDistribution(1, 1))


def test_diofant__stats__crv_types__CauchyDistribution():
    from diofant.stats.crv_types import CauchyDistribution
    assert _test_args(CauchyDistribution(0, 1))


def test_diofant__stats__crv_types__ChiDistribution():
    from diofant.stats.crv_types import ChiDistribution
    assert _test_args(ChiDistribution(1))


def test_diofant__stats__crv_types__ChiNoncentralDistribution():
    from diofant.stats.crv_types import ChiNoncentralDistribution
    assert _test_args(ChiNoncentralDistribution(1, 1))


def test_diofant__stats__crv_types__ChiSquaredDistribution():
    from diofant.stats.crv_types import ChiSquaredDistribution
    assert _test_args(ChiSquaredDistribution(1))


def test_diofant__stats__crv_types__DagumDistribution():
    from diofant.stats.crv_types import DagumDistribution
    assert _test_args(DagumDistribution(1, 1, 1))


def test_diofant__stats__crv_types__ExponentialDistribution():
    from diofant.stats.crv_types import ExponentialDistribution
    assert _test_args(ExponentialDistribution(1))


def test_diofant__stats__crv_types__FDistributionDistribution():
    from diofant.stats.crv_types import FDistributionDistribution
    assert _test_args(FDistributionDistribution(1, 1))


def test_diofant__stats__crv_types__FisherZDistribution():
    from diofant.stats.crv_types import FisherZDistribution
    assert _test_args(FisherZDistribution(1, 1))


def test_diofant__stats__crv_types__FrechetDistribution():
    from diofant.stats.crv_types import FrechetDistribution
    assert _test_args(FrechetDistribution(1, 1, 1))


def test_diofant__stats__crv_types__GammaInverseDistribution():
    from diofant.stats.crv_types import GammaInverseDistribution
    assert _test_args(GammaInverseDistribution(1, 1))


def test_diofant__stats__crv_types__GammaDistribution():
    from diofant.stats.crv_types import GammaDistribution
    assert _test_args(GammaDistribution(1, 1))


def test_diofant__stats__crv_types__KumaraswamyDistribution():
    from diofant.stats.crv_types import KumaraswamyDistribution
    assert _test_args(KumaraswamyDistribution(1, 1))


def test_diofant__stats__crv_types__LaplaceDistribution():
    from diofant.stats.crv_types import LaplaceDistribution
    assert _test_args(LaplaceDistribution(0, 1))


def test_diofant__stats__crv_types__LogisticDistribution():
    from diofant.stats.crv_types import LogisticDistribution
    assert _test_args(LogisticDistribution(0, 1))


def test_diofant__stats__crv_types__LogNormalDistribution():
    from diofant.stats.crv_types import LogNormalDistribution
    assert _test_args(LogNormalDistribution(0, 1))


def test_diofant__stats__crv_types__MaxwellDistribution():
    from diofant.stats.crv_types import MaxwellDistribution
    assert _test_args(MaxwellDistribution(1))


def test_diofant__stats__crv_types__NakagamiDistribution():
    from diofant.stats.crv_types import NakagamiDistribution
    assert _test_args(NakagamiDistribution(1, 1))


def test_diofant__stats__crv_types__NormalDistribution():
    from diofant.stats.crv_types import NormalDistribution
    assert _test_args(NormalDistribution(0, 1))


def test_diofant__stats__crv_types__ParetoDistribution():
    from diofant.stats.crv_types import ParetoDistribution
    assert _test_args(ParetoDistribution(1, 1))


def test_diofant__stats__crv_types__QuadraticUDistribution():
    from diofant.stats.crv_types import QuadraticUDistribution
    assert _test_args(QuadraticUDistribution(1, 2))


def test_diofant__stats__crv_types__RaisedCosineDistribution():
    from diofant.stats.crv_types import RaisedCosineDistribution
    assert _test_args(RaisedCosineDistribution(1, 1))


def test_diofant__stats__crv_types__RayleighDistribution():
    from diofant.stats.crv_types import RayleighDistribution
    assert _test_args(RayleighDistribution(1))


def test_diofant__stats__crv_types__StudentTDistribution():
    from diofant.stats.crv_types import StudentTDistribution
    assert _test_args(StudentTDistribution(1))


def test_diofant__stats__crv_types__TriangularDistribution():
    from diofant.stats.crv_types import TriangularDistribution
    assert _test_args(TriangularDistribution(-1, 0, 1))


def test_diofant__stats__crv_types__UniformDistribution():
    from diofant.stats.crv_types import UniformDistribution
    assert _test_args(UniformDistribution(0, 1))


def test_diofant__stats__crv_types__UniformSumDistribution():
    from diofant.stats.crv_types import UniformSumDistribution
    assert _test_args(UniformSumDistribution(1))


def test_diofant__stats__crv_types__VonMisesDistribution():
    from diofant.stats.crv_types import VonMisesDistribution
    assert _test_args(VonMisesDistribution(1, 1))


def test_diofant__stats__crv_types__WeibullDistribution():
    from diofant.stats.crv_types import WeibullDistribution
    assert _test_args(WeibullDistribution(1, 1))


def test_diofant__stats__crv_types__WignerSemicircleDistribution():
    from diofant.stats.crv_types import WignerSemicircleDistribution
    assert _test_args(WignerSemicircleDistribution(1))


def test_diofant__stats__drv_types__PoissonDistribution():
    from diofant.stats.drv_types import PoissonDistribution
    assert _test_args(PoissonDistribution(1))


def test_diofant__stats__drv_types__GeometricDistribution():
    from diofant.stats.drv_types import GeometricDistribution
    assert _test_args(GeometricDistribution(.5))


def test_diofant__core__symbol__BaseSymbol():
    pass


def test_diofant__core__symbol__Dummy():
    from diofant.core.symbol import Dummy
    assert _test_args(Dummy('t'))


def test_diofant__core__symbol__Symbol():
    from diofant.core.symbol import Symbol
    assert _test_args(Symbol('t'))


def test_diofant__core__symbol__Wild():
    from diofant.core.symbol import Wild
    assert _test_args(Wild('x', exclude=[x]))


def test_diofant__functions__combinatorial__factorials__CombinatorialFunction():
    pass


def test_diofant__functions__combinatorial__factorials__FallingFactorial():
    from diofant.functions.combinatorial.factorials import FallingFactorial
    assert _test_args(FallingFactorial(2, x))


def test_diofant__functions__combinatorial__factorials__MultiFactorial():
    from diofant.functions.combinatorial.factorials import MultiFactorial
    assert _test_args(MultiFactorial(x))


def test_diofant__functions__combinatorial__factorials__RisingFactorial():
    from diofant.functions.combinatorial.factorials import RisingFactorial
    assert _test_args(RisingFactorial(2, x))


def test_diofant__functions__combinatorial__factorials__binomial():
    from diofant.functions.combinatorial.factorials import binomial
    assert _test_args(binomial(2, x))


def test_diofant__functions__combinatorial__factorials__subfactorial():
    from diofant.functions.combinatorial.factorials import subfactorial
    assert _test_args(subfactorial(1))


def test_diofant__functions__combinatorial__factorials__factorial():
    from diofant.functions.combinatorial.factorials import factorial
    assert _test_args(factorial(x))


def test_diofant__functions__combinatorial__factorials__factorial2():
    from diofant.functions.combinatorial.factorials import factorial2
    assert _test_args(factorial2(x))


def test_diofant__functions__combinatorial__numbers__bell():
    from diofant.functions.combinatorial.numbers import bell
    assert _test_args(bell(x, y))


def test_diofant__functions__combinatorial__numbers__bernoulli():
    from diofant.functions.combinatorial.numbers import bernoulli
    assert _test_args(bernoulli(x))


def test_diofant__functions__combinatorial__numbers__catalan():
    from diofant.functions.combinatorial.numbers import catalan
    assert _test_args(catalan(x))


def test_diofant__functions__combinatorial__numbers__genocchi():
    from diofant.functions.combinatorial.numbers import genocchi
    assert _test_args(genocchi(x))


def test_diofant__functions__combinatorial__numbers__euler():
    from diofant.functions.combinatorial.numbers import euler
    assert _test_args(euler(x))


def test_diofant__functions__combinatorial__numbers__fibonacci():
    from diofant.functions.combinatorial.numbers import fibonacci
    assert _test_args(fibonacci(x))


def test_diofant__functions__combinatorial__numbers__harmonic():
    from diofant.functions.combinatorial.numbers import harmonic
    assert _test_args(harmonic(x, 2))


def test_diofant__functions__combinatorial__numbers__lucas():
    from diofant.functions.combinatorial.numbers import lucas
    assert _test_args(lucas(x))


def test_diofant__functions__elementary__complexes__Abs():
    from diofant.functions.elementary.complexes import Abs
    assert _test_args(Abs(x))


def test_diofant__functions__elementary__complexes__adjoint():
    from diofant.functions.elementary.complexes import adjoint
    assert _test_args(adjoint(x))


def test_diofant__functions__elementary__complexes__arg():
    from diofant.functions.elementary.complexes import arg
    assert _test_args(arg(x))


def test_diofant__functions__elementary__complexes__conjugate():
    from diofant.functions.elementary.complexes import conjugate
    assert _test_args(conjugate(x))


def test_diofant__functions__elementary__complexes__im():
    from diofant.functions.elementary.complexes import im
    assert _test_args(im(x))


def test_diofant__functions__elementary__complexes__re():
    from diofant.functions.elementary.complexes import re
    assert _test_args(re(x))


def test_diofant__functions__elementary__complexes__sign():
    from diofant.functions.elementary.complexes import sign
    assert _test_args(sign(x))


def test_diofant__functions__elementary__complexes__polar_lift():
    from diofant.functions.elementary.complexes import polar_lift
    assert _test_args(polar_lift(x))


def test_diofant__functions__elementary__complexes__periodic_argument():
    from diofant.functions.elementary.complexes import periodic_argument
    assert _test_args(periodic_argument(x, y))


def test_diofant__functions__elementary__complexes__principal_branch():
    from diofant.functions.elementary.complexes import principal_branch
    assert _test_args(principal_branch(x, y))


def test_diofant__functions__elementary__complexes__transpose():
    from diofant.functions.elementary.complexes import transpose
    assert _test_args(transpose(x))


def test_diofant__functions__elementary__exponential__LambertW():
    from diofant.functions.elementary.exponential import LambertW
    assert _test_args(LambertW(2))


def test_diofant__functions__elementary__exponential__exp():
    from diofant.functions.elementary.exponential import exp
    assert _test_args(exp(2))


def test_diofant__functions__elementary__exponential__exp_polar():
    from diofant.functions.elementary.exponential import exp_polar
    assert _test_args(exp_polar(2))


def test_diofant__functions__elementary__exponential__log():
    from diofant.functions.elementary.exponential import log
    assert _test_args(log(2))


def test_diofant__functions__elementary__hyperbolic__HyperbolicFunction():
    pass


def test_diofant__functions__elementary__hyperbolic__ReciprocalHyperbolicFunction():
    pass


def test_diofant__functions__elementary__hyperbolic__acosh():
    from diofant.functions.elementary.hyperbolic import acosh
    assert _test_args(acosh(2))


def test_diofant__functions__elementary__hyperbolic__acoth():
    from diofant.functions.elementary.hyperbolic import acoth
    assert _test_args(acoth(2))


def test_diofant__functions__elementary__hyperbolic__asinh():
    from diofant.functions.elementary.hyperbolic import asinh
    assert _test_args(asinh(2))


def test_diofant__functions__elementary__hyperbolic__atanh():
    from diofant.functions.elementary.hyperbolic import atanh
    assert _test_args(atanh(2))


def test_diofant__functions__elementary__hyperbolic__cosh():
    from diofant.functions.elementary.hyperbolic import cosh
    assert _test_args(cosh(2))


def test_diofant__functions__elementary__hyperbolic__coth():
    from diofant.functions.elementary.hyperbolic import coth
    assert _test_args(coth(2))


def test_diofant__functions__elementary__hyperbolic__csch():
    from diofant.functions.elementary.hyperbolic import csch
    assert _test_args(csch(2))


def test_diofant__functions__elementary__hyperbolic__sech():
    from diofant.functions.elementary.hyperbolic import sech
    assert _test_args(sech(2))


def test_diofant__functions__elementary__hyperbolic__sinh():
    from diofant.functions.elementary.hyperbolic import sinh
    assert _test_args(sinh(2))


def test_diofant__functions__elementary__hyperbolic__tanh():
    from diofant.functions.elementary.hyperbolic import tanh
    assert _test_args(tanh(2))


@pytest.mark.xfail
def test_diofant__functions__elementary__integers__RoundFunction():
    from diofant.functions.elementary.integers import RoundFunction
    assert _test_args(RoundFunction())


def test_diofant__functions__elementary__integers__ceiling():
    from diofant.functions.elementary.integers import ceiling
    assert _test_args(ceiling(x))


def test_diofant__functions__elementary__integers__floor():
    from diofant.functions.elementary.integers import floor
    assert _test_args(floor(x))


def test_diofant__functions__elementary__miscellaneous__IdentityFunction():
    from diofant.functions.elementary.miscellaneous import IdentityFunction
    assert _test_args(IdentityFunction())


def test_diofant__functions__elementary__miscellaneous__Max():
    from diofant.functions.elementary.miscellaneous import Max
    assert _test_args(Max(x, 2))


def test_diofant__functions__elementary__miscellaneous__Min():
    from diofant.functions.elementary.miscellaneous import Min
    assert _test_args(Min(x, 2))


def test_diofant__functions__elementary__miscellaneous__MinMaxBase():
    pass


def test_diofant__functions__elementary__piecewise__ExprCondPair():
    from diofant.functions.elementary.piecewise import ExprCondPair
    assert _test_args(ExprCondPair(1, True))


def test_diofant__functions__elementary__piecewise__Piecewise():
    from diofant.functions.elementary.piecewise import Piecewise
    assert _test_args(Piecewise((1, x >= 0), (0, True)))


def test_diofant__functions__elementary__trigonometric__TrigonometricFunction():
    pass


def test_diofant__functions__elementary__trigonometric__ReciprocalTrigonometricFunction():
    pass


def test_diofant__functions__elementary__trigonometric__InverseTrigonometricFunction():
    pass


def test_diofant__functions__elementary__trigonometric__acos():
    from diofant.functions.elementary.trigonometric import acos
    assert _test_args(acos(2))


def test_diofant__functions__elementary__trigonometric__acot():
    from diofant.functions.elementary.trigonometric import acot
    assert _test_args(acot(2))


def test_diofant__functions__elementary__trigonometric__asin():
    from diofant.functions.elementary.trigonometric import asin
    assert _test_args(asin(2))


def test_diofant__functions__elementary__trigonometric__asec():
    from diofant.functions.elementary.trigonometric import asec
    assert _test_args(asec(2))


def test_diofant__functions__elementary__trigonometric__acsc():
    from diofant.functions.elementary.trigonometric import acsc
    assert _test_args(acsc(2))


def test_diofant__functions__elementary__trigonometric__atan():
    from diofant.functions.elementary.trigonometric import atan
    assert _test_args(atan(2))


def test_diofant__functions__elementary__trigonometric__atan2():
    from diofant.functions.elementary.trigonometric import atan2
    assert _test_args(atan2(2, 3))


def test_diofant__functions__elementary__trigonometric__cos():
    from diofant.functions.elementary.trigonometric import cos
    assert _test_args(cos(2))


def test_diofant__functions__elementary__trigonometric__csc():
    from diofant.functions.elementary.trigonometric import csc
    assert _test_args(csc(2))


def test_diofant__functions__elementary__trigonometric__cot():
    from diofant.functions.elementary.trigonometric import cot
    assert _test_args(cot(2))


def test_diofant__functions__elementary__trigonometric__sin():
    assert _test_args(sin(2))


def test_diofant__functions__elementary__trigonometric__sec():
    from diofant.functions.elementary.trigonometric import sec
    assert _test_args(sec(2))


def test_diofant__functions__elementary__trigonometric__tan():
    from diofant.functions.elementary.trigonometric import tan
    assert _test_args(tan(2))


def test_diofant__functions__special__bessel__BesselBase():
    pass


def test_diofant__functions__special__bessel__SphericalBesselBase():
    pass


def test_diofant__functions__special__bessel__besseli():
    from diofant.functions.special.bessel import besseli
    assert _test_args(besseli(x, 1))


def test_diofant__functions__special__bessel__besselj():
    from diofant.functions.special.bessel import besselj
    assert _test_args(besselj(x, 1))


def test_diofant__functions__special__bessel__besselk():
    from diofant.functions.special.bessel import besselk
    assert _test_args(besselk(x, 1))


def test_diofant__functions__special__bessel__bessely():
    from diofant.functions.special.bessel import bessely
    assert _test_args(bessely(x, 1))


def test_diofant__functions__special__bessel__hankel1():
    from diofant.functions.special.bessel import hankel1
    assert _test_args(hankel1(x, 1))


def test_diofant__functions__special__bessel__hankel2():
    from diofant.functions.special.bessel import hankel2
    assert _test_args(hankel2(x, 1))


def test_diofant__functions__special__bessel__jn():
    from diofant.functions.special.bessel import jn
    assert _test_args(jn(0, x))


def test_diofant__functions__special__bessel__yn():
    from diofant.functions.special.bessel import yn
    assert _test_args(yn(0, x))


def test_diofant__functions__special__bessel__AiryBase():
    pass


def test_diofant__functions__special__bessel__airyai():
    from diofant.functions.special.bessel import airyai
    assert _test_args(airyai(2))


def test_diofant__functions__special__bessel__airybi():
    from diofant.functions.special.bessel import airybi
    assert _test_args(airybi(2))


def test_diofant__functions__special__bessel__airyaiprime():
    from diofant.functions.special.bessel import airyaiprime
    assert _test_args(airyaiprime(2))


def test_diofant__functions__special__bessel__airybiprime():
    from diofant.functions.special.bessel import airybiprime
    assert _test_args(airybiprime(2))


def test_diofant__functions__special__elliptic_integrals__elliptic_k():
    from diofant.functions.special.elliptic_integrals import elliptic_k as K
    assert _test_args(K(x))


def test_diofant__functions__special__elliptic_integrals__elliptic_f():
    from diofant.functions.special.elliptic_integrals import elliptic_f as F
    assert _test_args(F(x, y))


def test_diofant__functions__special__elliptic_integrals__elliptic_e():
    from diofant.functions.special.elliptic_integrals import elliptic_e as E
    assert _test_args(E(x))
    assert _test_args(E(x, y))


def test_diofant__functions__special__elliptic_integrals__elliptic_pi():
    from diofant.functions.special.elliptic_integrals import elliptic_pi as P
    assert _test_args(P(x, y))
    assert _test_args(P(x, y, z))


def test_diofant__functions__special__delta_functions__DiracDelta():
    from diofant.functions.special.delta_functions import DiracDelta
    assert _test_args(DiracDelta(x, 1))


def test_diofant__functions__special__delta_functions__Heaviside():
    from diofant.functions.special.delta_functions import Heaviside
    assert _test_args(Heaviside(x))


def test_diofant__functions__special__error_functions__erf():
    from diofant.functions.special.error_functions import erf
    assert _test_args(erf(2))


def test_diofant__functions__special__error_functions__erfc():
    from diofant.functions.special.error_functions import erfc
    assert _test_args(erfc(2))


def test_diofant__functions__special__error_functions__erfi():
    from diofant.functions.special.error_functions import erfi
    assert _test_args(erfi(2))


def test_diofant__functions__special__error_functions__erf2():
    from diofant.functions.special.error_functions import erf2
    assert _test_args(erf2(2, 3))


def test_diofant__functions__special__error_functions__erfinv():
    from diofant.functions.special.error_functions import erfinv
    assert _test_args(erfinv(2))


def test_diofant__functions__special__error_functions__erfcinv():
    from diofant.functions.special.error_functions import erfcinv
    assert _test_args(erfcinv(2))


def test_diofant__functions__special__error_functions__erf2inv():
    from diofant.functions.special.error_functions import erf2inv
    assert _test_args(erf2inv(2, 3))


def test_diofant__functions__special__error_functions__FresnelIntegral():
    pass


def test_diofant__functions__special__error_functions__fresnels():
    from diofant.functions.special.error_functions import fresnels
    assert _test_args(fresnels(2))


def test_diofant__functions__special__error_functions__fresnelc():
    from diofant.functions.special.error_functions import fresnelc
    assert _test_args(fresnelc(2))


def test_diofant__functions__special__error_functions__erfs():
    from diofant.functions.special.error_functions import _erfs
    assert _test_args(_erfs(2))


def test_diofant__functions__special__error_functions__Ei():
    from diofant.functions.special.error_functions import Ei
    assert _test_args(Ei(2))


def test_diofant__functions__special__error_functions__li():
    from diofant.functions.special.error_functions import li
    assert _test_args(li(2))


def test_diofant__functions__special__error_functions__Li():
    from diofant.functions.special.error_functions import Li
    assert _test_args(Li(2))


def test_diofant__functions__special__error_functions__TrigonometricIntegral():
    pass


def test_diofant__functions__special__error_functions__Si():
    from diofant.functions.special.error_functions import Si
    assert _test_args(Si(2))


def test_diofant__functions__special__error_functions__Ci():
    from diofant.functions.special.error_functions import Ci
    assert _test_args(Ci(2))


def test_diofant__functions__special__error_functions__Shi():
    from diofant.functions.special.error_functions import Shi
    assert _test_args(Shi(2))


def test_diofant__functions__special__error_functions__Chi():
    from diofant.functions.special.error_functions import Chi
    assert _test_args(Chi(2))


def test_diofant__functions__special__error_functions__expint():
    from diofant.functions.special.error_functions import expint
    assert _test_args(expint(y, x))


def test_diofant__functions__special__gamma_functions__gamma():
    from diofant.functions.special.gamma_functions import gamma
    assert _test_args(gamma(x))


def test_diofant__functions__special__gamma_functions__loggamma():
    from diofant.functions.special.gamma_functions import loggamma
    assert _test_args(loggamma(2))


def test_diofant__functions__special__gamma_functions__lowergamma():
    from diofant.functions.special.gamma_functions import lowergamma
    assert _test_args(lowergamma(x, 2))


def test_diofant__functions__special__gamma_functions__polygamma():
    from diofant.functions.special.gamma_functions import polygamma
    assert _test_args(polygamma(x, 2))


def test_diofant__functions__special__gamma_functions__uppergamma():
    from diofant.functions.special.gamma_functions import uppergamma
    assert _test_args(uppergamma(x, 2))


def test_diofant__functions__special__beta_functions__beta():
    from diofant.functions.special.beta_functions import beta
    assert _test_args(beta(x, x))


def test_diofant__functions__special__hyper__TupleParametersBase():
    pass


def test_diofant__functions__special__hyper__TupleArg():
    pass


def test_diofant__functions__special__hyper__hyper():
    from diofant.functions.special.hyper import hyper
    assert _test_args(hyper([1, 2, 3], [4, 5], x))


def test_diofant__functions__special__hyper__meijerg():
    from diofant.functions.special.hyper import meijerg
    assert _test_args(meijerg([1, 2, 3], [4, 5], [6], [], x))


def test_diofant__functions__special__hyper__HyperRep():
    pass


def test_diofant__functions__special__hyper__HyperRep_power1():
    from diofant.functions.special.hyper import HyperRep_power1
    assert _test_args(HyperRep_power1(x, y))


def test_diofant__functions__special__hyper__HyperRep_power2():
    from diofant.functions.special.hyper import HyperRep_power2
    assert _test_args(HyperRep_power2(x, y))


def test_diofant__functions__special__hyper__HyperRep_log1():
    from diofant.functions.special.hyper import HyperRep_log1
    assert _test_args(HyperRep_log1(x))


def test_diofant__functions__special__hyper__HyperRep_atanh():
    from diofant.functions.special.hyper import HyperRep_atanh
    assert _test_args(HyperRep_atanh(x))


def test_diofant__functions__special__hyper__HyperRep_asin1():
    from diofant.functions.special.hyper import HyperRep_asin1
    assert _test_args(HyperRep_asin1(x))


def test_diofant__functions__special__hyper__HyperRep_asin2():
    from diofant.functions.special.hyper import HyperRep_asin2
    assert _test_args(HyperRep_asin2(x))


def test_diofant__functions__special__hyper__HyperRep_sqrts1():
    from diofant.functions.special.hyper import HyperRep_sqrts1
    assert _test_args(HyperRep_sqrts1(x, y))


def test_diofant__functions__special__hyper__HyperRep_sqrts2():
    from diofant.functions.special.hyper import HyperRep_sqrts2
    assert _test_args(HyperRep_sqrts2(x, y))


def test_diofant__functions__special__hyper__HyperRep_log2():
    from diofant.functions.special.hyper import HyperRep_log2
    assert _test_args(HyperRep_log2(x))


def test_diofant__functions__special__hyper__HyperRep_cosasin():
    from diofant.functions.special.hyper import HyperRep_cosasin
    assert _test_args(HyperRep_cosasin(x, y))


def test_diofant__functions__special__hyper__HyperRep_sinasin():
    from diofant.functions.special.hyper import HyperRep_sinasin
    assert _test_args(HyperRep_sinasin(x, y))


def test_diofant__functions__special__polynomials__OrthogonalPolynomial():
    pass


def test_diofant__functions__special__polynomials__jacobi():
    from diofant.functions.special.polynomials import jacobi
    assert _test_args(jacobi(x, 2, 2, 2))


def test_diofant__functions__special__polynomials__gegenbauer():
    from diofant.functions.special.polynomials import gegenbauer
    assert _test_args(gegenbauer(x, 2, 2))


def test_diofant__functions__special__polynomials__chebyshevt():
    from diofant.functions.special.polynomials import chebyshevt
    assert _test_args(chebyshevt(x, 2))


def test_diofant__functions__special__polynomials__chebyshevt_root():
    from diofant.functions.special.polynomials import chebyshevt_root
    assert _test_args(chebyshevt_root(3, 2))


def test_diofant__functions__special__polynomials__chebyshevu():
    from diofant.functions.special.polynomials import chebyshevu
    assert _test_args(chebyshevu(x, 2))


def test_diofant__functions__special__polynomials__chebyshevu_root():
    from diofant.functions.special.polynomials import chebyshevu_root
    assert _test_args(chebyshevu_root(3, 2))


def test_diofant__functions__special__polynomials__hermite():
    from diofant.functions.special.polynomials import hermite
    assert _test_args(hermite(x, 2))


def test_diofant__functions__special__polynomials__legendre():
    from diofant.functions.special.polynomials import legendre
    assert _test_args(legendre(x, 2))


def test_diofant__functions__special__polynomials__assoc_legendre():
    from diofant.functions.special.polynomials import assoc_legendre
    assert _test_args(assoc_legendre(x, 0, y))


def test_diofant__functions__special__polynomials__laguerre():
    from diofant.functions.special.polynomials import laguerre
    assert _test_args(laguerre(x, 2))


def test_diofant__functions__special__polynomials__assoc_laguerre():
    from diofant.functions.special.polynomials import assoc_laguerre
    assert _test_args(assoc_laguerre(x, 0, y))


def test_diofant__functions__special__spherical_harmonics__Ynm():
    from diofant.functions.special.spherical_harmonics import Ynm
    assert _test_args(Ynm(1, 1, x, y))


def test_diofant__functions__special__spherical_harmonics__Znm():
    from diofant.functions.special.spherical_harmonics import Znm
    assert _test_args(Znm(1, 1, x, y))


def test_diofant__functions__special__tensor_functions__LeviCivita():
    from diofant.functions.special.tensor_functions import LeviCivita
    assert _test_args(LeviCivita(x, y, 2))


def test_diofant__functions__special__tensor_functions__KroneckerDelta():
    from diofant.functions.special.tensor_functions import KroneckerDelta
    assert _test_args(KroneckerDelta(x, y))


def test_diofant__functions__special__zeta_functions__dirichlet_eta():
    from diofant.functions.special.zeta_functions import dirichlet_eta
    assert _test_args(dirichlet_eta(x))


def test_diofant__functions__special__zeta_functions__zeta():
    from diofant.functions.special.zeta_functions import zeta
    assert _test_args(zeta(101))


def test_diofant__functions__special__zeta_functions__lerchphi():
    from diofant.functions.special.zeta_functions import lerchphi
    assert _test_args(lerchphi(x, y, z))


def test_diofant__functions__special__zeta_functions__polylog():
    from diofant.functions.special.zeta_functions import polylog
    assert _test_args(polylog(x, y))


def test_diofant__integrals__integrals__Integral():
    from diofant.integrals.integrals import Integral
    assert _test_args(Integral(2, (x, 0, 1)))


def test_diofant__integrals__risch__NonElementaryIntegral():
    from diofant.integrals.risch import NonElementaryIntegral
    assert _test_args(NonElementaryIntegral(exp(-x**2), x))


def test_diofant__integrals__transforms__IntegralTransform():
    pass


def test_diofant__integrals__transforms__MellinTransform():
    from diofant.integrals.transforms import MellinTransform
    assert _test_args(MellinTransform(2, x, y))


def test_diofant__integrals__transforms__InverseMellinTransform():
    from diofant.integrals.transforms import InverseMellinTransform
    assert _test_args(InverseMellinTransform(2, x, y, 0, 1))


def test_diofant__integrals__transforms__LaplaceTransform():
    from diofant.integrals.transforms import LaplaceTransform
    assert _test_args(LaplaceTransform(2, x, y))


def test_diofant__integrals__transforms__InverseLaplaceTransform():
    from diofant.integrals.transforms import InverseLaplaceTransform
    assert _test_args(InverseLaplaceTransform(2, x, y, 0))


def test_diofant__integrals__transforms__FourierTypeTransform():
    pass


def test_diofant__integrals__transforms__InverseFourierTransform():
    from diofant.integrals.transforms import InverseFourierTransform
    assert _test_args(InverseFourierTransform(2, x, y))


def test_diofant__integrals__transforms__FourierTransform():
    from diofant.integrals.transforms import FourierTransform
    assert _test_args(FourierTransform(2, x, y))


def test_diofant__integrals__transforms__SineCosineTypeTransform():
    pass


def test_diofant__integrals__transforms__InverseSineTransform():
    from diofant.integrals.transforms import InverseSineTransform
    assert _test_args(InverseSineTransform(2, x, y))


def test_diofant__integrals__transforms__SineTransform():
    from diofant.integrals.transforms import SineTransform
    assert _test_args(SineTransform(2, x, y))


def test_diofant__integrals__transforms__InverseCosineTransform():
    from diofant.integrals.transforms import InverseCosineTransform
    assert _test_args(InverseCosineTransform(2, x, y))


def test_diofant__integrals__transforms__CosineTransform():
    from diofant.integrals.transforms import CosineTransform
    assert _test_args(CosineTransform(2, x, y))


def test_diofant__integrals__transforms__HankelTypeTransform():
    pass


def test_diofant__integrals__transforms__InverseHankelTransform():
    from diofant.integrals.transforms import InverseHankelTransform
    assert _test_args(InverseHankelTransform(2, x, y, 0))


def test_diofant__integrals__transforms__HankelTransform():
    from diofant.integrals.transforms import HankelTransform
    assert _test_args(HankelTransform(2, x, y, 0))


def test_diofant__logic__boolalg__And():
    from diofant.logic.boolalg import And
    assert _test_args(And(x, y, 2))


def test_diofant__logic__boolalg__Boolean():
    pass


def test_diofant__logic__boolalg__BooleanFunction():
    from diofant.logic.boolalg import BooleanFunction
    assert _test_args(BooleanFunction(1, 2, 3))


def test_diofant__logic__boolalg__BooleanAtom():
    pass


def test_diofant__logic__boolalg__BooleanTrue():
    from diofant.logic.boolalg import true
    assert _test_args(true)


def test_diofant__logic__boolalg__BooleanFalse():
    from diofant.logic.boolalg import false
    assert _test_args(false)


def test_diofant__logic__boolalg__Equivalent():
    from diofant.logic.boolalg import Equivalent
    assert _test_args(Equivalent(x, 2))


def test_diofant__logic__boolalg__ITE():
    from diofant.logic.boolalg import ITE
    assert _test_args(ITE(x, y, 2))


def test_diofant__logic__boolalg__Implies():
    from diofant.logic.boolalg import Implies
    assert _test_args(Implies(x, y))


def test_diofant__logic__boolalg__Nand():
    from diofant.logic.boolalg import Nand
    assert _test_args(Nand(x, y, 2))


def test_diofant__logic__boolalg__Nor():
    from diofant.logic.boolalg import Nor
    assert _test_args(Nor(x, y))


def test_diofant__logic__boolalg__Not():
    from diofant.logic.boolalg import Not
    assert _test_args(Not(x))


def test_diofant__logic__boolalg__Or():
    from diofant.logic.boolalg import Or
    assert _test_args(Or(x, y))


def test_diofant__logic__boolalg__Xor():
    from diofant.logic.boolalg import Xor
    assert _test_args(Xor(x, y, 2))


def test_diofant__matrices__matrices__DeferredVector():
    from diofant.matrices.matrices import DeferredVector
    assert _test_args(DeferredVector("X"))


def test_diofant__matrices__expressions__matexpr__MatrixBase():
    pass


def test_diofant__matrices__immutable__ImmutableMatrix():
    from diofant.matrices.immutable import ImmutableMatrix
    m = ImmutableMatrix([[1, 2], [3, 4]])
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))
    m = ImmutableMatrix(1, 1, [1])
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))
    m = ImmutableMatrix(2, 2, lambda i, j: 1)
    assert m[0, 0] is S.One
    m = ImmutableMatrix(2, 2, lambda i, j: 1/(1 + i) + 1/(1 + j))
    assert m[1, 1] is S.One  # true div. will give 1.0 if i,j not sympified
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))


def test_diofant__matrices__immutable__ImmutableSparseMatrix():
    from diofant.matrices.immutable import ImmutableSparseMatrix
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
    assert m[0, 0] is S.One
    m = ImmutableSparseMatrix(2, 2, lambda i, j: 1/(1 + i) + 1/(1 + j))
    assert m[1, 1] is S.One  # true div. will give 1.0 if i,j not sympified
    assert _test_args(m)
    assert _test_args(Basic(*list(m)))


def test_diofant__matrices__expressions__slice__MatrixSlice():
    from diofant.matrices.expressions.slice import MatrixSlice
    from diofant.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', 4, 4)
    assert _test_args(MatrixSlice(X, (0, 2), (0, 2)))


def test_diofant__matrices__expressions__blockmatrix__BlockDiagMatrix():
    from diofant.matrices.expressions.blockmatrix import BlockDiagMatrix
    from diofant.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, x)
    Y = MatrixSymbol('Y', y, y)
    assert _test_args(BlockDiagMatrix(X, Y))


def test_diofant__matrices__expressions__blockmatrix__BlockMatrix():
    from diofant.matrices.expressions.blockmatrix import BlockMatrix
    from diofant.matrices.expressions import MatrixSymbol, ZeroMatrix
    X = MatrixSymbol('X', x, x)
    Y = MatrixSymbol('Y', y, y)
    Z = MatrixSymbol('Z', x, y)
    O = ZeroMatrix(y, x)
    assert _test_args(BlockMatrix([[X, Z], [O, Y]]))


def test_diofant__matrices__expressions__inverse__Inverse():
    from diofant.matrices.expressions.inverse import Inverse
    from diofant.matrices.expressions import MatrixSymbol
    assert _test_args(Inverse(MatrixSymbol('A', 3, 3)))


def test_diofant__matrices__expressions__matadd__MatAdd():
    from diofant.matrices.expressions.matadd import MatAdd
    from diofant.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', x, y)
    assert _test_args(MatAdd(X, Y))


def test_diofant__matrices__expressions__matexpr__Identity():
    from diofant.matrices.expressions.matexpr import Identity
    assert _test_args(Identity(3))


def test_diofant__matrices__expressions__matexpr__MatrixExpr():
    pass


def test_diofant__matrices__expressions__matexpr__MatrixElement():
    from diofant.matrices.expressions.matexpr import MatrixSymbol, MatrixElement
    assert _test_args(MatrixElement(MatrixSymbol('A', 3, 5), Integer(2), Integer(3)))


@pytest.mark.xfail
def test_diofant__matrices__expressions__matexpr__MatrixSymbol():
    from diofant.matrices.expressions.matexpr import MatrixSymbol
    assert _test_args(MatrixSymbol('A', 3, 5))


def test_diofant__matrices__expressions__matexpr__ZeroMatrix():
    from diofant.matrices.expressions.matexpr import ZeroMatrix
    assert _test_args(ZeroMatrix(3, 5))


def test_diofant__matrices__expressions__matmul__MatMul():
    from diofant.matrices.expressions.matmul import MatMul
    from diofant.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', y, x)
    assert _test_args(MatMul(X, Y))


def test_diofant__matrices__expressions__diagonal__DiagonalMatrix():
    from diofant.matrices.expressions.diagonal import DiagonalMatrix
    from diofant.matrices.expressions import MatrixSymbol
    x = MatrixSymbol('x', 10, 1)
    assert _test_args(DiagonalMatrix(x))


def test_diofant__matrices__expressions__diagonal__DiagonalOf():
    from diofant.matrices.expressions.diagonal import DiagonalOf
    from diofant.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('x', 10, 10)
    assert _test_args(DiagonalOf(X))


def test_diofant__matrices__expressions__hadamard__HadamardProduct():
    from diofant.matrices.expressions.hadamard import HadamardProduct
    from diofant.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', x, y)
    assert _test_args(HadamardProduct(X, Y))


def test_diofant__matrices__expressions__matpow__MatPow():
    from diofant.matrices.expressions.matpow import MatPow
    from diofant.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, x)
    assert _test_args(MatPow(X, 2))


def test_diofant__matrices__expressions__transpose__Transpose():
    from diofant.matrices.expressions.transpose import Transpose
    from diofant.matrices.expressions import MatrixSymbol
    assert _test_args(Transpose(MatrixSymbol('A', 3, 5)))


def test_diofant__matrices__expressions__adjoint__Adjoint():
    from diofant.matrices.expressions.adjoint import Adjoint
    from diofant.matrices.expressions import MatrixSymbol
    assert _test_args(Adjoint(MatrixSymbol('A', 3, 5)))


def test_diofant__matrices__expressions__trace__Trace():
    from diofant.matrices.expressions.trace import Trace
    from diofant.matrices.expressions import MatrixSymbol
    assert _test_args(Trace(MatrixSymbol('A', 3, 3)))


def test_diofant__matrices__expressions__determinant__Determinant():
    from diofant.matrices.expressions.determinant import Determinant
    from diofant.matrices.expressions import MatrixSymbol
    assert _test_args(Determinant(MatrixSymbol('A', 3, 3)))


def test_diofant__matrices__expressions__funcmatrix__FunctionMatrix():
    from diofant.matrices.expressions.funcmatrix import FunctionMatrix
    from diofant import Lambda, symbols
    i, j = symbols('i,j')
    assert _test_args(FunctionMatrix(3, 3, Lambda((i, j), i - j) ))


def test_diofant__matrices__expressions__fourier__DFT():
    from diofant.matrices.expressions.fourier import DFT
    assert _test_args(DFT(Integer(2)))


def test_diofant__matrices__expressions__fourier__IDFT():
    from diofant.matrices.expressions.fourier import IDFT
    assert _test_args(IDFT(Integer(2)))

from diofant.matrices.expressions import MatrixSymbol
X = MatrixSymbol('X', 10, 10)


def test_diofant__matrices__expressions__factorizations__LofLU():
    from diofant.matrices.expressions.factorizations import LofLU
    assert _test_args(LofLU(X))


def test_diofant__matrices__expressions__factorizations__UofLU():
    from diofant.matrices.expressions.factorizations import UofLU
    assert _test_args(UofLU(X))


def test_diofant__matrices__expressions__factorizations__QofQR():
    from diofant.matrices.expressions.factorizations import QofQR
    assert _test_args(QofQR(X))


def test_diofant__matrices__expressions__factorizations__RofQR():
    from diofant.matrices.expressions.factorizations import RofQR
    assert _test_args(RofQR(X))


def test_diofant__matrices__expressions__factorizations__LofCholesky():
    from diofant.matrices.expressions.factorizations import LofCholesky
    assert _test_args(LofCholesky(X))


def test_diofant__matrices__expressions__factorizations__UofCholesky():
    from diofant.matrices.expressions.factorizations import UofCholesky
    assert _test_args(UofCholesky(X))


def test_diofant__matrices__expressions__factorizations__EigenVectors():
    from diofant.matrices.expressions.factorizations import EigenVectors
    assert _test_args(EigenVectors(X))


def test_diofant__matrices__expressions__factorizations__EigenValues():
    from diofant.matrices.expressions.factorizations import EigenValues
    assert _test_args(EigenValues(X))


def test_diofant__matrices__expressions__factorizations__UofSVD():
    from diofant.matrices.expressions.factorizations import UofSVD
    assert _test_args(UofSVD(X))


def test_diofant__matrices__expressions__factorizations__VofSVD():
    from diofant.matrices.expressions.factorizations import VofSVD
    assert _test_args(VofSVD(X))


def test_diofant__matrices__expressions__factorizations__SofSVD():
    from diofant.matrices.expressions.factorizations import SofSVD
    assert _test_args(SofSVD(X))


def test_diofant__matrices__expressions__factorizations__Factorization():
    pass


def test_diofant__core__numbers__AlgebraicNumber():
    from diofant.core.numbers import AlgebraicNumber
    assert _test_args(AlgebraicNumber(sqrt(2), [1, 2, 3]))


def test_diofant__polys__polytools__GroebnerBasis():
    from diofant.polys.polytools import GroebnerBasis
    assert _test_args(GroebnerBasis([x, y, z], x, y, z))


def test_diofant__polys__polytools__Poly():
    from diofant.polys.polytools import Poly
    assert _test_args(Poly(2, x, y))


def test_diofant__polys__polytools__PurePoly():
    from diofant.polys.polytools import PurePoly
    assert _test_args(PurePoly(2, x, y))


def test_diofant__polys__rootoftools__RootOf():
    from diofant.polys.rootoftools import RootOf
    assert _test_args(RootOf(x**3 + x + 1, 0))


def test_diofant__polys__rootoftools__RootSum():
    from diofant.polys.rootoftools import RootSum
    assert _test_args(RootSum(x**3 + x + 1, sin))


def test_diofant__series__limits__Limit():
    from diofant.series.limits import Limit
    assert _test_args(Limit(x, x, 0, dir='-'))


def test_diofant__series__order__Order():
    from diofant.series.order import Order
    assert _test_args(Order(1, x, y))


def test_diofant__simplify__hyperexpand__Hyper_Function():
    from diofant.simplify.hyperexpand import Hyper_Function
    assert _test_args(Hyper_Function([2], [1]))


def test_diofant__simplify__hyperexpand__G_Function():
    from diofant.simplify.hyperexpand import G_Function
    assert _test_args(G_Function([2], [1], [], []))


def test_diofant__tensor__array__dense_ndim_array__ImmutableDenseNDimArray():
    from diofant.tensor.array.dense_ndim_array import ImmutableDenseNDimArray
    densarr = ImmutableDenseNDimArray(range(10, 34), (2, 3, 4))
    assert _test_args(densarr)


def test_diofant__tensor__array__sparse_ndim_array__ImmutableSparseNDimArray():
    from diofant.tensor.array.sparse_ndim_array import ImmutableSparseNDimArray
    sparr = ImmutableSparseNDimArray(range(10, 34), (2, 3, 4))
    assert _test_args(sparr)


def test_diofant__tensor__indexed__Idx():
    from diofant.tensor.indexed import Idx
    assert _test_args(Idx('test'))
    assert _test_args(Idx(1, (0, 10)))


def test_diofant__tensor__indexed__Indexed():
    from diofant.tensor.indexed import Indexed, Idx
    assert _test_args(Indexed('A', Idx('i'), Idx('j')))


def test_diofant__tensor__indexed__IndexedBase():
    from diofant.tensor.indexed import IndexedBase
    assert _test_args(IndexedBase('A', shape=(x, y)))
    assert _test_args(IndexedBase('A', 1))
    assert _test_args(IndexedBase('A')[0, 1])


def test_diofant__tensor__tensor__TensorIndexType():
    from diofant.tensor.tensor import TensorIndexType
    assert _test_args(TensorIndexType('Lorentz', metric=False))


def test_diofant__tensor__tensor__TensorSymmetry():
    from diofant.tensor.tensor import TensorSymmetry, get_symmetric_group_sgs
    assert _test_args(TensorSymmetry(get_symmetric_group_sgs(2)))


def test_diofant__tensor__tensor__TensorType():
    from diofant.tensor.tensor import TensorIndexType, TensorSymmetry, get_symmetric_group_sgs, TensorType
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    assert _test_args(TensorType([Lorentz], sym))


def test_diofant__tensor__tensor__TensorHead():
    from diofant.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, TensorHead
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    assert _test_args(TensorHead('p', S1, 0))


def test_diofant__tensor__tensor__TensorIndex():
    from diofant.tensor.tensor import TensorIndexType, TensorIndex
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    assert _test_args(TensorIndex('i', Lorentz))


def test_diofant__tensor__tensor__TensExpr():
    pass


def test_diofant__tensor__tensor__TensAdd():
    from diofant.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, tensor_indices, TensAdd
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p, q = S1('p,q')
    t1 = p(a)
    t2 = q(a)
    assert _test_args(TensAdd(t1, t2))


def test_diofant__tensor__tensor__Tensor():
    from diofant.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, tensor_indices
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p = S1('p')
    assert _test_args(p(a))


def test_diofant__tensor__tensor__TensMul():
    from diofant.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, tensor_indices
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
    from diofant.geometry.curve import Curve
    assert _test_args(Curve((x, 1), (x, 0, 1)))


def test_diofant__geometry__point__Point():
    from diofant.geometry.point import Point
    assert _test_args(Point(0, 1))


def test_diofant__geometry__point__Point2D():
    from diofant.geometry.point import Point2D
    assert _test_args(Point2D(0, 1))


def test_diofant__geometry__point__Point3D():
    from diofant.geometry.point import Point3D
    assert _test_args(Point3D(0, 1, 2))


def test_diofant__geometry__ellipse__Ellipse():
    from diofant.geometry.ellipse import Ellipse
    assert _test_args(Ellipse((0, 1), 2, 3))


def test_diofant__geometry__ellipse__Circle():
    from diofant.geometry.ellipse import Circle
    assert _test_args(Circle((0, 1), 2))


def test_diofant__geometry__line__LinearEntity():
    pass


def test_diofant__geometry__line__Line():
    from diofant.geometry.line import Line
    assert _test_args(Line((0, 1), (2, 3)))


def test_diofant__geometry__line__Ray():
    from diofant.geometry.line import Ray
    assert _test_args(Ray((0, 1), (2, 3)))


def test_diofant__geometry__line__Segment():
    from diofant.geometry.line import Segment
    assert _test_args(Segment((0, 1), (2, 3)))


def test_diofant__geometry__line3d__LinearEntity3D():
    pass


def test_diofant__geometry__line3d__Line3D():
    from diofant.geometry.line3d import Line3D
    assert _test_args(Line3D((0, 1, 1), (2, 3, 4)))


def test_diofant__geometry__line3d__Segment3D():
    from diofant.geometry.line3d import Segment3D
    assert _test_args(Segment3D((0, 1, 1), (2, 3, 4)))


def test_diofant__geometry__line3d__Ray3D():
    from diofant.geometry.line3d import Ray3D
    assert _test_args(Ray3D((0, 1, 1), (2, 3, 4)))


def test_diofant__geometry__plane__Plane():
    from diofant.geometry.plane import Plane
    assert _test_args(Plane((1, 1, 1), (-3, 4, -2), (1, 2, 3)))


def test_diofant__geometry__polygon__Polygon():
    from diofant.geometry.polygon import Polygon
    assert _test_args(Polygon((0, 1), (2, 3), (4, 5), (6, 7)))


def test_diofant__geometry__polygon__RegularPolygon():
    from diofant.geometry.polygon import RegularPolygon
    assert _test_args(RegularPolygon((0, 1), 2, 3, 4))


def test_diofant__geometry__polygon__Triangle():
    from diofant.geometry.polygon import Triangle
    assert _test_args(Triangle((0, 1), (2, 3), (4, 5)))


def test_diofant__geometry__entity__GeometryEntity():
    from diofant.geometry.entity import GeometryEntity
    from diofant.geometry.point import Point
    assert _test_args(GeometryEntity(Point(1, 0), 1, [1, 2]))


def test_diofant__geometry__entity__GeometrySet():
    pass


def test_diofant__diffgeom__diffgeom__Manifold():
    from diofant.diffgeom import Manifold
    assert _test_args(Manifold('name', 3))


def test_diofant__diffgeom__diffgeom__Patch():
    from diofant.diffgeom import Manifold, Patch
    assert _test_args(Patch('name', Manifold('name', 3)))


def test_diofant__diffgeom__diffgeom__CoordSystem():
    from diofant.diffgeom import Manifold, Patch, CoordSystem
    assert _test_args(CoordSystem('name', Patch('name', Manifold('name', 3))))


@pytest.mark.xfail
def test_diofant__diffgeom__diffgeom__Point():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, Point
    assert _test_args(Point(
        CoordSystem('name', Patch('name', Manifold('name', 3))), [x, y]))


def test_diofant__diffgeom__diffgeom__BaseScalarField():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseScalarField(cs, 0))


def test_diofant__diffgeom__diffgeom__BaseVectorField():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseVectorField
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseVectorField(cs, 0))


def test_diofant__diffgeom__diffgeom__Differential():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(Differential(BaseScalarField(cs, 0)))


def test_diofant__diffgeom__diffgeom__Commutator():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseVectorField, Commutator
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    cs1 = CoordSystem('name1', Patch('name', Manifold('name', 3)))
    v = BaseVectorField(cs, 0)
    v1 = BaseVectorField(cs1, 0)
    assert _test_args(Commutator(v, v1))


def test_diofant__diffgeom__diffgeom__TensorProduct():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential, TensorProduct
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    assert _test_args(TensorProduct(d, d))


def test_diofant__diffgeom__diffgeom__WedgeProduct():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential, WedgeProduct
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    d1 = Differential(BaseScalarField(cs, 1))
    assert _test_args(WedgeProduct(d, d1))


def test_diofant__diffgeom__diffgeom__LieDerivative():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential, BaseVectorField, LieDerivative
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    v = BaseVectorField(cs, 0)
    assert _test_args(LieDerivative(v, d))


@pytest.mark.xfail
def test_diofant__diffgeom__diffgeom__BaseCovarDerivativeOp():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseCovarDerivativeOp
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseCovarDerivativeOp(cs, 0, [[[0, ]*3, ]*3, ]*3))


def test_diofant__diffgeom__diffgeom__CovarDerivativeOp():
    from diofant.diffgeom import Manifold, Patch, CoordSystem, BaseVectorField, CovarDerivativeOp
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    v = BaseVectorField(cs, 0)
    _test_args(CovarDerivativeOp(v, [[[0, ]*3, ]*3, ]*3))


def test_diofant__ntheory__factor___totient():
    from diofant.ntheory.factor_ import totient
    k = symbols('k', integer=True)
    t = totient(k)
    assert _test_args(t)


def test_diofant__ntheory__factor___divisor_sigma():
    from diofant.ntheory.factor_ import divisor_sigma
    k = symbols('k', integer=True)
    n = symbols('n', integer=True)
    t = divisor_sigma(n, k)
    assert _test_args(t)


def test_diofant__ntheory__residue_ntheory__mobius():
    from diofant.ntheory import mobius
    assert _test_args(mobius(2))


def test_diofant__printing__codeprinter__Assignment():
    from diofant.printing.codeprinter import Assignment
    assert _test_args(Assignment(x, y))


def test_diofant__vector__coordsysrect__CoordSysCartesian():
    from diofant.vector.coordsysrect import CoordSysCartesian
    assert _test_args(CoordSysCartesian('C'))


def test_diofant__vector__point__Point():
    from diofant.vector.point import Point
    assert _test_args(Point('P'))


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
    from diofant.vector.vector import BaseVector
    from diofant.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(BaseVector('Ci', 0, C, ' ', ' '))


def test_diofant__vector__vector__VectorAdd():
    from diofant.vector.vector import VectorAdd, VectorMul
    from diofant.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    from diofant.abc import a, b, c, x, y, z
    v1 = a*C.i + b*C.j + c*C.k
    v2 = x*C.i + y*C.j + z*C.k
    assert _test_args(VectorAdd(v1, v2))
    assert _test_args(VectorMul(x, v1))


def test_diofant__vector__vector__VectorMul():
    from diofant.vector.vector import VectorMul
    from diofant.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    from diofant.abc import a
    assert _test_args(VectorMul(a, C.i))


def test_diofant__vector__vector__VectorZero():
    from diofant.vector.vector import VectorZero
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
    from diofant.vector.dyadic import BaseDyadic
    from diofant.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(BaseDyadic(C.i, C.j))


def test_diofant__vector__dyadic__DyadicMul():
    from diofant.vector.dyadic import BaseDyadic, DyadicMul
    from diofant.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(DyadicMul(3, BaseDyadic(C.i, C.j)))


def test_diofant__vector__dyadic__DyadicAdd():
    from diofant.vector.dyadic import BaseDyadic, DyadicAdd
    from diofant.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(2 * DyadicAdd(BaseDyadic(C.i, C.i),
                                    BaseDyadic(C.i, C.j)))


def test_diofant__vector__dyadic__DyadicZero():
    from diofant.vector.dyadic import DyadicZero
    assert _test_args(DyadicZero())


def test_diofant__vector__deloperator__Del():
    from diofant.vector.deloperator import Del
    from diofant.vector.coordsysrect import CoordSysCartesian
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
    from diofant.vector.orienters import AxisOrienter
    from diofant.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(AxisOrienter(x, C.i))


def test_diofant__vector__orienters__BodyOrienter():
    from diofant.vector.orienters import BodyOrienter
    assert _test_args(BodyOrienter(x, y, z, '123'))


def test_diofant__vector__orienters__SpaceOrienter():
    from diofant.vector.orienters import SpaceOrienter
    assert _test_args(SpaceOrienter(x, y, z, '123'))


def test_diofant__vector__orienters__QuaternionOrienter():
    from diofant.vector.orienters import QuaternionOrienter
    a, b, c, d = symbols('a b c d')
    assert _test_args(QuaternionOrienter(a, b, c, d))


def test_diofant__vector__scalar__BaseScalar():
    from diofant.vector.scalar import BaseScalar
    from diofant.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(BaseScalar('Cx', 0, C, ' ', ' '))
