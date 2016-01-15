"""Test whether all elements of cls.args are instances of Basic. """

# NOTE: keep tests sorted by (module, class name) key.

import os
import re
import warnings
import io

import pytest

from sympy import Basic, S, symbols, sqrt, sin, oo, Interval, exp, Integer
from sympy.utilities.exceptions import SymPyDeprecationWarning

x, y, z = symbols('x,y,z')


def test_all_classes_are_tested():
    this = os.path.split(__file__)[0]
    path = os.path.join(this, os.pardir, os.pardir)
    sympy_path = os.path.abspath(path)
    prefix = os.path.split(sympy_path)[0] + os.sep

    re_cls = re.compile("^class ([A-Za-z][A-Za-z0-9_]*)\s*\(", re.MULTILINE)

    modules = {}

    for root, dirs, files in os.walk(sympy_path):
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

    # reset all SymPyDeprecationWarning into errors
    warnings.simplefilter("error", category=SymPyDeprecationWarning)

    assert not failed, "Missing classes: %s.  Please add tests for these to sympy/core/tests/test_args.py." % ", ".join(failed)


def _test_args(obj):
    return all(isinstance(arg, Basic) for arg in obj.args)


def test_sympy__assumptions__assume__AppliedPredicate():
    from sympy.assumptions.assume import AppliedPredicate, Predicate
    assert _test_args(AppliedPredicate(Predicate("test"), 2))


def test_sympy__assumptions__assume__Predicate():
    from sympy.assumptions.assume import Predicate
    assert _test_args(Predicate("test"))


@pytest.mark.xfail
def test_sympy__combinatorics__graycode__GrayCode():
    from sympy.combinatorics.graycode import GrayCode
    # an integer is given and returned from GrayCode as the arg
    assert _test_args(GrayCode(3, start='100'))
    assert _test_args(GrayCode(3, rank=1))


def test_sympy__combinatorics__subsets__Subset():
    from sympy.combinatorics.subsets import Subset
    assert _test_args(Subset([0, 1], [0, 1, 2, 3]))
    assert _test_args(Subset(['c', 'd'], ['a', 'b', 'c', 'd']))


@pytest.mark.xfail
def test_sympy__combinatorics__permutations__Permutation():
    from sympy.combinatorics.permutations import Permutation
    assert _test_args(Permutation([0, 1, 2, 3]))


def test_sympy__combinatorics__perm_groups__PermutationGroup():
    from sympy.combinatorics.permutations import Permutation
    from sympy.combinatorics.perm_groups import PermutationGroup
    assert _test_args(PermutationGroup([Permutation([0, 1])]))


def test_sympy__combinatorics__polyhedron__Polyhedron():
    from sympy.combinatorics.permutations import Permutation
    from sympy.combinatorics.polyhedron import Polyhedron
    from sympy.abc import w, x, y, z
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
def test_sympy__combinatorics__prufer__Prufer():
    from sympy.combinatorics.prufer import Prufer
    assert _test_args(Prufer([[0, 1], [0, 2], [0, 3]], 4))


def test_sympy__combinatorics__partitions__Partition():
    from sympy.combinatorics.partitions import Partition
    assert _test_args(Partition([1]))


@pytest.mark.xfail
def test_sympy__combinatorics__partitions__IntegerPartition():
    from sympy.combinatorics.partitions import IntegerPartition
    assert _test_args(IntegerPartition([1]))


def test_sympy__concrete__products__Product():
    from sympy.concrete.products import Product
    assert _test_args(Product(x, (x, 0, 10)))
    assert _test_args(Product(x, (x, 0, y), (y, 0, 10)))


def test_sympy__concrete__expr_with_limits__ExprWithLimits():
    from sympy.concrete.expr_with_limits import ExprWithLimits
    assert _test_args(ExprWithLimits(x, (x, 0, 10)))
    assert _test_args(ExprWithLimits(x*y, (x, 0, 10.),(y,1.,3)))


def test_sympy__concrete__expr_with_limits__AddWithLimits():
    from sympy.concrete.expr_with_limits import AddWithLimits
    assert _test_args(AddWithLimits(x, (x, 0, 10)))
    assert _test_args(AddWithLimits(x*y, (x, 0, 10),(y,1,3)))


def test_sympy__concrete__expr_with_intlimits__ExprWithIntLimits():
    from sympy.concrete.expr_with_intlimits import ExprWithIntLimits
    assert _test_args(ExprWithIntLimits(x, (x, 0, 10)))
    assert _test_args(ExprWithIntLimits(x*y, (x, 0, 10),(y,1,3)))


def test_sympy__concrete__summations__Sum():
    from sympy.concrete.summations import Sum
    assert _test_args(Sum(x, (x, 0, 10)))
    assert _test_args(Sum(x, (x, 0, y), (y, 0, 10)))


def test_sympy__core__add__Add():
    from sympy.core.add import Add
    assert _test_args(Add(x, y, z, 2))


def test_sympy__core__basic__Atom():
    from sympy.core.basic import Atom
    assert _test_args(Atom())


def test_sympy__core__basic__Basic():
    from sympy.core.basic import Basic
    assert _test_args(Basic())


def test_sympy__core__containers__Dict():
    from sympy.core.containers import Dict
    assert _test_args(Dict({x: y, y: z}))


def test_sympy__core__containers__Tuple():
    from sympy.core.containers import Tuple
    assert _test_args(Tuple(x, y, z, 2))


def test_sympy__core__expr__AtomicExpr():
    from sympy.core.expr import AtomicExpr
    assert _test_args(AtomicExpr())


def test_sympy__core__expr__Expr():
    from sympy.core.expr import Expr
    assert _test_args(Expr())


def test_sympy__core__function__Application():
    from sympy.core.function import Application
    assert _test_args(Application(1, 2, 3))


def test_sympy__core__function__AppliedUndef():
    from sympy.core.function import AppliedUndef
    assert _test_args(AppliedUndef(1, 2, 3))


def test_sympy__core__function__Derivative():
    from sympy.core.function import Derivative
    assert _test_args(Derivative(2, x, y, 3))


def test_sympy__core__function__Function():
    pass


def test_sympy__core__function__Lambda():
    from sympy.core.function import Lambda
    assert _test_args(Lambda((x, y), x + y + z))


def test_sympy__core__function__Subs():
    from sympy.core.function import Subs
    assert _test_args(Subs(x + y, x, 2))


def test_sympy__core__function__WildFunction():
    from sympy.core.function import WildFunction
    assert _test_args(WildFunction('f'))


def test_sympy__core__mod__Mod():
    from sympy.core.mod import Mod
    assert _test_args(Mod(x, 2))


def test_sympy__core__mul__Mul():
    from sympy.core.mul import Mul
    assert _test_args(Mul(2, x, y, z))


def test_sympy__core__numbers__Catalan():
    from sympy.core.numbers import Catalan
    assert _test_args(Catalan())


def test_sympy__core__numbers__ComplexInfinity():
    from sympy.core.numbers import ComplexInfinity
    assert _test_args(ComplexInfinity())


def test_sympy__core__numbers__EulerGamma():
    from sympy.core.numbers import EulerGamma
    assert _test_args(EulerGamma())


def test_sympy__core__numbers__Exp1():
    from sympy.core.numbers import Exp1
    assert _test_args(Exp1())


def test_sympy__core__numbers__Float():
    from sympy.core.numbers import Float
    assert _test_args(Float(1.23))


def test_sympy__core__numbers__GoldenRatio():
    from sympy.core.numbers import GoldenRatio
    assert _test_args(GoldenRatio())


def test_sympy__core__numbers__Half():
    from sympy.core.numbers import Half
    assert _test_args(Half())


def test_sympy__core__numbers__ImaginaryUnit():
    from sympy.core.numbers import ImaginaryUnit
    assert _test_args(ImaginaryUnit())


def test_sympy__core__numbers__Infinity():
    from sympy.core.numbers import Infinity
    assert _test_args(Infinity())


def test_sympy__core__numbers__Integer():
    from sympy.core.numbers import Integer
    assert _test_args(Integer(7))


def test_sympy__core__numbers__IntegerConstant():
    pass


def test_sympy__core__numbers__NaN():
    from sympy.core.numbers import NaN
    assert _test_args(NaN())


def test_sympy__core__numbers__NegativeInfinity():
    from sympy.core.numbers import NegativeInfinity
    assert _test_args(NegativeInfinity())


def test_sympy__core__numbers__NegativeOne():
    from sympy.core.numbers import NegativeOne
    assert _test_args(NegativeOne())


def test_sympy__core__numbers__Number():
    from sympy.core.numbers import Number
    assert _test_args(Number(1, 7))


def test_sympy__core__numbers__NumberSymbol():
    from sympy.core.numbers import NumberSymbol
    assert _test_args(NumberSymbol())


def test_sympy__core__numbers__One():
    from sympy.core.numbers import One
    assert _test_args(One())


def test_sympy__core__numbers__Pi():
    from sympy.core.numbers import Pi
    assert _test_args(Pi())


def test_sympy__core__numbers__Rational():
    from sympy.core.numbers import Rational
    assert _test_args(Rational(1, 7))


def test_sympy__core__numbers__RationalConstant():
    pass


def test_sympy__core__numbers__Zero():
    from sympy.core.numbers import Zero
    assert _test_args(Zero())


def test_sympy__core__operations__AssocOp():
    pass


def test_sympy__core__operations__LatticeOp():
    pass


def test_sympy__core__power__Pow():
    from sympy.core.power import Pow
    assert _test_args(Pow(x, 2))


def test_sympy__core__relational__Equality():
    from sympy.core.relational import Equality
    assert _test_args(Equality(x, 2))


def test_sympy__core__relational__GreaterThan():
    from sympy.core.relational import GreaterThan
    assert _test_args(GreaterThan(x, 2))


def test_sympy__core__relational__LessThan():
    from sympy.core.relational import LessThan
    assert _test_args(LessThan(x, 2))


def test_sympy__core__relational__Relational():
    pass


def test_sympy__core__relational__StrictGreaterThan():
    from sympy.core.relational import StrictGreaterThan
    assert _test_args(StrictGreaterThan(x, 2))


def test_sympy__core__relational__StrictLessThan():
    from sympy.core.relational import StrictLessThan
    assert _test_args(StrictLessThan(x, 2))


def test_sympy__core__relational__Unequality():
    from sympy.core.relational import Unequality
    assert _test_args(Unequality(x, 2))


def test_sympy__sets__sets__EmptySet():
    from sympy.sets.sets import EmptySet
    assert _test_args(EmptySet())


def test_sympy__sets__sets__UniversalSet():
    from sympy.sets.sets import UniversalSet
    assert _test_args(UniversalSet())


def test_sympy__sets__sets__FiniteSet():
    from sympy.sets.sets import FiniteSet
    assert _test_args(FiniteSet(x, y, z))


def test_sympy__sets__sets__Interval():
    from sympy.sets.sets import Interval
    assert _test_args(Interval(0, 1))


def test_sympy__sets__sets__ProductSet():
    from sympy.sets.sets import ProductSet, Interval
    assert _test_args(ProductSet(Interval(0, 1), Interval(0, 1)))


def test_sympy__sets__sets__Set():
    from sympy.sets.sets import Set
    assert _test_args(Set())


def test_sympy__sets__sets__Intersection():
    from sympy.sets.sets import Intersection, Interval
    assert _test_args(Intersection(Interval(0, 3), Interval(2, 4),
        evaluate=False))


def test_sympy__sets__sets__Union():
    from sympy.sets.sets import Union, Interval
    assert _test_args(Union(Interval(0, 1), Interval(2, 3)))


def test_sympy__sets__sets__Complement():
    from sympy.sets.sets import Complement
    assert _test_args(Complement(Interval(0, 2), Interval(0, 1)))


def test_sympy__sets__sets__SymmetricDifference():
    from sympy.sets.sets import FiniteSet, SymmetricDifference
    assert _test_args(SymmetricDifference(FiniteSet(1, 2, 3),
           FiniteSet(2, 3, 4)))


def test_sympy__core__trace__Tr():
    from sympy.core.trace import Tr
    a, b = symbols('a b')
    assert _test_args(Tr(a + b))


def test_sympy__sets__fancysets__Naturals():
    from sympy.sets.fancysets import Naturals
    assert _test_args(Naturals())


def test_sympy__sets__fancysets__Naturals0():
    from sympy.sets.fancysets import Naturals0
    assert _test_args(Naturals0())


def test_sympy__sets__fancysets__Integers():
    from sympy.sets.fancysets import Integers
    assert _test_args(Integers())


def test_sympy__sets__fancysets__Reals():
    from sympy.sets.fancysets import Reals
    assert _test_args(Reals())


def test_sympy__sets__fancysets__ImageSet():
    from sympy.sets.fancysets import ImageSet
    from sympy import S, Lambda, Symbol
    x = Symbol('x')
    assert _test_args(ImageSet(Lambda(x, x**2), S.Naturals))


def test_sympy__sets__fancysets__Range():
    from sympy.sets.fancysets import Range
    assert _test_args(Range(1, 5, 1))


def test_sympy__sets__contains__Contains():
    from sympy.sets.fancysets import Range
    from sympy.sets.contains import Contains
    assert _test_args(Contains(x, Range(0, 10, 2)))


# STATS


from sympy.stats.crv_types import NormalDistribution
nd = NormalDistribution(0, 1)
from sympy.stats.frv_types import DieDistribution
die = DieDistribution(6)


def test_sympy__stats__crv__ContinuousDomain():
    from sympy.stats.crv import ContinuousDomain
    assert _test_args(ContinuousDomain({x}, Interval(-oo, oo)))


def test_sympy__stats__crv__SingleContinuousDomain():
    from sympy.stats.crv import SingleContinuousDomain
    assert _test_args(SingleContinuousDomain(x, Interval(-oo, oo)))


def test_sympy__stats__crv__ProductContinuousDomain():
    from sympy.stats.crv import SingleContinuousDomain, ProductContinuousDomain
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    E = SingleContinuousDomain(y, Interval(0, oo))
    assert _test_args(ProductContinuousDomain(D, E))


def test_sympy__stats__crv__ConditionalContinuousDomain():
    from sympy.stats.crv import (SingleContinuousDomain,
            ConditionalContinuousDomain)
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    assert _test_args(ConditionalContinuousDomain(D, x > 0))


def test_sympy__stats__crv__ContinuousPSpace():
    from sympy.stats.crv import ContinuousPSpace, SingleContinuousDomain
    D = SingleContinuousDomain(x, Interval(-oo, oo))
    assert _test_args(ContinuousPSpace(D, nd))


def test_sympy__stats__crv__SingleContinuousPSpace():
    from sympy.stats.crv import SingleContinuousPSpace
    assert _test_args(SingleContinuousPSpace(x, nd))


def test_sympy__stats__crv__ProductContinuousPSpace():
    from sympy.stats.crv import ProductContinuousPSpace, SingleContinuousPSpace
    A = SingleContinuousPSpace(x, nd)
    B = SingleContinuousPSpace(y, nd)
    assert _test_args(ProductContinuousPSpace(A, B))


def test_sympy__stats__crv__SingleContinuousDistribution():
    pass


def test_sympy__stats__drv__SingleDiscreteDomain():
    from sympy.stats.drv import SingleDiscreteDomain
    assert _test_args(SingleDiscreteDomain(x, S.Naturals))


def test_sympy__stats__drv__SingleDiscretePSpace():
    from sympy.stats.drv import SingleDiscretePSpace
    from sympy.stats.drv_types import PoissonDistribution
    assert _test_args(SingleDiscretePSpace(x, PoissonDistribution(1)))


def test_sympy__stats__drv__SingleDiscreteDistribution():
    pass


def test_sympy__stats__rv__RandomDomain():
    from sympy.stats.rv import RandomDomain
    from sympy.sets.sets import FiniteSet
    assert _test_args(RandomDomain(FiniteSet(x), FiniteSet(1, 2, 3)))


def test_sympy__stats__rv__SingleDomain():
    from sympy.stats.rv import SingleDomain
    from sympy.sets.sets import FiniteSet
    assert _test_args(SingleDomain(x, FiniteSet(1, 2, 3)))


def test_sympy__stats__rv__ConditionalDomain():
    from sympy.stats.rv import ConditionalDomain, RandomDomain
    from sympy.sets.sets import FiniteSet
    D = RandomDomain(FiniteSet(x), FiniteSet(1, 2))
    assert _test_args(ConditionalDomain(D, x > 1))


def test_sympy__stats__rv__PSpace():
    from sympy.stats.rv import PSpace, RandomDomain
    from sympy import FiniteSet
    D = RandomDomain(FiniteSet(x), FiniteSet(1, 2, 3, 4, 5, 6))
    assert _test_args(PSpace(D, die))


def test_sympy__stats__rv__SinglePSpace():
    pass


def test_sympy__stats__rv__RandomSymbol():
    from sympy.stats.rv import RandomSymbol
    from sympy.stats.crv import SingleContinuousPSpace
    A = SingleContinuousPSpace(x, nd)
    assert _test_args(RandomSymbol(A, x))


def test_sympy__stats__rv__ProductPSpace():
    from sympy.stats.rv import ProductPSpace
    from sympy.stats.crv import SingleContinuousPSpace
    A = SingleContinuousPSpace(x, nd)
    B = SingleContinuousPSpace(y, nd)
    assert _test_args(ProductPSpace(A, B))


def test_sympy__stats__rv__ProductDomain():
    from sympy.stats.rv import ProductDomain, SingleDomain
    D = SingleDomain(x, Interval(-oo, oo))
    E = SingleDomain(y, Interval(0, oo))
    assert _test_args(ProductDomain(D, E))


def test_sympy__stats__frv_types__DiscreteUniformDistribution():
    from sympy.stats.frv_types import DiscreteUniformDistribution
    from sympy.core.containers import Tuple
    assert _test_args(DiscreteUniformDistribution(Tuple(*range(6))))


def test_sympy__stats__frv_types__DieDistribution():
    from sympy.stats.frv_types import DieDistribution
    assert _test_args(DieDistribution(6))


def test_sympy__stats__frv_types__BernoulliDistribution():
    from sympy.stats.frv_types import BernoulliDistribution
    assert _test_args(BernoulliDistribution(S.Half, 0, 1))


def test_sympy__stats__frv_types__BinomialDistribution():
    from sympy.stats.frv_types import BinomialDistribution
    assert _test_args(BinomialDistribution(5, S.Half, 1, 0))


def test_sympy__stats__frv_types__HypergeometricDistribution():
    from sympy.stats.frv_types import HypergeometricDistribution
    assert _test_args(HypergeometricDistribution(10, 5, 3))


def test_sympy__stats__frv_types__RademacherDistribution():
    from sympy.stats.frv_types import RademacherDistribution
    assert _test_args(RademacherDistribution())


def test_sympy__stats__frv__FiniteDomain():
    from sympy.stats.frv import FiniteDomain
    assert _test_args(FiniteDomain({(x, 1), (x, 2)}))  # x can be 1 or 2


def test_sympy__stats__frv__SingleFiniteDomain():
    from sympy.stats.frv import SingleFiniteDomain
    assert _test_args(SingleFiniteDomain(x, {1, 2}))  # x can be 1 or 2


def test_sympy__stats__frv__ProductFiniteDomain():
    from sympy.stats.frv import SingleFiniteDomain, ProductFiniteDomain
    xd = SingleFiniteDomain(x, {1, 2})
    yd = SingleFiniteDomain(y, {1, 2})
    assert _test_args(ProductFiniteDomain(xd, yd))


def test_sympy__stats__frv__ConditionalFiniteDomain():
    from sympy.stats.frv import SingleFiniteDomain, ConditionalFiniteDomain
    xd = SingleFiniteDomain(x, {1, 2})
    assert _test_args(ConditionalFiniteDomain(xd, x > 1))


def test_sympy__stats__frv__FinitePSpace():
    from sympy.stats.frv import FinitePSpace, SingleFiniteDomain
    xd = SingleFiniteDomain(x, {1, 2, 3, 4, 5, 6})
    p = 1.0/6
    xd = SingleFiniteDomain(x, {1, 2})
    assert _test_args(FinitePSpace(xd, {(x, 1): S.Half, (x, 2): S.Half}))


def test_sympy__stats__frv__SingleFinitePSpace():
    from sympy.stats.frv import SingleFinitePSpace
    from sympy import Symbol

    assert _test_args(SingleFinitePSpace(Symbol('x'), die))


def test_sympy__stats__frv__ProductFinitePSpace():
    from sympy.stats.frv import SingleFinitePSpace, ProductFinitePSpace
    from sympy import Symbol
    xp = SingleFinitePSpace(Symbol('x'), die)
    yp = SingleFinitePSpace(Symbol('y'), die)
    assert _test_args(ProductFinitePSpace(xp, yp))


def test_sympy__stats__frv__SingleFiniteDistribution():
    pass


def test_sympy__stats__crv__ContinuousDistribution():
    pass


def test_sympy__stats__frv_types__FiniteDistributionHandmade():
    from sympy.stats.frv_types import FiniteDistributionHandmade
    assert _test_args(FiniteDistributionHandmade({1: 1}))


def test_sympy__stats__crv__ContinuousDistributionHandmade():
    from sympy.stats.crv import ContinuousDistributionHandmade
    from sympy import Symbol, Interval
    assert _test_args(ContinuousDistributionHandmade(Symbol('x'),
                                                     Interval(0, 2)))


def test_sympy__stats__rv__Density():
    from sympy.stats.rv import Density
    from sympy.stats.crv_types import Normal
    assert _test_args(Density(Normal('x', 0, 1)))


def test_sympy__stats__crv_types__ArcsinDistribution():
    from sympy.stats.crv_types import ArcsinDistribution
    assert _test_args(ArcsinDistribution(0, 1))


def test_sympy__stats__crv_types__BeniniDistribution():
    from sympy.stats.crv_types import BeniniDistribution
    assert _test_args(BeniniDistribution(1, 1, 1))


def test_sympy__stats__crv_types__BetaDistribution():
    from sympy.stats.crv_types import BetaDistribution
    assert _test_args(BetaDistribution(1, 1))


def test_sympy__stats__crv_types__BetaPrimeDistribution():
    from sympy.stats.crv_types import BetaPrimeDistribution
    assert _test_args(BetaPrimeDistribution(1, 1))


def test_sympy__stats__crv_types__CauchyDistribution():
    from sympy.stats.crv_types import CauchyDistribution
    assert _test_args(CauchyDistribution(0, 1))


def test_sympy__stats__crv_types__ChiDistribution():
    from sympy.stats.crv_types import ChiDistribution
    assert _test_args(ChiDistribution(1))


def test_sympy__stats__crv_types__ChiNoncentralDistribution():
    from sympy.stats.crv_types import ChiNoncentralDistribution
    assert _test_args(ChiNoncentralDistribution(1,1))


def test_sympy__stats__crv_types__ChiSquaredDistribution():
    from sympy.stats.crv_types import ChiSquaredDistribution
    assert _test_args(ChiSquaredDistribution(1))


def test_sympy__stats__crv_types__DagumDistribution():
    from sympy.stats.crv_types import DagumDistribution
    assert _test_args(DagumDistribution(1, 1, 1))


def test_sympy__stats__crv_types__ExponentialDistribution():
    from sympy.stats.crv_types import ExponentialDistribution
    assert _test_args(ExponentialDistribution(1))


def test_sympy__stats__crv_types__FDistributionDistribution():
    from sympy.stats.crv_types import FDistributionDistribution
    assert _test_args(FDistributionDistribution(1, 1))


def test_sympy__stats__crv_types__FisherZDistribution():
    from sympy.stats.crv_types import FisherZDistribution
    assert _test_args(FisherZDistribution(1, 1))


def test_sympy__stats__crv_types__FrechetDistribution():
    from sympy.stats.crv_types import FrechetDistribution
    assert _test_args(FrechetDistribution(1, 1, 1))


def test_sympy__stats__crv_types__GammaInverseDistribution():
    from sympy.stats.crv_types import GammaInverseDistribution
    assert _test_args(GammaInverseDistribution(1, 1))


def test_sympy__stats__crv_types__GammaDistribution():
    from sympy.stats.crv_types import GammaDistribution
    assert _test_args(GammaDistribution(1, 1))


def test_sympy__stats__crv_types__KumaraswamyDistribution():
    from sympy.stats.crv_types import KumaraswamyDistribution
    assert _test_args(KumaraswamyDistribution(1, 1))


def test_sympy__stats__crv_types__LaplaceDistribution():
    from sympy.stats.crv_types import LaplaceDistribution
    assert _test_args(LaplaceDistribution(0, 1))


def test_sympy__stats__crv_types__LogisticDistribution():
    from sympy.stats.crv_types import LogisticDistribution
    assert _test_args(LogisticDistribution(0, 1))


def test_sympy__stats__crv_types__LogNormalDistribution():
    from sympy.stats.crv_types import LogNormalDistribution
    assert _test_args(LogNormalDistribution(0, 1))


def test_sympy__stats__crv_types__MaxwellDistribution():
    from sympy.stats.crv_types import MaxwellDistribution
    assert _test_args(MaxwellDistribution(1))


def test_sympy__stats__crv_types__NakagamiDistribution():
    from sympy.stats.crv_types import NakagamiDistribution
    assert _test_args(NakagamiDistribution(1, 1))


def test_sympy__stats__crv_types__NormalDistribution():
    from sympy.stats.crv_types import NormalDistribution
    assert _test_args(NormalDistribution(0, 1))


def test_sympy__stats__crv_types__ParetoDistribution():
    from sympy.stats.crv_types import ParetoDistribution
    assert _test_args(ParetoDistribution(1, 1))


def test_sympy__stats__crv_types__QuadraticUDistribution():
    from sympy.stats.crv_types import QuadraticUDistribution
    assert _test_args(QuadraticUDistribution(1, 2))


def test_sympy__stats__crv_types__RaisedCosineDistribution():
    from sympy.stats.crv_types import RaisedCosineDistribution
    assert _test_args(RaisedCosineDistribution(1, 1))


def test_sympy__stats__crv_types__RayleighDistribution():
    from sympy.stats.crv_types import RayleighDistribution
    assert _test_args(RayleighDistribution(1))


def test_sympy__stats__crv_types__StudentTDistribution():
    from sympy.stats.crv_types import StudentTDistribution
    assert _test_args(StudentTDistribution(1))


def test_sympy__stats__crv_types__TriangularDistribution():
    from sympy.stats.crv_types import TriangularDistribution
    assert _test_args(TriangularDistribution(-1, 0, 1))


def test_sympy__stats__crv_types__UniformDistribution():
    from sympy.stats.crv_types import UniformDistribution
    assert _test_args(UniformDistribution(0, 1))


def test_sympy__stats__crv_types__UniformSumDistribution():
    from sympy.stats.crv_types import UniformSumDistribution
    assert _test_args(UniformSumDistribution(1))


def test_sympy__stats__crv_types__VonMisesDistribution():
    from sympy.stats.crv_types import VonMisesDistribution
    assert _test_args(VonMisesDistribution(1, 1))


def test_sympy__stats__crv_types__WeibullDistribution():
    from sympy.stats.crv_types import WeibullDistribution
    assert _test_args(WeibullDistribution(1, 1))


def test_sympy__stats__crv_types__WignerSemicircleDistribution():
    from sympy.stats.crv_types import WignerSemicircleDistribution
    assert _test_args(WignerSemicircleDistribution(1))


def test_sympy__stats__drv_types__PoissonDistribution():
    from sympy.stats.drv_types import PoissonDistribution
    assert _test_args(PoissonDistribution(1))


def test_sympy__stats__drv_types__GeometricDistribution():
    from sympy.stats.drv_types import GeometricDistribution
    assert _test_args(GeometricDistribution(.5))


def test_sympy__core__symbol__Dummy():
    from sympy.core.symbol import Dummy
    assert _test_args(Dummy('t'))


def test_sympy__core__symbol__Symbol():
    from sympy.core.symbol import Symbol
    assert _test_args(Symbol('t'))


def test_sympy__core__symbol__Wild():
    from sympy.core.symbol import Wild
    assert _test_args(Wild('x', exclude=[x]))


def test_sympy__functions__combinatorial__factorials__CombinatorialFunction():
    pass


def test_sympy__functions__combinatorial__factorials__FallingFactorial():
    from sympy.functions.combinatorial.factorials import FallingFactorial
    assert _test_args(FallingFactorial(2, x))


def test_sympy__functions__combinatorial__factorials__MultiFactorial():
    from sympy.functions.combinatorial.factorials import MultiFactorial
    assert _test_args(MultiFactorial(x))


def test_sympy__functions__combinatorial__factorials__RisingFactorial():
    from sympy.functions.combinatorial.factorials import RisingFactorial
    assert _test_args(RisingFactorial(2, x))


def test_sympy__functions__combinatorial__factorials__binomial():
    from sympy.functions.combinatorial.factorials import binomial
    assert _test_args(binomial(2, x))


def test_sympy__functions__combinatorial__factorials__subfactorial():
    from sympy.functions.combinatorial.factorials import subfactorial
    assert _test_args(subfactorial(1))


def test_sympy__functions__combinatorial__factorials__factorial():
    from sympy.functions.combinatorial.factorials import factorial
    assert _test_args(factorial(x))


def test_sympy__functions__combinatorial__factorials__factorial2():
    from sympy.functions.combinatorial.factorials import factorial2
    assert _test_args(factorial2(x))


def test_sympy__functions__combinatorial__numbers__bell():
    from sympy.functions.combinatorial.numbers import bell
    assert _test_args(bell(x, y))


def test_sympy__functions__combinatorial__numbers__bernoulli():
    from sympy.functions.combinatorial.numbers import bernoulli
    assert _test_args(bernoulli(x))


def test_sympy__functions__combinatorial__numbers__catalan():
    from sympy.functions.combinatorial.numbers import catalan
    assert _test_args(catalan(x))


def test_sympy__functions__combinatorial__numbers__genocchi():
    from sympy.functions.combinatorial.numbers import genocchi
    assert _test_args(genocchi(x))


def test_sympy__functions__combinatorial__numbers__euler():
    from sympy.functions.combinatorial.numbers import euler
    assert _test_args(euler(x))


def test_sympy__functions__combinatorial__numbers__fibonacci():
    from sympy.functions.combinatorial.numbers import fibonacci
    assert _test_args(fibonacci(x))


def test_sympy__functions__combinatorial__numbers__harmonic():
    from sympy.functions.combinatorial.numbers import harmonic
    assert _test_args(harmonic(x, 2))


def test_sympy__functions__combinatorial__numbers__lucas():
    from sympy.functions.combinatorial.numbers import lucas
    assert _test_args(lucas(x))


def test_sympy__functions__elementary__complexes__Abs():
    from sympy.functions.elementary.complexes import Abs
    assert _test_args(Abs(x))


def test_sympy__functions__elementary__complexes__adjoint():
    from sympy.functions.elementary.complexes import adjoint
    assert _test_args(adjoint(x))


def test_sympy__functions__elementary__complexes__arg():
    from sympy.functions.elementary.complexes import arg
    assert _test_args(arg(x))


def test_sympy__functions__elementary__complexes__conjugate():
    from sympy.functions.elementary.complexes import conjugate
    assert _test_args(conjugate(x))


def test_sympy__functions__elementary__complexes__im():
    from sympy.functions.elementary.complexes import im
    assert _test_args(im(x))


def test_sympy__functions__elementary__complexes__re():
    from sympy.functions.elementary.complexes import re
    assert _test_args(re(x))


def test_sympy__functions__elementary__complexes__sign():
    from sympy.functions.elementary.complexes import sign
    assert _test_args(sign(x))


def test_sympy__functions__elementary__complexes__polar_lift():
    from sympy.functions.elementary.complexes import polar_lift
    assert _test_args(polar_lift(x))


def test_sympy__functions__elementary__complexes__periodic_argument():
    from sympy.functions.elementary.complexes import periodic_argument
    assert _test_args(periodic_argument(x, y))


def test_sympy__functions__elementary__complexes__principal_branch():
    from sympy.functions.elementary.complexes import principal_branch
    assert _test_args(principal_branch(x, y))


def test_sympy__functions__elementary__complexes__transpose():
    from sympy.functions.elementary.complexes import transpose
    assert _test_args(transpose(x))


def test_sympy__functions__elementary__exponential__LambertW():
    from sympy.functions.elementary.exponential import LambertW
    assert _test_args(LambertW(2))


def test_sympy__functions__elementary__exponential__exp():
    from sympy.functions.elementary.exponential import exp
    assert _test_args(exp(2))


def test_sympy__functions__elementary__exponential__exp_polar():
    from sympy.functions.elementary.exponential import exp_polar
    assert _test_args(exp_polar(2))


def test_sympy__functions__elementary__exponential__log():
    from sympy.functions.elementary.exponential import log
    assert _test_args(log(2))


def test_sympy__functions__elementary__hyperbolic__HyperbolicFunction():
    pass


def test_sympy__functions__elementary__hyperbolic__ReciprocalHyperbolicFunction():
    pass


def test_sympy__functions__elementary__hyperbolic__acosh():
    from sympy.functions.elementary.hyperbolic import acosh
    assert _test_args(acosh(2))


def test_sympy__functions__elementary__hyperbolic__acoth():
    from sympy.functions.elementary.hyperbolic import acoth
    assert _test_args(acoth(2))


def test_sympy__functions__elementary__hyperbolic__asinh():
    from sympy.functions.elementary.hyperbolic import asinh
    assert _test_args(asinh(2))


def test_sympy__functions__elementary__hyperbolic__atanh():
    from sympy.functions.elementary.hyperbolic import atanh
    assert _test_args(atanh(2))


def test_sympy__functions__elementary__hyperbolic__cosh():
    from sympy.functions.elementary.hyperbolic import cosh
    assert _test_args(cosh(2))


def test_sympy__functions__elementary__hyperbolic__coth():
    from sympy.functions.elementary.hyperbolic import coth
    assert _test_args(coth(2))


def test_sympy__functions__elementary__hyperbolic__csch():
    from sympy.functions.elementary.hyperbolic import csch
    assert _test_args(csch(2))


def test_sympy__functions__elementary__hyperbolic__sech():
    from sympy.functions.elementary.hyperbolic import sech
    assert _test_args(sech(2))


def test_sympy__functions__elementary__hyperbolic__sinh():
    from sympy.functions.elementary.hyperbolic import sinh
    assert _test_args(sinh(2))


def test_sympy__functions__elementary__hyperbolic__tanh():
    from sympy.functions.elementary.hyperbolic import tanh
    assert _test_args(tanh(2))


@pytest.mark.xfail
def test_sympy__functions__elementary__integers__RoundFunction():
    from sympy.functions.elementary.integers import RoundFunction
    assert _test_args(RoundFunction())


def test_sympy__functions__elementary__integers__ceiling():
    from sympy.functions.elementary.integers import ceiling
    assert _test_args(ceiling(x))


def test_sympy__functions__elementary__integers__floor():
    from sympy.functions.elementary.integers import floor
    assert _test_args(floor(x))


def test_sympy__functions__elementary__miscellaneous__IdentityFunction():
    from sympy.functions.elementary.miscellaneous import IdentityFunction
    assert _test_args(IdentityFunction())


def test_sympy__functions__elementary__miscellaneous__Max():
    from sympy.functions.elementary.miscellaneous import Max
    assert _test_args(Max(x, 2))


def test_sympy__functions__elementary__miscellaneous__Min():
    from sympy.functions.elementary.miscellaneous import Min
    assert _test_args(Min(x, 2))


def test_sympy__functions__elementary__miscellaneous__MinMaxBase():
    pass


def test_sympy__functions__elementary__piecewise__ExprCondPair():
    from sympy.functions.elementary.piecewise import ExprCondPair
    assert _test_args(ExprCondPair(1, True))


def test_sympy__functions__elementary__piecewise__Piecewise():
    from sympy.functions.elementary.piecewise import Piecewise
    assert _test_args(Piecewise((1, x >= 0), (0, True)))


def test_sympy__functions__elementary__trigonometric__TrigonometricFunction():
    pass


def test_sympy__functions__elementary__trigonometric__ReciprocalTrigonometricFunction():
    pass


def test_sympy__functions__elementary__trigonometric__InverseTrigonometricFunction():
    pass


def test_sympy__functions__elementary__trigonometric__acos():
    from sympy.functions.elementary.trigonometric import acos
    assert _test_args(acos(2))


def test_sympy__functions__elementary__trigonometric__acot():
    from sympy.functions.elementary.trigonometric import acot
    assert _test_args(acot(2))


def test_sympy__functions__elementary__trigonometric__asin():
    from sympy.functions.elementary.trigonometric import asin
    assert _test_args(asin(2))


def test_sympy__functions__elementary__trigonometric__asec():
    from sympy.functions.elementary.trigonometric import asec
    assert _test_args(asec(2))


def test_sympy__functions__elementary__trigonometric__acsc():
    from sympy.functions.elementary.trigonometric import acsc
    assert _test_args(acsc(2))


def test_sympy__functions__elementary__trigonometric__atan():
    from sympy.functions.elementary.trigonometric import atan
    assert _test_args(atan(2))


def test_sympy__functions__elementary__trigonometric__atan2():
    from sympy.functions.elementary.trigonometric import atan2
    assert _test_args(atan2(2, 3))


def test_sympy__functions__elementary__trigonometric__cos():
    from sympy.functions.elementary.trigonometric import cos
    assert _test_args(cos(2))


def test_sympy__functions__elementary__trigonometric__csc():
    from sympy.functions.elementary.trigonometric import csc
    assert _test_args(csc(2))


def test_sympy__functions__elementary__trigonometric__cot():
    from sympy.functions.elementary.trigonometric import cot
    assert _test_args(cot(2))


def test_sympy__functions__elementary__trigonometric__sin():
    assert _test_args(sin(2))


def test_sympy__functions__elementary__trigonometric__sec():
    from sympy.functions.elementary.trigonometric import sec
    assert _test_args(sec(2))


def test_sympy__functions__elementary__trigonometric__tan():
    from sympy.functions.elementary.trigonometric import tan
    assert _test_args(tan(2))


def test_sympy__functions__special__bessel__BesselBase():
    pass


def test_sympy__functions__special__bessel__SphericalBesselBase():
    pass


def test_sympy__functions__special__bessel__besseli():
    from sympy.functions.special.bessel import besseli
    assert _test_args(besseli(x, 1))


def test_sympy__functions__special__bessel__besselj():
    from sympy.functions.special.bessel import besselj
    assert _test_args(besselj(x, 1))


def test_sympy__functions__special__bessel__besselk():
    from sympy.functions.special.bessel import besselk
    assert _test_args(besselk(x, 1))


def test_sympy__functions__special__bessel__bessely():
    from sympy.functions.special.bessel import bessely
    assert _test_args(bessely(x, 1))


def test_sympy__functions__special__bessel__hankel1():
    from sympy.functions.special.bessel import hankel1
    assert _test_args(hankel1(x, 1))


def test_sympy__functions__special__bessel__hankel2():
    from sympy.functions.special.bessel import hankel2
    assert _test_args(hankel2(x, 1))


def test_sympy__functions__special__bessel__jn():
    from sympy.functions.special.bessel import jn
    assert _test_args(jn(0, x))


def test_sympy__functions__special__bessel__yn():
    from sympy.functions.special.bessel import yn
    assert _test_args(yn(0, x))


def test_sympy__functions__special__bessel__AiryBase():
    pass


def test_sympy__functions__special__bessel__airyai():
    from sympy.functions.special.bessel import airyai
    assert _test_args(airyai(2))


def test_sympy__functions__special__bessel__airybi():
    from sympy.functions.special.bessel import airybi
    assert _test_args(airybi(2))


def test_sympy__functions__special__bessel__airyaiprime():
    from sympy.functions.special.bessel import airyaiprime
    assert _test_args(airyaiprime(2))


def test_sympy__functions__special__bessel__airybiprime():
    from sympy.functions.special.bessel import airybiprime
    assert _test_args(airybiprime(2))


def test_sympy__functions__special__elliptic_integrals__elliptic_k():
    from sympy.functions.special.elliptic_integrals import elliptic_k as K
    assert _test_args(K(x))


def test_sympy__functions__special__elliptic_integrals__elliptic_f():
    from sympy.functions.special.elliptic_integrals import elliptic_f as F
    assert _test_args(F(x, y))


def test_sympy__functions__special__elliptic_integrals__elliptic_e():
    from sympy.functions.special.elliptic_integrals import elliptic_e as E
    assert _test_args(E(x))
    assert _test_args(E(x, y))


def test_sympy__functions__special__elliptic_integrals__elliptic_pi():
    from sympy.functions.special.elliptic_integrals import elliptic_pi as P
    assert _test_args(P(x, y))
    assert _test_args(P(x, y, z))


def test_sympy__functions__special__delta_functions__DiracDelta():
    from sympy.functions.special.delta_functions import DiracDelta
    assert _test_args(DiracDelta(x, 1))


def test_sympy__functions__special__delta_functions__Heaviside():
    from sympy.functions.special.delta_functions import Heaviside
    assert _test_args(Heaviside(x))


def test_sympy__functions__special__error_functions__erf():
    from sympy.functions.special.error_functions import erf
    assert _test_args(erf(2))


def test_sympy__functions__special__error_functions__erfc():
    from sympy.functions.special.error_functions import erfc
    assert _test_args(erfc(2))


def test_sympy__functions__special__error_functions__erfi():
    from sympy.functions.special.error_functions import erfi
    assert _test_args(erfi(2))


def test_sympy__functions__special__error_functions__erf2():
    from sympy.functions.special.error_functions import erf2
    assert _test_args(erf2(2, 3))


def test_sympy__functions__special__error_functions__erfinv():
    from sympy.functions.special.error_functions import erfinv
    assert _test_args(erfinv(2))


def test_sympy__functions__special__error_functions__erfcinv():
    from sympy.functions.special.error_functions import erfcinv
    assert _test_args(erfcinv(2))


def test_sympy__functions__special__error_functions__erf2inv():
    from sympy.functions.special.error_functions import erf2inv
    assert _test_args(erf2inv(2, 3))


def test_sympy__functions__special__error_functions__FresnelIntegral():
    pass


def test_sympy__functions__special__error_functions__fresnels():
    from sympy.functions.special.error_functions import fresnels
    assert _test_args(fresnels(2))


def test_sympy__functions__special__error_functions__fresnelc():
    from sympy.functions.special.error_functions import fresnelc
    assert _test_args(fresnelc(2))


def test_sympy__functions__special__error_functions__erfs():
    from sympy.functions.special.error_functions import _erfs
    assert _test_args(_erfs(2))


def test_sympy__functions__special__error_functions__Ei():
    from sympy.functions.special.error_functions import Ei
    assert _test_args(Ei(2))


def test_sympy__functions__special__error_functions__li():
    from sympy.functions.special.error_functions import li
    assert _test_args(li(2))


def test_sympy__functions__special__error_functions__Li():
    from sympy.functions.special.error_functions import Li
    assert _test_args(Li(2))


def test_sympy__functions__special__error_functions__TrigonometricIntegral():
    pass


def test_sympy__functions__special__error_functions__Si():
    from sympy.functions.special.error_functions import Si
    assert _test_args(Si(2))


def test_sympy__functions__special__error_functions__Ci():
    from sympy.functions.special.error_functions import Ci
    assert _test_args(Ci(2))


def test_sympy__functions__special__error_functions__Shi():
    from sympy.functions.special.error_functions import Shi
    assert _test_args(Shi(2))


def test_sympy__functions__special__error_functions__Chi():
    from sympy.functions.special.error_functions import Chi
    assert _test_args(Chi(2))


def test_sympy__functions__special__error_functions__expint():
    from sympy.functions.special.error_functions import expint
    assert _test_args(expint(y, x))


def test_sympy__functions__special__gamma_functions__gamma():
    from sympy.functions.special.gamma_functions import gamma
    assert _test_args(gamma(x))


def test_sympy__functions__special__gamma_functions__loggamma():
    from sympy.functions.special.gamma_functions import loggamma
    assert _test_args(loggamma(2))


def test_sympy__functions__special__gamma_functions__lowergamma():
    from sympy.functions.special.gamma_functions import lowergamma
    assert _test_args(lowergamma(x, 2))


def test_sympy__functions__special__gamma_functions__polygamma():
    from sympy.functions.special.gamma_functions import polygamma
    assert _test_args(polygamma(x, 2))


def test_sympy__functions__special__gamma_functions__uppergamma():
    from sympy.functions.special.gamma_functions import uppergamma
    assert _test_args(uppergamma(x, 2))


def test_sympy__functions__special__beta_functions__beta():
    from sympy.functions.special.beta_functions import beta
    assert _test_args(beta(x, x))


def test_sympy__functions__special__hyper__TupleParametersBase():
    pass


def test_sympy__functions__special__hyper__TupleArg():
    pass


def test_sympy__functions__special__hyper__hyper():
    from sympy.functions.special.hyper import hyper
    assert _test_args(hyper([1, 2, 3], [4, 5], x))


def test_sympy__functions__special__hyper__meijerg():
    from sympy.functions.special.hyper import meijerg
    assert _test_args(meijerg([1, 2, 3], [4, 5], [6], [], x))


def test_sympy__functions__special__hyper__HyperRep():
    pass


def test_sympy__functions__special__hyper__HyperRep_power1():
    from sympy.functions.special.hyper import HyperRep_power1
    assert _test_args(HyperRep_power1(x, y))


def test_sympy__functions__special__hyper__HyperRep_power2():
    from sympy.functions.special.hyper import HyperRep_power2
    assert _test_args(HyperRep_power2(x, y))


def test_sympy__functions__special__hyper__HyperRep_log1():
    from sympy.functions.special.hyper import HyperRep_log1
    assert _test_args(HyperRep_log1(x))


def test_sympy__functions__special__hyper__HyperRep_atanh():
    from sympy.functions.special.hyper import HyperRep_atanh
    assert _test_args(HyperRep_atanh(x))


def test_sympy__functions__special__hyper__HyperRep_asin1():
    from sympy.functions.special.hyper import HyperRep_asin1
    assert _test_args(HyperRep_asin1(x))


def test_sympy__functions__special__hyper__HyperRep_asin2():
    from sympy.functions.special.hyper import HyperRep_asin2
    assert _test_args(HyperRep_asin2(x))


def test_sympy__functions__special__hyper__HyperRep_sqrts1():
    from sympy.functions.special.hyper import HyperRep_sqrts1
    assert _test_args(HyperRep_sqrts1(x, y))


def test_sympy__functions__special__hyper__HyperRep_sqrts2():
    from sympy.functions.special.hyper import HyperRep_sqrts2
    assert _test_args(HyperRep_sqrts2(x, y))


def test_sympy__functions__special__hyper__HyperRep_log2():
    from sympy.functions.special.hyper import HyperRep_log2
    assert _test_args(HyperRep_log2(x))


def test_sympy__functions__special__hyper__HyperRep_cosasin():
    from sympy.functions.special.hyper import HyperRep_cosasin
    assert _test_args(HyperRep_cosasin(x, y))


def test_sympy__functions__special__hyper__HyperRep_sinasin():
    from sympy.functions.special.hyper import HyperRep_sinasin
    assert _test_args(HyperRep_sinasin(x, y))


def test_sympy__functions__special__polynomials__OrthogonalPolynomial():
    pass


def test_sympy__functions__special__polynomials__jacobi():
    from sympy.functions.special.polynomials import jacobi
    assert _test_args(jacobi(x, 2, 2, 2))


def test_sympy__functions__special__polynomials__gegenbauer():
    from sympy.functions.special.polynomials import gegenbauer
    assert _test_args(gegenbauer(x, 2, 2))


def test_sympy__functions__special__polynomials__chebyshevt():
    from sympy.functions.special.polynomials import chebyshevt
    assert _test_args(chebyshevt(x, 2))


def test_sympy__functions__special__polynomials__chebyshevt_root():
    from sympy.functions.special.polynomials import chebyshevt_root
    assert _test_args(chebyshevt_root(3, 2))


def test_sympy__functions__special__polynomials__chebyshevu():
    from sympy.functions.special.polynomials import chebyshevu
    assert _test_args(chebyshevu(x, 2))


def test_sympy__functions__special__polynomials__chebyshevu_root():
    from sympy.functions.special.polynomials import chebyshevu_root
    assert _test_args(chebyshevu_root(3, 2))


def test_sympy__functions__special__polynomials__hermite():
    from sympy.functions.special.polynomials import hermite
    assert _test_args(hermite(x, 2))


def test_sympy__functions__special__polynomials__legendre():
    from sympy.functions.special.polynomials import legendre
    assert _test_args(legendre(x, 2))


def test_sympy__functions__special__polynomials__assoc_legendre():
    from sympy.functions.special.polynomials import assoc_legendre
    assert _test_args(assoc_legendre(x, 0, y))


def test_sympy__functions__special__polynomials__laguerre():
    from sympy.functions.special.polynomials import laguerre
    assert _test_args(laguerre(x, 2))


def test_sympy__functions__special__polynomials__assoc_laguerre():
    from sympy.functions.special.polynomials import assoc_laguerre
    assert _test_args(assoc_laguerre(x, 0, y))


def test_sympy__functions__special__spherical_harmonics__Ynm():
    from sympy.functions.special.spherical_harmonics import Ynm
    assert _test_args(Ynm(1, 1, x, y))


def test_sympy__functions__special__spherical_harmonics__Znm():
    from sympy.functions.special.spherical_harmonics import Znm
    assert _test_args(Znm(1, 1, x, y))


def test_sympy__functions__special__tensor_functions__LeviCivita():
    from sympy.functions.special.tensor_functions import LeviCivita
    assert _test_args(LeviCivita(x, y, 2))


def test_sympy__functions__special__tensor_functions__KroneckerDelta():
    from sympy.functions.special.tensor_functions import KroneckerDelta
    assert _test_args(KroneckerDelta(x, y))


def test_sympy__functions__special__zeta_functions__dirichlet_eta():
    from sympy.functions.special.zeta_functions import dirichlet_eta
    assert _test_args(dirichlet_eta(x))


def test_sympy__functions__special__zeta_functions__zeta():
    from sympy.functions.special.zeta_functions import zeta
    assert _test_args(zeta(101))


def test_sympy__functions__special__zeta_functions__lerchphi():
    from sympy.functions.special.zeta_functions import lerchphi
    assert _test_args(lerchphi(x, y, z))


def test_sympy__functions__special__zeta_functions__polylog():
    from sympy.functions.special.zeta_functions import polylog
    assert _test_args(polylog(x, y))


def test_sympy__integrals__integrals__Integral():
    from sympy.integrals.integrals import Integral
    assert _test_args(Integral(2, (x, 0, 1)))


def test_sympy__integrals__risch__NonElementaryIntegral():
    from sympy.integrals.risch import NonElementaryIntegral
    assert _test_args(NonElementaryIntegral(exp(-x**2), x))


def test_sympy__integrals__transforms__IntegralTransform():
    pass


def test_sympy__integrals__transforms__MellinTransform():
    from sympy.integrals.transforms import MellinTransform
    assert _test_args(MellinTransform(2, x, y))


def test_sympy__integrals__transforms__InverseMellinTransform():
    from sympy.integrals.transforms import InverseMellinTransform
    assert _test_args(InverseMellinTransform(2, x, y, 0, 1))


def test_sympy__integrals__transforms__LaplaceTransform():
    from sympy.integrals.transforms import LaplaceTransform
    assert _test_args(LaplaceTransform(2, x, y))


def test_sympy__integrals__transforms__InverseLaplaceTransform():
    from sympy.integrals.transforms import InverseLaplaceTransform
    assert _test_args(InverseLaplaceTransform(2, x, y, 0))


def test_sympy__integrals__transforms__FourierTypeTransform():
    pass


def test_sympy__integrals__transforms__InverseFourierTransform():
    from sympy.integrals.transforms import InverseFourierTransform
    assert _test_args(InverseFourierTransform(2, x, y))


def test_sympy__integrals__transforms__FourierTransform():
    from sympy.integrals.transforms import FourierTransform
    assert _test_args(FourierTransform(2, x, y))


def test_sympy__integrals__transforms__SineCosineTypeTransform():
    pass


def test_sympy__integrals__transforms__InverseSineTransform():
    from sympy.integrals.transforms import InverseSineTransform
    assert _test_args(InverseSineTransform(2, x, y))


def test_sympy__integrals__transforms__SineTransform():
    from sympy.integrals.transforms import SineTransform
    assert _test_args(SineTransform(2, x, y))


def test_sympy__integrals__transforms__InverseCosineTransform():
    from sympy.integrals.transforms import InverseCosineTransform
    assert _test_args(InverseCosineTransform(2, x, y))


def test_sympy__integrals__transforms__CosineTransform():
    from sympy.integrals.transforms import CosineTransform
    assert _test_args(CosineTransform(2, x, y))


def test_sympy__integrals__transforms__HankelTypeTransform():
    pass


def test_sympy__integrals__transforms__InverseHankelTransform():
    from sympy.integrals.transforms import InverseHankelTransform
    assert _test_args(InverseHankelTransform(2, x, y, 0))


def test_sympy__integrals__transforms__HankelTransform():
    from sympy.integrals.transforms import HankelTransform
    assert _test_args(HankelTransform(2, x, y, 0))


@pytest.mark.xfail
def test_sympy__liealgebras__cartan_type__CartanType_generator():
    from sympy.liealgebras.cartan_type import CartanType_generator
    assert _test_args(CartanType_generator("A2"))


@pytest.mark.xfail
def test_sympy__liealgebras__cartan_type__Standard_Cartan():
    from sympy.liealgebras.cartan_type import Standard_Cartan
    assert _test_args(Standard_Cartan("A", 2))


@pytest.mark.xfail
def test_sympy__liealgebras__weyl_group__WeylGroup():
    from sympy.liealgebras.weyl_group import WeylGroup
    assert _test_args(WeylGroup("B4"))


@pytest.mark.xfail
def test_sympy__liealgebras__root_system__RootSystem():
    from sympy.liealgebras.root_system import RootSystem
    assert _test_args(RootSystem("A2"))


@pytest.mark.xfail
def test_sympy__liealgebras__type_a__TypeA():
    from sympy.liealgebras.type_a import TypeA
    assert _test_args(TypeA(2))


@pytest.mark.xfail
def test_sympy__liealgebras__type_b__TypeB():
    from sympy.liealgebras.type_b import TypeB
    assert _test_args(TypeB(4))


@pytest.mark.xfail
def test_sympy__liealgebras__type_c__TypeC():
    from sympy.liealgebras.type_c import TypeC
    assert _test_args(TypeC(4))


@pytest.mark.xfail
def test_sympy__liealgebras__type_d__TypeD():
    from sympy.liealgebras.type_d import TypeD
    assert _test_args(TypeD(4))


@pytest.mark.xfail
def test_sympy__liealgebras__type_e__TypeE():
    from sympy.liealgebras.type_e import TypeE
    assert _test_args(TypeE(6))


@pytest.mark.xfail
def test_sympy__liealgebras__type_f__TypeF():
    from sympy.liealgebras.type_f import TypeF
    assert _test_args(TypeF(4))


@pytest.mark.xfail
def test_sympy__liealgebras__type_g__TypeG():
    from sympy.liealgebras.type_g import TypeG
    assert _test_args(TypeG(2))


def test_sympy__logic__boolalg__And():
    from sympy.logic.boolalg import And
    assert _test_args(And(x, y, 2))


def test_sympy__logic__boolalg__Boolean():
    pass


def test_sympy__logic__boolalg__BooleanFunction():
    from sympy.logic.boolalg import BooleanFunction
    assert _test_args(BooleanFunction(1, 2, 3))


def test_sympy__logic__boolalg__BooleanAtom():
    pass


def test_sympy__logic__boolalg__BooleanTrue():
    from sympy.logic.boolalg import true
    assert _test_args(true)


def test_sympy__logic__boolalg__BooleanFalse():
    from sympy.logic.boolalg import false
    assert _test_args(false)


def test_sympy__logic__boolalg__Equivalent():
    from sympy.logic.boolalg import Equivalent
    assert _test_args(Equivalent(x, 2))


def test_sympy__logic__boolalg__ITE():
    from sympy.logic.boolalg import ITE
    assert _test_args(ITE(x, y, 2))


def test_sympy__logic__boolalg__Implies():
    from sympy.logic.boolalg import Implies
    assert _test_args(Implies(x, y))


def test_sympy__logic__boolalg__Nand():
    from sympy.logic.boolalg import Nand
    assert _test_args(Nand(x, y, 2))


def test_sympy__logic__boolalg__Nor():
    from sympy.logic.boolalg import Nor
    assert _test_args(Nor(x, y))


def test_sympy__logic__boolalg__Not():
    from sympy.logic.boolalg import Not
    assert _test_args(Not(x))


def test_sympy__logic__boolalg__Or():
    from sympy.logic.boolalg import Or
    assert _test_args(Or(x, y))


def test_sympy__logic__boolalg__Xor():
    from sympy.logic.boolalg import Xor
    assert _test_args(Xor(x, y, 2))


def test_sympy__matrices__matrices__DeferredVector():
    from sympy.matrices.matrices import DeferredVector
    assert _test_args(DeferredVector("X"))


def test_sympy__matrices__expressions__matexpr__MatrixBase():
    pass


def test_sympy__matrices__immutable__ImmutableMatrix():
    from sympy.matrices.immutable import ImmutableMatrix
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


def test_sympy__matrices__immutable__ImmutableSparseMatrix():
    from sympy.matrices.immutable import ImmutableSparseMatrix
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


def test_sympy__matrices__expressions__slice__MatrixSlice():
    from sympy.matrices.expressions.slice import MatrixSlice
    from sympy.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', 4, 4)
    assert _test_args(MatrixSlice(X, (0, 2), (0, 2)))


def test_sympy__matrices__expressions__blockmatrix__BlockDiagMatrix():
    from sympy.matrices.expressions.blockmatrix import BlockDiagMatrix
    from sympy.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, x)
    Y = MatrixSymbol('Y', y, y)
    assert _test_args(BlockDiagMatrix(X, Y))


def test_sympy__matrices__expressions__blockmatrix__BlockMatrix():
    from sympy.matrices.expressions.blockmatrix import BlockMatrix
    from sympy.matrices.expressions import MatrixSymbol, ZeroMatrix
    X = MatrixSymbol('X', x, x)
    Y = MatrixSymbol('Y', y, y)
    Z = MatrixSymbol('Z', x, y)
    O = ZeroMatrix(y, x)
    assert _test_args(BlockMatrix([[X, Z], [O, Y]]))


def test_sympy__matrices__expressions__inverse__Inverse():
    from sympy.matrices.expressions.inverse import Inverse
    from sympy.matrices.expressions import MatrixSymbol
    assert _test_args(Inverse(MatrixSymbol('A', 3, 3)))


def test_sympy__matrices__expressions__matadd__MatAdd():
    from sympy.matrices.expressions.matadd import MatAdd
    from sympy.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', x, y)
    assert _test_args(MatAdd(X, Y))


def test_sympy__matrices__expressions__matexpr__Identity():
    from sympy.matrices.expressions.matexpr import Identity
    assert _test_args(Identity(3))


def test_sympy__matrices__expressions__matexpr__MatrixExpr():
    pass


def test_sympy__matrices__expressions__matexpr__MatrixElement():
    from sympy.matrices.expressions.matexpr import MatrixSymbol, MatrixElement
    from sympy import S
    assert _test_args(MatrixElement(MatrixSymbol('A', 3, 5), Integer(2), Integer(3)))


@pytest.mark.xfail
def test_sympy__matrices__expressions__matexpr__MatrixSymbol():
    from sympy.matrices.expressions.matexpr import MatrixSymbol
    assert _test_args(MatrixSymbol('A', 3, 5))


def test_sympy__matrices__expressions__matexpr__ZeroMatrix():
    from sympy.matrices.expressions.matexpr import ZeroMatrix
    assert _test_args(ZeroMatrix(3, 5))


def test_sympy__matrices__expressions__matmul__MatMul():
    from sympy.matrices.expressions.matmul import MatMul
    from sympy.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', y, x)
    assert _test_args(MatMul(X, Y))


def test_sympy__matrices__expressions__diagonal__DiagonalMatrix():
    from sympy.matrices.expressions.diagonal import DiagonalMatrix
    from sympy.matrices.expressions import MatrixSymbol
    x = MatrixSymbol('x', 10, 1)
    assert _test_args(DiagonalMatrix(x))


def test_sympy__matrices__expressions__diagonal__DiagonalOf():
    from sympy.matrices.expressions.diagonal import DiagonalOf
    from sympy.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('x', 10, 10)
    assert _test_args(DiagonalOf(X))


def test_sympy__matrices__expressions__hadamard__HadamardProduct():
    from sympy.matrices.expressions.hadamard import HadamardProduct
    from sympy.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, y)
    Y = MatrixSymbol('Y', x, y)
    assert _test_args(HadamardProduct(X, Y))


def test_sympy__matrices__expressions__matpow__MatPow():
    from sympy.matrices.expressions.matpow import MatPow
    from sympy.matrices.expressions import MatrixSymbol
    X = MatrixSymbol('X', x, x)
    assert _test_args(MatPow(X, 2))


def test_sympy__matrices__expressions__transpose__Transpose():
    from sympy.matrices.expressions.transpose import Transpose
    from sympy.matrices.expressions import MatrixSymbol
    assert _test_args(Transpose(MatrixSymbol('A', 3, 5)))


def test_sympy__matrices__expressions__adjoint__Adjoint():
    from sympy.matrices.expressions.adjoint import Adjoint
    from sympy.matrices.expressions import MatrixSymbol
    assert _test_args(Adjoint(MatrixSymbol('A', 3, 5)))


def test_sympy__matrices__expressions__trace__Trace():
    from sympy.matrices.expressions.trace import Trace
    from sympy.matrices.expressions import MatrixSymbol
    assert _test_args(Trace(MatrixSymbol('A', 3, 3)))


def test_sympy__matrices__expressions__determinant__Determinant():
    from sympy.matrices.expressions.determinant import Determinant
    from sympy.matrices.expressions import MatrixSymbol
    assert _test_args(Determinant(MatrixSymbol('A', 3, 3)))


def test_sympy__matrices__expressions__funcmatrix__FunctionMatrix():
    from sympy.matrices.expressions.funcmatrix import FunctionMatrix
    from sympy import Lambda, symbols
    i, j = symbols('i,j')
    assert _test_args(FunctionMatrix(3, 3, Lambda((i, j), i - j) ))


def test_sympy__matrices__expressions__fourier__DFT():
    from sympy.matrices.expressions.fourier import DFT
    from sympy import S
    assert _test_args(DFT(Integer(2)))


def test_sympy__matrices__expressions__fourier__IDFT():
    from sympy.matrices.expressions.fourier import IDFT
    from sympy import S
    assert _test_args(IDFT(Integer(2)))

from sympy.matrices.expressions import MatrixSymbol
X = MatrixSymbol('X', 10, 10)


def test_sympy__matrices__expressions__factorizations__LofLU():
    from sympy.matrices.expressions.factorizations import LofLU
    assert _test_args(LofLU(X))


def test_sympy__matrices__expressions__factorizations__UofLU():
    from sympy.matrices.expressions.factorizations import UofLU
    assert _test_args(UofLU(X))


def test_sympy__matrices__expressions__factorizations__QofQR():
    from sympy.matrices.expressions.factorizations import QofQR
    assert _test_args(QofQR(X))


def test_sympy__matrices__expressions__factorizations__RofQR():
    from sympy.matrices.expressions.factorizations import RofQR
    assert _test_args(RofQR(X))


def test_sympy__matrices__expressions__factorizations__LofCholesky():
    from sympy.matrices.expressions.factorizations import LofCholesky
    assert _test_args(LofCholesky(X))


def test_sympy__matrices__expressions__factorizations__UofCholesky():
    from sympy.matrices.expressions.factorizations import UofCholesky
    assert _test_args(UofCholesky(X))


def test_sympy__matrices__expressions__factorizations__EigenVectors():
    from sympy.matrices.expressions.factorizations import EigenVectors
    assert _test_args(EigenVectors(X))


def test_sympy__matrices__expressions__factorizations__EigenValues():
    from sympy.matrices.expressions.factorizations import EigenValues
    assert _test_args(EigenValues(X))


def test_sympy__matrices__expressions__factorizations__UofSVD():
    from sympy.matrices.expressions.factorizations import UofSVD
    assert _test_args(UofSVD(X))


def test_sympy__matrices__expressions__factorizations__VofSVD():
    from sympy.matrices.expressions.factorizations import VofSVD
    assert _test_args(VofSVD(X))


def test_sympy__matrices__expressions__factorizations__SofSVD():
    from sympy.matrices.expressions.factorizations import SofSVD
    assert _test_args(SofSVD(X))


def test_sympy__matrices__expressions__factorizations__Factorization():
    pass


def test_sympy__core__numbers__AlgebraicNumber():
    from sympy.core.numbers import AlgebraicNumber
    assert _test_args(AlgebraicNumber(sqrt(2), [1, 2, 3]))


def test_sympy__polys__polytools__GroebnerBasis():
    from sympy.polys.polytools import GroebnerBasis
    assert _test_args(GroebnerBasis([x, y, z], x, y, z))


def test_sympy__polys__polytools__Poly():
    from sympy.polys.polytools import Poly
    assert _test_args(Poly(2, x, y))


def test_sympy__polys__polytools__PurePoly():
    from sympy.polys.polytools import PurePoly
    assert _test_args(PurePoly(2, x, y))


def test_sympy__polys__rootoftools__RootOf():
    from sympy.polys.rootoftools import RootOf
    assert _test_args(RootOf(x**3 + x + 1, 0))


def test_sympy__polys__rootoftools__RootSum():
    from sympy.polys.rootoftools import RootSum
    assert _test_args(RootSum(x**3 + x + 1, sin))


def test_sympy__series__limits__Limit():
    from sympy.series.limits import Limit
    assert _test_args(Limit(x, x, 0, dir='-'))


def test_sympy__series__order__Order():
    from sympy.series.order import Order
    assert _test_args(Order(1, x, y))


def test_sympy__simplify__hyperexpand__Hyper_Function():
    from sympy.simplify.hyperexpand import Hyper_Function
    assert _test_args(Hyper_Function([2], [1]))


def test_sympy__simplify__hyperexpand__G_Function():
    from sympy.simplify.hyperexpand import G_Function
    assert _test_args(G_Function([2], [1], [], []))


def test_sympy__tensor__indexed__Idx():
    from sympy.tensor.indexed import Idx
    assert _test_args(Idx('test'))
    assert _test_args(Idx(1, (0, 10)))


def test_sympy__tensor__indexed__Indexed():
    from sympy.tensor.indexed import Indexed, Idx
    assert _test_args(Indexed('A', Idx('i'), Idx('j')))


def test_sympy__tensor__indexed__IndexedBase():
    from sympy.tensor.indexed import IndexedBase
    assert _test_args(IndexedBase('A', shape=(x, y)))
    assert _test_args(IndexedBase('A', 1))
    assert _test_args(IndexedBase('A')[0, 1])


def test_sympy__tensor__tensor__TensorIndexType():
    from sympy.tensor.tensor import TensorIndexType
    assert _test_args(TensorIndexType('Lorentz', metric=False))


def test_sympy__tensor__tensor__TensorSymmetry():
    from sympy.tensor.tensor import TensorSymmetry, get_symmetric_group_sgs
    assert _test_args(TensorSymmetry(get_symmetric_group_sgs(2)))


def test_sympy__tensor__tensor__TensorType():
    from sympy.tensor.tensor import TensorIndexType, TensorSymmetry, get_symmetric_group_sgs, TensorType
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    assert _test_args(TensorType([Lorentz], sym))


def test_sympy__tensor__tensor__TensorHead():
    from sympy.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, TensorHead
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    assert _test_args(TensorHead('p', S1, 0))


def test_sympy__tensor__tensor__TensorIndex():
    from sympy.tensor.tensor import TensorIndexType, TensorIndex
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    assert _test_args(TensorIndex('i', Lorentz))


def test_sympy__tensor__tensor__TensExpr():
    pass


def test_sympy__tensor__tensor__TensAdd():
    from sympy.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, tensor_indices, TensAdd
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p, q = S1('p,q')
    t1 = p(a)
    t2 = q(a)
    assert _test_args(TensAdd(t1, t2))


def test_sympy__tensor__tensor__Tensor():
    from sympy.core import S
    from sympy.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, tensor_indices, TensMul, TIDS
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p = S1('p')
    assert _test_args(p(a))


def test_sympy__tensor__tensor__TensMul():
    from sympy.core import S
    from sympy.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, tensor_indices, TensMul, TIDS
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p = S1('p')
    q = S1('q')
    assert _test_args(3*p(a)*q(b))


def test_as_coeff_add():
    assert (7, (3*x, 4*x**2)) == (7 + 3*x + 4*x**2).as_coeff_add()


def test_sympy__geometry__curve__Curve():
    from sympy.geometry.curve import Curve
    assert _test_args(Curve((x, 1), (x, 0, 1)))


def test_sympy__geometry__point__Point():
    from sympy.geometry.point import Point
    assert _test_args(Point(0, 1))


def test_sympy__geometry__point__Point2D():
    from sympy.geometry.point import Point2D
    assert _test_args(Point2D(0, 1))


def test_sympy__geometry__point__Point3D():
    from sympy.geometry.point import Point3D
    assert _test_args(Point3D(0, 1, 2))


def test_sympy__geometry__ellipse__Ellipse():
    from sympy.geometry.ellipse import Ellipse
    assert _test_args(Ellipse((0, 1), 2, 3))


def test_sympy__geometry__ellipse__Circle():
    from sympy.geometry.ellipse import Circle
    assert _test_args(Circle((0, 1), 2))


def test_sympy__geometry__line__LinearEntity():
    pass


def test_sympy__geometry__line__Line():
    from sympy.geometry.line import Line
    assert _test_args(Line((0, 1), (2, 3)))


def test_sympy__geometry__line__Ray():
    from sympy.geometry.line import Ray
    assert _test_args(Ray((0, 1), (2, 3)))


def test_sympy__geometry__line__Segment():
    from sympy.geometry.line import Segment
    assert _test_args(Segment((0, 1), (2, 3)))


def test_sympy__geometry__line3d__LinearEntity3D():
    pass


def test_sympy__geometry__line3d__Line3D():
    from sympy.geometry.line3d import Line3D
    assert _test_args(Line3D((0, 1, 1), (2, 3, 4)))


def test_sympy__geometry__line3d__Segment3D():
    from sympy.geometry.line3d import Segment3D
    assert _test_args(Segment3D((0, 1, 1), (2, 3, 4)))


def test_sympy__geometry__line3d__Ray3D():
    from sympy.geometry.line3d import Ray3D
    assert _test_args(Ray3D((0, 1, 1), (2, 3, 4)))


def test_sympy__geometry__plane__Plane():
    from sympy.geometry.plane import Plane
    assert _test_args(Plane((1, 1, 1), (-3, 4, -2), (1, 2, 3)))


def test_sympy__geometry__polygon__Polygon():
    from sympy.geometry.polygon import Polygon
    assert _test_args(Polygon((0, 1), (2, 3), (4, 5), (6, 7)))


def test_sympy__geometry__polygon__RegularPolygon():
    from sympy.geometry.polygon import RegularPolygon
    assert _test_args(RegularPolygon((0, 1), 2, 3, 4))


def test_sympy__geometry__polygon__Triangle():
    from sympy.geometry.polygon import Triangle
    assert _test_args(Triangle((0, 1), (2, 3), (4, 5)))


def test_sympy__geometry__entity__GeometryEntity():
    from sympy.geometry.entity import GeometryEntity
    from sympy.geometry.point import Point
    assert _test_args(GeometryEntity(Point(1, 0), 1, [1, 2]))


def test_sympy__geometry__entity__GeometrySet():
    pass


def test_sympy__diffgeom__diffgeom__Manifold():
    from sympy.diffgeom import Manifold
    assert _test_args(Manifold('name', 3))


def test_sympy__diffgeom__diffgeom__Patch():
    from sympy.diffgeom import Manifold, Patch
    assert _test_args(Patch('name', Manifold('name', 3)))


def test_sympy__diffgeom__diffgeom__CoordSystem():
    from sympy.diffgeom import Manifold, Patch, CoordSystem
    assert _test_args(CoordSystem('name', Patch('name', Manifold('name', 3))))


@pytest.mark.xfail
def test_sympy__diffgeom__diffgeom__Point():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, Point
    assert _test_args(Point(
        CoordSystem('name', Patch('name', Manifold('name', 3))), [x, y]))


def test_sympy__diffgeom__diffgeom__BaseScalarField():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseScalarField(cs, 0))


def test_sympy__diffgeom__diffgeom__BaseVectorField():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseVectorField
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseVectorField(cs, 0))


def test_sympy__diffgeom__diffgeom__Differential():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(Differential(BaseScalarField(cs, 0)))


def test_sympy__diffgeom__diffgeom__Commutator():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseVectorField, Commutator
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    cs1 = CoordSystem('name1', Patch('name', Manifold('name', 3)))
    v = BaseVectorField(cs, 0)
    v1 = BaseVectorField(cs1, 0)
    assert _test_args(Commutator(v, v1))


def test_sympy__diffgeom__diffgeom__TensorProduct():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential, TensorProduct
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    assert _test_args(TensorProduct(d, d))


def test_sympy__diffgeom__diffgeom__WedgeProduct():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential, WedgeProduct
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    d1 = Differential(BaseScalarField(cs, 1))
    assert _test_args(WedgeProduct(d, d1))


def test_sympy__diffgeom__diffgeom__LieDerivative():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseScalarField, Differential, BaseVectorField, LieDerivative
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    d = Differential(BaseScalarField(cs, 0))
    v = BaseVectorField(cs, 0)
    assert _test_args(LieDerivative(v, d))


@pytest.mark.xfail
def test_sympy__diffgeom__diffgeom__BaseCovarDerivativeOp():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseCovarDerivativeOp
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    assert _test_args(BaseCovarDerivativeOp(cs, 0, [[[0, ]*3, ]*3, ]*3))


def test_sympy__diffgeom__diffgeom__CovarDerivativeOp():
    from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseVectorField, CovarDerivativeOp
    cs = CoordSystem('name', Patch('name', Manifold('name', 3)))
    v = BaseVectorField(cs, 0)
    _test_args(CovarDerivativeOp(v, [[[0, ]*3, ]*3, ]*3))


def test_sympy__categories__baseclasses__Class():
    from sympy.categories.baseclasses import Class
    assert _test_args(Class())


def test_sympy__categories__baseclasses__Object():
    from sympy.categories import Object
    assert _test_args(Object("A"))


@pytest.mark.xfail
def test_sympy__categories__baseclasses__Morphism():
    from sympy.categories import Object, Morphism
    assert _test_args(Morphism(Object("A"), Object("B")))


def test_sympy__categories__baseclasses__IdentityMorphism():
    from sympy.categories import Object, IdentityMorphism
    assert _test_args(IdentityMorphism(Object("A")))


def test_sympy__categories__baseclasses__NamedMorphism():
    from sympy.categories import Object, NamedMorphism
    assert _test_args(NamedMorphism(Object("A"), Object("B"), "f"))


def test_sympy__categories__baseclasses__CompositeMorphism():
    from sympy.categories import Object, NamedMorphism, CompositeMorphism
    A = Object("A")
    B = Object("B")
    C = Object("C")
    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    assert _test_args(CompositeMorphism(f, g))


def test_sympy__categories__baseclasses__Diagram():
    from sympy.categories import Object, NamedMorphism, Diagram
    A = Object("A")
    B = Object("B")
    C = Object("C")
    f = NamedMorphism(A, B, "f")
    d = Diagram([f])
    assert _test_args(d)


def test_sympy__categories__baseclasses__Category():
    from sympy.categories import Object, NamedMorphism, Diagram, Category
    A = Object("A")
    B = Object("B")
    C = Object("C")
    f = NamedMorphism(A, B, "f")
    g = NamedMorphism(B, C, "g")
    d1 = Diagram([f, g])
    d2 = Diagram([f])
    K = Category("K", commutative_diagrams=[d1, d2])
    assert _test_args(K)


def test_sympy__ntheory__factor___totient():
    from sympy.ntheory.factor_ import totient
    k = symbols('k', integer=True)
    t = totient(k)
    assert _test_args(t)


def test_sympy__ntheory__factor___divisor_sigma():
    from sympy.ntheory.factor_ import divisor_sigma
    k = symbols('k', integer=True)
    n = symbols('n', integer=True)
    t = divisor_sigma(n, k)
    assert _test_args(t)


def test_sympy__ntheory__residue_ntheory__mobius():
    from sympy.ntheory import mobius
    assert _test_args(mobius(2))


def test_sympy__printing__codeprinter__Assignment():
    from sympy.printing.codeprinter import Assignment
    assert _test_args(Assignment(x, y))


def test_sympy__vector__coordsysrect__CoordSysCartesian():
    from sympy.vector.coordsysrect import CoordSysCartesian
    assert _test_args(CoordSysCartesian('C'))


def test_sympy__vector__point__Point():
    from sympy.vector.point import Point
    assert _test_args(Point('P'))


def test_sympy__vector__basisdependent__BasisDependent():
    from sympy.vector.basisdependent import BasisDependent
    # These classes have been created to maintain an OOP hierarchy
    # for Vectors and Dyadics. Are NOT meant to be initialized


def test_sympy__vector__basisdependent__BasisDependentMul():
    from sympy.vector.basisdependent import BasisDependentMul
    # These classes have been created to maintain an OOP hierarchy
    # for Vectors and Dyadics. Are NOT meant to be initialized


def test_sympy__vector__basisdependent__BasisDependentAdd():
    from sympy.vector.basisdependent import BasisDependentAdd
    # These classes have been created to maintain an OOP hierarchy
    # for Vectors and Dyadics. Are NOT meant to be initialized


def test_sympy__vector__basisdependent__BasisDependentZero():
    from sympy.vector.basisdependent import BasisDependentZero
    # These classes have been created to maintain an OOP hierarchy
    # for Vectors and Dyadics. Are NOT meant to be initialized


def test_sympy__vector__vector__BaseVector():
    from sympy.vector.vector import BaseVector
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(BaseVector('Ci', 0, C, ' ', ' '))


def test_sympy__vector__vector__VectorAdd():
    from sympy.vector.vector import VectorAdd, VectorMul
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    from sympy.abc import a, b, c, x, y, z
    v1 = a*C.i + b*C.j + c*C.k
    v2 = x*C.i + y*C.j + z*C.k
    assert _test_args(VectorAdd(v1, v2))
    assert _test_args(VectorMul(x, v1))


def test_sympy__vector__vector__VectorMul():
    from sympy.vector.vector import VectorMul
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    from sympy.abc import a
    assert _test_args(VectorMul(a, C.i))


def test_sympy__vector__vector__VectorZero():
    from sympy.vector.vector import VectorZero
    assert _test_args(VectorZero())


def test_sympy__vector__vector__Vector():
    from sympy.vector.vector import Vector
    # Vector is never to be initialized using args
    pass


def test_sympy__vector__dyadic__Dyadic():
    from sympy.vector.dyadic import Dyadic
    # Dyadic is never to be initialized using args
    pass


def test_sympy__vector__dyadic__BaseDyadic():
    from sympy.vector.dyadic import BaseDyadic
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(BaseDyadic(C.i, C.j))


def test_sympy__vector__dyadic__DyadicMul():
    from sympy.vector.dyadic import BaseDyadic, DyadicMul
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(DyadicMul(3, BaseDyadic(C.i, C.j)))


def test_sympy__vector__dyadic__DyadicAdd():
    from sympy.vector.dyadic import BaseDyadic, DyadicAdd
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(2 * DyadicAdd(BaseDyadic(C.i, C.i),
                                    BaseDyadic(C.i, C.j)))


def test_sympy__vector__dyadic__DyadicZero():
    from sympy.vector.dyadic import DyadicZero
    assert _test_args(DyadicZero())


def test_sympy__vector__deloperator__Del():
    from sympy.vector.deloperator import Del
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(Del(C))


def test_sympy__vector__orienters__Orienter():
    from sympy.vector.orienters import Orienter
    # Not to be initialized


def test_sympy__vector__orienters__ThreeAngleOrienter():
    from sympy.vector.orienters import ThreeAngleOrienter
    # Not to be initialized


def test_sympy__vector__orienters__AxisOrienter():
    from sympy.vector.orienters import AxisOrienter
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(AxisOrienter(x, C.i))


def test_sympy__vector__orienters__BodyOrienter():
    from sympy.vector.orienters import BodyOrienter
    assert _test_args(BodyOrienter(x, y, z, '123'))


def test_sympy__vector__orienters__SpaceOrienter():
    from sympy.vector.orienters import SpaceOrienter
    assert _test_args(SpaceOrienter(x, y, z, '123'))


def test_sympy__vector__orienters__QuaternionOrienter():
    from sympy.vector.orienters import QuaternionOrienter
    a, b, c, d = symbols('a b c d')
    assert _test_args(QuaternionOrienter(a, b, c, d))


def test_sympy__vector__scalar__BaseScalar():
    from sympy.vector.scalar import BaseScalar
    from sympy.vector.coordsysrect import CoordSysCartesian
    C = CoordSysCartesian('C')
    assert _test_args(BaseScalar('Cx', 0, C, ' ', ' '))
