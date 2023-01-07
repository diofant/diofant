"""Test whether all elements of cls.args are instances of Basic."""

# NOTE: keep tests sorted by (module, class name) key.

import inspect
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
                     Xor, ZeroMatrix, divisor_sigma, false, mobius, sin,
                     symbols, totient, true)
from diofant.abc import a, b, w, x, y, z
from diofant.concrete.expr_with_intlimits import ExprWithIntLimits
from diofant.concrete.expr_with_limits import AddWithLimits, ExprWithLimits
from diofant.core.function import Application, AppliedUndef
from diofant.core.numbers import (Catalan, ComplexInfinity, EulerGamma, Exp1,
                                  GoldenRatio, Half, ImaginaryUnit, Infinity,
                                  NaN, NegativeInfinity, NegativeOne,
                                  NumberSymbol, One, Pi, Zero)
from diofant.core.trace import Tr
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
from diofant.tensor import ImmutableDenseNDimArray, ImmutableSparseNDimArray
from diofant.tensor.tensor import (TensAdd, TensorHead, TensorIndex,
                                   TensorIndexType, TensorSymmetry, TensorType,
                                   get_symmetric_group_sgs, tensor_indices)
from diofant.utilities.exceptions import DiofantDeprecationWarning


__all__ = ()


def test_all_classes_are_tested():
    this = os.path.split(__file__)[0]
    path = os.path.join(this, os.pardir, os.pardir)
    diofant_path = os.path.abspath(path)
    prefix = os.path.split(diofant_path)[0] + os.sep

    re_cls = re.compile(r'^class ([A-Za-z][A-Za-z0-9_]*)\s*\(', re.MULTILINE)

    modules = {}

    for root, _, files in os.walk(diofant_path):
        module = root.replace(prefix, '').replace(os.sep, '.')

        for file in files:
            if file.startswith(('_', 'test_', 'bench_')):
                continue
            if not file.endswith('.py'):
                continue

            with open(os.path.join(root, file), encoding='utf-8') as f:
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
    assert _test_args(FunctionMatrix(3, 3, Lambda((i, j), i - j)))


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


def test_diofant__calculus__limits__Limit():
    assert _test_args(Limit(x, x, 0, dir=1))


def test_diofant__calculus__order__Order():
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
    a = tensor_indices('a', Lorentz)
    sym = TensorSymmetry(get_symmetric_group_sgs(1))
    S1 = TensorType([Lorentz], sym)
    p, q = S1('p,q')
    t1 = p(a)
    t2 = q(a)
    assert _test_args(TensAdd(t1, t2))


def test_diofant__tensor__tensor__Tensor():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a = tensor_indices('a', Lorentz)
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
