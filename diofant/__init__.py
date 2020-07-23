"""Diofant is a Python library for symbolic mathematics."""

import os
DIOFANT_DEBUG = os.getenv('DIOFANT_DEBUG', 'False') != 'False'
del os

import pkg_resources
__version__ = pkg_resources.get_distribution(__name__).version
del pkg_resources

from .core import (Add, Atom, AtomicExpr, Basic, Catalan, Derivative, Dict,
                   Dummy, E, Eq, Equality, EulerGamma, Expr, Float, Function,
                   FunctionClass, Ge, GoldenRatio, GreaterThan, Gt, I, Integer,
                   Lambda, Le, LessThan, Lt, Mod, Mul, N, Ne, Number,
                   NumberSymbol, PoleError, Pow, PrecisionExhausted, Rational,
                   Rel, Relational, S, StrictGreaterThan, StrictLessThan, Subs,
                   Symbol, SympifyError, Tuple, Unequality, Wild, WildFunction,
                   cacheit, comp, count_ops, diff, evaluate, expand,
                   expand_complex, expand_func, expand_log, expand_mul,
                   expand_multinomial, expand_power_base, expand_power_exp,
                   expand_trig, factor_nc, factor_terms, gcd_terms, igcd, ilcm,
                   integer_digits, integer_nthroot, mod_inverse, nan, nfloat,
                   oo, pi, preorder_traversal, prod, symbols, sympify, var,
                   vectorize, zoo)
from .logic import (ITE, And, Equivalent, Implies, Nand, Nor, Not, Or, POSform,
                    SOPform, Xor, bool_map, false, satisfiable, simplify_logic,
                    to_cnf, to_dnf, to_nnf, true)
from .polys import (LC, LM, LT, BasePolynomialError, CoercionFailed,
                    ComputationFailed, DomainError, EvaluationFailed,
                    ExactQuotientFailed, ExtraneousFactors, FlagError,
                    FractionField, GeneratorsError, GeneratorsNeeded,
                    GroebnerBasis, HeuristicGCDFailed, HomomorphismFailed,
                    IsomorphismFailed, Monomial, MultivariatePolynomialError,
                    NotAlgebraic, NotInvertible, NotReversible,
                    OperationNotSupported, OptionError, Options,
                    PolificationFailed, Poly, PolynomialDivisionFailed,
                    PolynomialError, PolynomialRing, PurePoly,
                    RefinementFailed, RootOf, RootSum, UnificationFailed,
                    UnivariatePolynomialError, UnivarPolynomialRing, apart,
                    apart_list, assemble_partfrac_list, cancel,
                    chebyshevt_poly, chebyshevu_poly, cofactors, compose,
                    construct_domain, content, count_roots, cyclotomic_poly,
                    decompose, degree, degree_list, discriminant, div, exquo,
                    factor, factor_list, field, field_isomorphism, gcd,
                    gcd_list, gcdex, grevlex, grlex, groebner, half_gcdex,
                    hermite_poly, horner, igrevlex, igrlex, ilex, interpolate,
                    interpolating_poly, invert, itermonomials, jacobi_poly,
                    laguerre_poly, lcm, lcm_list, legendre_poly, lex,
                    minimal_polynomial, monic, nroots, parallel_poly_from_expr,
                    poly, prem, primitive, primitive_element, quo, random_poly,
                    real_roots, reduced, rem, resultant, ring, roots,
                    spherical_bessel_fn, sqf, sqf_list, sqf_norm, sqf_part,
                    subresultants, swinnerton_dyer_poly, symmetric_poly,
                    symmetrize, terms_gcd, together, trunc, viete)
from .domains import (CC, EX, FF, GF, GROUND_TYPES, QQ, RR, ZZ, AlgebraicField,
                      ComplexAlgebraicField, ComplexField, Domain,
                      ExpressionDomain, FF_gmpy, FF_python, FiniteField,
                      IntegerRing, PythonRational, QQ_gmpy, QQ_python,
                      RationalField, RealAlgebraicField, RealField, ZZ_gmpy,
                      ZZ_python)
from .series import *
from .functions import *
from .ntheory import *
from .concrete import *
from .simplify import *
from .sets import *
from .solvers import *
from .matrices import *
from .geometry import *
from .utilities import *
from .integrals import *
from .tensor import *
from .calculus import *
from .combinatorics import *
from .plotting import *
from .printing import *
from .interactive import *
