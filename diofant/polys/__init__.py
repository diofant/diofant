"""Polynomial manipulation algorithms and algebraic objects."""

from .polytools import (Poly, PurePoly,
                        parallel_poly_from_expr, degree, LC,
                        LM, LT, div, rem, quo,
                        exquo, half_gcdex, gcdex, invert, subresultants,
                        resultant, discriminant, cofactors, gcd,
                        lcm, terms_gcd, trunc, monic, content,
                        primitive, compose, decompose,
                        sqf_norm, sqf_part, sqf_list, sqf, factor_list,
                        factor, count_roots, real_roots, nroots,
                        cancel, reduced, groebner, GroebnerBasis)
from .polyfuncs import symmetrize, horner, interpolate, viete
from .rationaltools import together
from .polyerrors import (BasePolynomialError,
                         ExactQuotientFailedError, PolynomialDivisionFailedError,
                         OperationNotSupportedError, HeuristicGCDFailedError,
                         HomomorphismFailedError, IsomorphismFailedError,
                         ExtraneousFactorsError, EvaluationFailedError,
                         RefinementFailedError, CoercionFailedError, NotInvertibleError,
                         NotReversibleError, NotAlgebraicError, DomainError,
                         PolynomialError, UnificationFailedError,
                         GeneratorsError, GeneratorsNeededError, ComputationFailedError,
                         UnivariatePolynomialError, MultivariatePolynomialError,
                         PolificationFailedError, OptionError, FlagError)
from .numberfields import (minimal_polynomial,
                           primitive_element, field_isomorphism)
from .monomials import Monomial, itermonomials
from .orderings import (lex, grlex, grevlex, ilex,
                        igrlex, igrevlex)
from .rootoftools import RootOf, RootSum
from .polyroots import roots
from .constructor import construct_domain
from .specialpolys import (swinnerton_dyer_poly,
                           cyclotomic_poly, symmetric_poly,
                           random_poly, interpolating_poly)
from .orthopolys import (jacobi_poly, chebyshevt_poly,
                         chebyshevu_poly, hermite_poly, legendre_poly,
                         laguerre_poly, spherical_bessel_fn)
from .partfrac import apart, apart_list, assemble_partfrac_list
from .polyoptions import Options
from .rings import PolynomialRing, ring
from .fields import FractionField, field
from .univar import UnivarPolynomialRing


__all__ = ('Poly', 'PurePoly',
           'parallel_poly_from_expr', 'degree', 'LC',
           'LM', 'LT', 'div', 'rem', 'quo', 'exquo', 'half_gcdex',
           'gcdex', 'invert', 'subresultants', 'resultant', 'discriminant',
           'cofactors', 'gcd', 'lcm', 'terms_gcd',
           'trunc', 'monic', 'content', 'primitive', 'compose', 'decompose',
           'sqf_norm', 'sqf_part', 'sqf_list', 'sqf', 'factor_list',
           'factor', 'count_roots', 'real_roots', 'nroots', 'cancel',
           'reduced', 'groebner', 'GroebnerBasis', 'symmetrize',
           'horner', 'interpolate', 'viete', 'together', 'BasePolynomialError',
           'ExactQuotientFailedError', 'PolynomialDivisionFailedError',
           'OperationNotSupportedError', 'HeuristicGCDFailedError',
           'HomomorphismFailedError', 'IsomorphismFailedError', 'ExtraneousFactorsError',
           'EvaluationFailedError', 'RefinementFailedError', 'CoercionFailedError',
           'NotInvertibleError', 'NotReversibleError', 'NotAlgebraicError', 'DomainError',
           'PolynomialError', 'UnificationFailedError', 'GeneratorsError',
           'GeneratorsNeededError', 'ComputationFailedError', 'UnivariatePolynomialError',
           'MultivariatePolynomialError', 'PolificationFailedError', 'OptionError',
           'FlagError', 'minimal_polynomial', 'primitive_element',
           'field_isomorphism', 'Monomial', 'itermonomials', 'lex', 'grlex',
           'grevlex', 'ilex', 'igrlex', 'igrevlex', 'RootOf', 'RootSum',
           'roots', 'construct_domain', 'swinnerton_dyer_poly',
           'cyclotomic_poly', 'symmetric_poly', 'random_poly',
           'interpolating_poly', 'jacobi_poly', 'chebyshevt_poly',
           'chebyshevu_poly', 'hermite_poly', 'legendre_poly',
           'laguerre_poly', 'spherical_bessel_fn', 'apart', 'apart_list',
           'assemble_partfrac_list', 'Options', 'PolynomialRing', 'ring',
           'FractionField', 'field', 'UnivarPolynomialRing')
