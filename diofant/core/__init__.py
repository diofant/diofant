"""Core module. Provides the basic operations needed in diofant.
"""

from .sympify import sympify, SympifyError
from .cache import cacheit
from .basic import Basic, Atom, preorder_traversal
from .singleton import S
from .expr import Expr, AtomicExpr
from .symbol import Symbol, Wild, Dummy, symbols, var
from .numbers import (Number, Float, Rational, Integer,
                      NumberSymbol, E, I, nan, oo, pi, zoo, comp,
                      mod_inverse, integer_digits)
from .power import Pow, integer_nthroot
from .mul import Mul
from .add import Add
from .mod import Mod
from .relational import (Rel, Eq, Ne, Lt, Le, Gt, Ge, Equality, Relational,
                         GreaterThan, LessThan, Unequality, StrictGreaterThan,
                         StrictLessThan)
from .multidimensional import vectorize
from .function import (Lambda, WildFunction, Derivative, diff,
                       FunctionClass, Function, Subs, expand, PoleError,
                       count_ops, expand_mul, expand_log, expand_func,
                       expand_trig, expand_complex, expand_multinomial, nfloat,
                       expand_power_base, expand_power_exp)
from .evalf import PrecisionExhausted, N
from .containers import Tuple, Dict
from .exprtools import gcd_terms, factor_terms, factor_nc
from .evaluate import evaluate

# expose singletons
Catalan = S.Catalan
EulerGamma = S.EulerGamma
GoldenRatio = S.GoldenRatio


__all__ = ('sympify', 'SympifyError', 'cacheit', 'Basic', 'Atom',
           'preorder_traversal', 'S', 'Expr', 'AtomicExpr', 'Symbol',
           'Wild', 'Dummy', 'symbols', 'var', 'Number', 'Float', 'Rational',
           'Integer', 'NumberSymbol', 'E', 'I', 'nan', 'oo',
           'pi', 'zoo', 'comp', 'mod_inverse', 'integer_digits', 'Pow',
           'integer_nthroot', 'Mul', 'Add', 'Mod', 'Rel', 'Eq', 'Ne',
           'Lt', 'Le', 'Gt', 'Ge', 'Equality', 'GreaterThan', 'LessThan',
           'Unequality', 'StrictGreaterThan', 'StrictLessThan', 'vectorize',
           'Lambda', 'WildFunction', 'Derivative', 'diff', 'FunctionClass',
           'Function', 'Subs', 'expand', 'PoleError', 'count_ops',
           'expand_mul', 'expand_log', 'expand_func', 'expand_trig',
           'expand_complex', 'expand_multinomial', 'nfloat',
           'expand_power_base', 'expand_power_exp', 'PrecisionExhausted', 'N',
           'Tuple', 'Dict', 'gcd_terms', 'factor_terms', 'factor_nc',
           'evaluate', 'Catalan', 'EulerGamma', 'GoldenRatio', 'Relational')
