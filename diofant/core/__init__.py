"""Core module. Provides the basic operations needed in diofant.
"""

from .sympify import sympify, SympifyError  # noqa: F401
from .cache import cacheit  # noqa: F401
from .basic import Basic, Atom, preorder_traversal  # noqa: F401
from .singleton import S  # noqa: F401
from .expr import Expr, AtomicExpr  # noqa: F401
from .symbol import Symbol, Wild, Dummy, symbols, var  # noqa: F401
from .numbers import (Number, Float, Rational, Integer,  # noqa: F401
                      NumberSymbol, igcd, ilcm, seterr, E, I, nan, oo,
                      pi, zoo, comp, mod_inverse)
from .power import Pow, integer_nthroot  # noqa: F401
from .mul import Mul, prod  # noqa: F401
from .add import Add  # noqa: F401
from .mod import Mod  # noqa: F401
from .relational import (Rel, Eq, Ne, Lt, Le, Gt, Ge, Equality,  # noqa: F401
                         GreaterThan, LessThan, Unequality, StrictGreaterThan,
                         StrictLessThan)
from .multidimensional import vectorize  # noqa: F401
from .function import (Lambda, WildFunction, Derivative, diff,  # noqa: F401
                       FunctionClass, Function, Subs, expand, PoleError,
                       count_ops, expand_mul, expand_log, expand_func,
                       expand_trig, expand_complex, expand_multinomial, nfloat,
                       expand_power_base, expand_power_exp)
from .evalf import PrecisionExhausted, N  # noqa: F401
from .containers import Tuple, Dict  # noqa: F401
from .exprtools import gcd_terms, factor_terms, factor_nc  # noqa: F401
from .evaluate import evaluate  # noqa: F401

# expose singletons
Catalan = S.Catalan
EulerGamma = S.EulerGamma
GoldenRatio = S.GoldenRatio
