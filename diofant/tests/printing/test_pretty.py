import sys
import warnings
from io import StringIO

import pytest

from diofant import (FF, QQ, RR, ZZ, Add, Adjoint, And, Basic, Chi, Ci,
                     Complement, Contains, Derivative, Dict, DiracDelta, E, Ei,
                     Eq, Equivalent, EulerGamma, FiniteSet, Float, Function,
                     Ge, GoldenRatio, Gt, I, Implies, Integer, Integral,
                     Intersection, Interval, Inverse, KroneckerDelta, Lambda,
                     Le, Limit, Lt, Matrix, MatrixSymbol, Mod, Mul, Nand, Ne,
                     Nor, Not, O, Or, Piecewise, Pow, Product, Range, Rational,
                     Ray, RealField, RootOf, RootSum, S, Segment, Shi, Si,
                     Subs, Sum, Symbol, SymmetricDifference, Trace, Transpose,
                     Tuple, Union, Xor, atan2, binomial, catalan, cbrt,
                     ceiling, conjugate, cos, elliptic_e, elliptic_f,
                     elliptic_k, elliptic_pi, euler, exp, expint, factorial,
                     factorial2, floor, gamma, grlex, groebner, hyper, ilex,
                     log, lowergamma, meijerg, oo, pi, pprint, root, sin, sqrt,
                     subfactorial, symbols, tan, uppergamma)
from diofant.abc import (a, b, c, d, e, f, k, l, lamda, m, n, phi, t, theta, w,
                         x, y, z)
from diofant.core.trace import Tr
from diofant.diffgeom import BaseVectorField
from diofant.diffgeom.rn import R2, R2_r
from diofant.printing.pretty import pretty as xpretty
from diofant.printing.pretty.pretty_symbology import U, xobj
from diofant.stats import Die, Exponential, Normal, pspace, where
from diofant.tensor import (ImmutableDenseNDimArray, ImmutableSparseNDimArray,
                            MutableDenseNDimArray, MutableSparseNDimArray,
                            tensorproduct)


__all__ = ()


"""
Expressions whose pretty-printing is tested here:
(A '#' to the right of an expression indicates that its various acceptable
orderings are accounted for by the tests.)


BASIC EXPRESSIONS:

oo
(x**2)
1/x
y*x**-2
x**Rational(-5,2)
(-2)**x
Pow(3, 1, evaluate=False)
(x**2 + x + 1)  #
1-x  #
1-2*x  #
x/y
-x/y
(x+2)/y  #
(1+x)*y  #3
-5*x/(x+10)  # correct placement of negative sign
1 - Rational(3,2)*(x+1)
-(-x + 5)*(-x - 2*sqrt(2) + 5) - (-y + 5)*(-y + 5) # issue sympy/sympy#5524


ORDERING:

x**2 + x + 1
1 - x
1 - 2*x
2*x**4 + y**2 - x**2 + y**3


RELATIONAL:

Eq(x, y)
Lt(x, y)
Gt(x, y)
Le(x, y)
Ge(x, y)
Ne(x/(y+1), y**2)  #


RATIONAL NUMBERS:

y*x**-2
y**Rational(3,2) * x**Rational(-5,2)
sin(x)**3/tan(x)**2


FUNCTIONS (ABS, CONJ, EXP, FUNCTION BRACES, FACTORIAL, FLOOR, CEILING):

(2*x + exp(x))  #
abs(x)
abs(x/(x**2+1)) #
abs(1 / (y - abs(x)))
factorial(n)
factorial(2*n)
subfactorial(n)
subfactorial(2*n)
factorial(factorial(factorial(n)))
factorial(n+1) #
conjugate(x)
conjugate(f(x+1)) #
f(x)
f(x, y)
f(x/(y+1), y) #
f(x**x**x**x**x**x)
sin(x)**2
conjugate(a+b*I)
conjugate(exp(a+b*I))
conjugate( f(1 + conjugate(f(x))) ) #
f(x/(y+1), y)  # denom of first arg
floor(1 / (y - floor(x)))
ceiling(1 / (y - ceiling(x)))


SQRT:

sqrt(2)
cbrt(2)
root(2, 1000)
sqrt(x**2 + 1)
cbrt(1 + sqrt(5))
2**(1/x)
sqrt(2+pi)
root(2+(1+x**2)/(2+x), 4)+(1+root(x,1000))/sqrt(3+x**2)


DERIVATIVES:

Derivative(log(x), x, evaluate=False)
Derivative(log(x), x, evaluate=False) + x  #
Derivative(log(x) + x**2, x, y, evaluate=False)
Derivative(2*x*y, y, x, evaluate=False) + x**2  #
beta(alpha).diff(alpha)


INTEGRALS:

Integral(log(x), x)
Integral(x**2, x)
Integral((sin(x))**2 / (tan(x))**2)
Integral(x**(2**x), x)
Integral(x**2, (x,1,2))
Integral(x**2, (x,Rational(1,2),10))
Integral(x**2*y**2, x,y)
Integral(x**2, (x, None, 1))
Integral(x**2, (x, 1, None))
Integral(sin(theta)/cos(phi), (theta,0,phi), (phi, 0, 2*pi))


MATRICES:

Matrix([[x**2+1, 1], [y, x+y]])  #
Matrix([[x/y, y, theta], [0, exp(I*k*phi), 1]])


PIECEWISE:

Piecewise((x,x<1),(x**2,True))


SEQUENCES (TUPLES, LISTS, DICTIONARIES):

()
[]
{}
(1/x,)
[x**2, 1/x, x, y, sin(theta)**2/cos(phi)**2]
(x**2, 1/x, x, y, sin(theta)**2/cos(phi)**2)
{x: sin(x)}
{1/x: 1/y, x: sin(x)**2}  #
[x**2]
(x**2,)
{x**2: 1}


LIMITS:

Limit(x, x, oo)
Limit(x**2, x, 0)
Limit(1/x, x, 0)
Limit(sin(x)/x, x, 0)


SUBS:

Subs(f(x), (x, ph**2))
Subs(f(x).diff(x), (x, 0))
Subs(f(x).diff(x)/y, (x, 0), (y, Rational(1, 2)))


ORDER:

O(1)
O(1/x)
O(x**2 + y**2)
"""


def pretty(expr, order=None):
    """ASCII pretty-printing"""
    return xpretty(expr, order=order, use_unicode=False, wrap_line=False)


def upretty(expr, order=None):
    """Unicode pretty-printing"""
    return xpretty(expr, order=order, use_unicode=True, wrap_line=False)


def test_pretty_ascii_str():
    assert pretty( 'xxx' ) == 'xxx'
    assert pretty( "xxx'xxx" ) == "xxx'xxx"
    assert pretty( 'xxx"xxx' ) == 'xxx"xxx'
    assert pretty( 'xxx\nxxx' ) == 'xxx\nxxx'


def test_pretty_unicode_str():
    assert pretty( 'xxx' ) == 'xxx'
    assert pretty( 'xxx' ) == 'xxx'
    assert pretty( "xxx'xxx" ) == "xxx'xxx"
    assert pretty( 'xxx"xxx' ) == 'xxx"xxx'
    assert pretty( 'xxx\nxxx' ) == 'xxx\nxxx'


def test_upretty_greek():
    assert upretty( oo ) == '∞'
    assert upretty( Symbol('alpha^+_1') ) == 'α⁺₁'
    assert upretty( Symbol('beta') ) == 'β'
    assert upretty(Symbol('lambda')) == 'λ'


def test_upretty_multiindex():
    assert upretty( Symbol('beta12') ) == 'β₁₂'
    assert upretty( Symbol('Y00') ) == 'Y₀₀'
    assert upretty( Symbol('Y_00') ) == 'Y₀₀'
    assert upretty( Symbol('F^+-') ) == 'F⁺⁻'


def test_upretty_sub_super():
    assert upretty( Symbol('beta_1_2') ) == 'β₁ ₂'
    assert upretty( Symbol('beta^1^2') ) == 'β¹ ²'
    assert upretty( Symbol('beta_1^2') ) == 'β²₁'
    assert upretty( Symbol('beta_10_20') ) == 'β₁₀ ₂₀'
    assert upretty( Symbol('beta_ax_gamma^i') ) == 'βⁱₐₓ ᵧ'
    assert upretty( Symbol('F^1^2_3_4') ) == 'F¹ ²₃ ₄'
    assert upretty( Symbol('F_1_2^3^4') ) == 'F³ ⁴₁ ₂'
    assert upretty( Symbol('F_1_2_3_4') ) == 'F₁ ₂ ₃ ₄'
    assert upretty( Symbol('F^1^2^3^4') ) == 'F¹ ² ³ ⁴'


def test_upretty_subs_missing_in_24():
    assert upretty( Symbol('F_beta') ) == 'Fᵦ'
    assert upretty( Symbol('F_gamma') ) == 'Fᵧ'
    assert upretty( Symbol('F_rho') ) == 'Fᵨ'
    assert upretty( Symbol('F_phi') ) == 'Fᵩ'
    assert upretty( Symbol('F_chi') ) == 'Fᵪ'

    assert upretty( Symbol('F_a') ) == 'Fₐ'
    assert upretty( Symbol('F_e') ) == 'Fₑ'
    assert upretty( Symbol('F_i') ) == 'Fᵢ'
    assert upretty( Symbol('F_o') ) == 'Fₒ'
    assert upretty( Symbol('F_u') ) == 'Fᵤ'
    assert upretty( Symbol('F_r') ) == 'Fᵣ'
    assert upretty( Symbol('F_v') ) == 'Fᵥ'
    assert upretty( Symbol('F_x') ) == 'Fₓ'


def test_missing_in_2X_sympyissue_9047():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert upretty( Symbol('F_h') ) == 'Fₕ'
        assert upretty( Symbol('F_k') ) == 'Fₖ'
        assert upretty( Symbol('F_l') ) == 'Fₗ'
        assert upretty( Symbol('F_m') ) == 'Fₘ'
        assert upretty( Symbol('F_n') ) == 'Fₙ'
        assert upretty( Symbol('F_p') ) == 'Fₚ'
        assert upretty( Symbol('F_s') ) == 'Fₛ'
        assert upretty( Symbol('F_t') ) == 'Fₜ'


def test_upretty_modifiers():
    # Accents
    assert upretty( Symbol('Fmathring') ) == 'F̊'
    assert upretty( Symbol('Fddddot') ) == 'F̈̈'
    assert upretty( Symbol('Fdddot') ) == 'F̈̇'
    assert upretty( Symbol('Fddot') ) == 'F̈'
    assert upretty( Symbol('Fdot') ) == 'Ḟ'
    assert upretty( Symbol('Fcheck') ) == 'F̌'
    assert upretty( Symbol('Fbreve') ) == 'F̆'
    assert upretty( Symbol('Facute') ) == 'F́'
    assert upretty( Symbol('Fgrave') ) == 'F̀'
    assert upretty( Symbol('Ftilde') ) == 'F̃'
    assert upretty( Symbol('Fhat') ) == 'F̂'
    assert upretty( Symbol('Fbar') ) == 'F̅'
    assert upretty( Symbol('Fvec') ) == 'F⃗'
    assert upretty( Symbol('Fprime') ) == 'F′'
    assert upretty( Symbol('Fprm') ) == 'F′'
    # No faces are actually implemented, but test to make sure the modifiers are stripped
    assert upretty( Symbol('Fbold') ) == 'Fbold'
    assert upretty( Symbol('Fbm') ) == 'Fbm'
    assert upretty( Symbol('Fcal') ) == 'Fcal'
    assert upretty( Symbol('Fscr') ) == 'Fscr'
    assert upretty( Symbol('Ffrak') ) == 'Ffrak'
    # Brackets
    assert upretty( Symbol('Fnorm') ) == '‖F‖'
    assert upretty( Symbol('Favg') ) == '⟨F⟩'
    assert upretty( Symbol('Fabs') ) == '|F|'
    assert upretty( Symbol('Fmag') ) == '|F|'
    # Combinations
    assert upretty( Symbol('xvecdot') ) == 'x⃗̇'
    assert upretty( Symbol('xDotVec') ) == 'ẋ⃗'
    assert upretty( Symbol('xHATNorm') ) == '‖x̂‖'
    assert upretty( Symbol('xMathring_yCheckPRM__zbreveAbs') ) == 'x̊_y̌′__|z̆|'
    assert upretty( Symbol('alphadothat_nVECDOT__tTildePrime') ) == 'α̇̂_n⃗̇__t̃′'
    assert upretty( Symbol('x_dot') ) == 'x_dot'
    assert upretty( Symbol('x__dot') ) == 'x__dot'


def test_symbology():
    pytest.raises(ValueError, lambda: xobj('', 0))
    assert U('NO SUCH SYMBOL') is None


def test_pretty_atom():
    assert upretty(S.Rationals) == 'ℚ'


def test_pretty_basic():
    assert pretty(-Rational(1, 2)) == '-1/2'
    assert pretty(-Rational(13, 22)) == \
        """\
-13 \n\
----\n\
 22 \
"""
    expr = oo
    ascii_str = \
        """\
oo\
"""
    ucode_str = \
        """\
∞\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x**2)
    ascii_str = \
        """\
 2\n\
x \
"""
    ucode_str = \
        """\
 2\n\
x \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 1/x
    ascii_str = \
        """\
1\n\
-\n\
x\
"""
    ucode_str = \
        """\
1\n\
─\n\
x\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # not the same as 1/x
    expr = x**-1.0
    ascii_str = \
        """\
 -1.0\n\
x    \
"""
    ucode_str = \
        """\
 -1.0\n\
x    \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # see issue sympy/sympy#2860
    expr = Pow(Integer(2), -1.0, evaluate=False)
    ascii_str = \
        """\
 -1.0\n\
2    \
"""
    ucode_str = \
        """\
 -1.0\n\
2    \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = y*x**-2
    ascii_str = \
        """\
y \n\
--\n\
 2\n\
x \
"""
    ucode_str = \
        """\
y \n\
──\n\
 2\n\
x \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x**Rational(-5, 2)
    ascii_str = \
        """\
 1  \n\
----\n\
 5/2\n\
x   \
"""
    ucode_str = \
        """\
 1  \n\
────\n\
 5/2\n\
x   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (-2)**x
    ascii_str = \
        """\
    x\n\
(-2) \
"""
    ucode_str = \
        """\
    x\n\
(-2) \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # See issue sympy/sympy#4923
    expr = Pow(3, 1, evaluate=False)
    ascii_str = \
        """\
 1\n\
3 \
"""
    ucode_str = \
        """\
 1\n\
3 \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x**2 + x + 1)
    ascii_str_1 = \
        """\
         2\n\
1 + x + x \
"""
    ascii_str_2 = \
        """\
 2        \n\
x  + x + 1\
"""
    ascii_str_3 = \
        """\
 2        \n\
x  + 1 + x\
"""
    ucode_str_1 = \
        """\
         2\n\
1 + x + x \
"""
    ucode_str_2 = \
        """\
 2        \n\
x  + x + 1\
"""
    ucode_str_3 = \
        """\
 2        \n\
x  + 1 + x\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2, ascii_str_3]
    assert upretty(expr) in [ucode_str_1, ucode_str_2, ucode_str_3]

    expr = 1 - x
    ascii_str_1 = \
        """\
1 - x\
"""
    ascii_str_2 = \
        """\
-x + 1\
"""
    ucode_str_1 = \
        """\
1 - x\
"""
    ucode_str_2 = \
        """\
-x + 1\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = 1 - 2*x
    ascii_str_1 = \
        """\
1 - 2*x\
"""
    ascii_str_2 = \
        """\
-2*x + 1\
"""
    ucode_str_1 = \
        """\
1 - 2⋅x\
"""
    ucode_str_2 = \
        """\
-2⋅x + 1\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = x/y
    ascii_str = \
        """\
x\n\
-\n\
y\
"""
    ucode_str = \
        """\
x\n\
─\n\
y\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = -x/y
    ascii_str = \
        """\
-x \n\
---\n\
 y \
"""
    ucode_str = \
        """\
-x \n\
───\n\
 y \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x + 2)/y
    ascii_str_1 = \
        """\
2 + x\n\
-----\n\
  y  \
"""
    ascii_str_2 = \
        """\
x + 2\n\
-----\n\
  y  \
"""
    ucode_str_1 = \
        """\
2 + x\n\
─────\n\
  y  \
"""
    ucode_str_2 = \
        """\
x + 2\n\
─────\n\
  y  \
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = (1 + x)*y
    ascii_str_1 = \
        """\
y*(1 + x)\
"""
    ascii_str_2 = \
        """\
(1 + x)*y\
"""
    ascii_str_3 = \
        """\
y*(x + 1)\
"""
    ucode_str_1 = \
        """\
y⋅(1 + x)\
"""
    ucode_str_2 = \
        """\
(1 + x)⋅y\
"""
    ucode_str_3 = \
        """\
y⋅(x + 1)\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2, ascii_str_3]
    assert upretty(expr) in [ucode_str_1, ucode_str_2, ucode_str_3]

    # Test for correct placement of the negative sign
    expr = -5*x/(x + 10)
    ascii_str_1 = \
        """\
-5*x  \n\
------\n\
10 + x\
"""
    ascii_str_2 = \
        """\
-5*x  \n\
------\n\
x + 10\
"""
    ucode_str_1 = \
        """\
-5⋅x  \n\
──────\n\
10 + x\
"""
    ucode_str_2 = \
        """\
-5⋅x  \n\
──────\n\
x + 10\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = -Rational(1, 2) - 3*x
    ascii_str = \
        """\
-3*x - 1/2\
"""
    ucode_str = \
        """\
-3⋅x - 1/2\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Rational(1, 2) - 3*x
    ascii_str = \
        """\
-3*x + 1/2\
"""
    ucode_str = \
        """\
-3⋅x + 1/2\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = -Rational(1, 2) - 3*x/2
    ascii_str = \
        """\
  3*x   1\n\
- --- - -\n\
   2    2\
"""
    ucode_str = \
        """\
  3⋅x   1\n\
- ─── - ─\n\
   2    2\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Rational(1, 2) - 3*x/2
    ascii_str = \
        """\
  3*x   1\n\
- --- + -\n\
   2    2\
"""
    ucode_str = \
        """\
  3⋅x   1\n\
- ─── + ─\n\
   2    2\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # issue sympy/sympy#7927
    e = sin(x/2)**cos(x/2)
    ucode_str = \
        """\
           ⎛x⎞\n\
        cos⎜─⎟\n\
           ⎝2⎠\n\
⎛   ⎛x⎞⎞      \n\
⎜sin⎜─⎟⎟      \n\
⎝   ⎝2⎠⎠      \
"""
    assert upretty(e) == ucode_str
    e = sin(x)**Rational(11, 13)
    ucode_str = \
        """\
        11\n\
        ──\n\
        13\n\
(sin(x))  \
"""
    assert upretty(e) == ucode_str

    # issue sympy/sympy#7117
    # See also issue sympy/sympy#5031 (hence the evaluate=False in these).
    e = Eq(x + 1, x/2)
    q = Mul(2, e, evaluate=False)
    assert upretty(q) == """\
  ⎛        x⎞\n\
2⋅⎜x + 1 = ─⎟\n\
  ⎝        2⎠\
"""
    q = Add(e, 6, evaluate=False)
    assert upretty(q) == """\
    ⎛        x⎞\n\
6 + ⎜x + 1 = ─⎟\n\
    ⎝        2⎠\
"""
    q = Pow(e, 2, evaluate=False)
    assert upretty(q) == """\
           2\n\
⎛        x⎞ \n\
⎜x + 1 = ─⎟ \n\
⎝        2⎠ \
"""
    e2 = Eq(x, 2)
    q = Mul(e, e2, evaluate=False)
    assert upretty(q) == """\
⎛        x⎞        \n\
⎜x + 1 = ─⎟⋅(x = 2)\n\
⎝        2⎠        \
"""


def test_negative_fractions():
    expr = -x/y
    ascii_str =\
        """\
-x \n\
---\n\
 y \
"""
    ucode_str =\
        """\
-x \n\
───\n\
 y \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -x*z/y
    ascii_str =\
        """\
-x*z \n\
-----\n\
  y  \
"""
    ucode_str =\
        """\
-x⋅z \n\
─────\n\
  y  \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = x**2/y
    ascii_str =\
        """\
 2\n\
x \n\
--\n\
y \
"""
    ucode_str =\
        """\
 2\n\
x \n\
──\n\
y \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -x**2/y
    ascii_str =\
        """\
  2 \n\
-x  \n\
----\n\
 y  \
"""
    ucode_str =\
        """\
  2 \n\
-x  \n\
────\n\
 y  \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -x/(y*z)
    ascii_str =\
        """\
-x \n\
---\n\
y*z\
"""
    ucode_str =\
        """\
-x \n\
───\n\
y⋅z\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -a/y**2
    ascii_str =\
        """\
-a \n\
---\n\
  2\n\
 y \
"""
    ucode_str =\
        """\
-a \n\
───\n\
  2\n\
 y \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = y**(-a/b)
    ascii_str =\
        """\
 -a \n\
 ---\n\
  b \n\
y   \
"""
    ucode_str =\
        """\
 -a \n\
 ───\n\
  b \n\
y   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -1/y**2
    ascii_str =\
        """\
-1 \n\
---\n\
  2\n\
 y \
"""
    ucode_str =\
        """\
-1 \n\
───\n\
  2\n\
 y \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = -10/b**2
    ascii_str =\
        """\
-10 \n\
----\n\
  2 \n\
 b  \
"""
    ucode_str =\
        """\
-10 \n\
────\n\
  2 \n\
 b  \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str
    expr = Rational(-200, 37)
    ascii_str =\
        """\
-200 \n\
-----\n\
  37 \
"""
    ucode_str =\
        """\
-200 \n\
─────\n\
  37 \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_sympyissue_5524():
    assert pretty(-(-x + 5)*(-x - 2*sqrt(2) + 5) - (-y + 5)*(-y + 5)) == \
        """\
        /         ___    \\           2\n\
(x - 5)*\\-x - 2*\\/ 2  + 5/ - (-y + 5) \
"""

    assert upretty(-(-x + 5)*(-x - 2*sqrt(2) + 5) - (-y + 5)*(-y + 5)) == \
        """\
        ⎛         ___    ⎞           2\n\
(x - 5)⋅⎝-x - 2⋅╲╱ 2  + 5⎠ - (-y + 5) \
"""


def test_pretty_ordering():
    assert pretty(x**2 + x + 1, order='lex') == \
        """\
 2        \n\
x  + x + 1\
"""
    assert pretty(x**2 + x + 1, order='rev-lex') == \
        """\
         2\n\
1 + x + x \
"""
    assert pretty(1 - x, order='lex') == '-x + 1'
    assert pretty(1 - x, order='rev-lex') == '1 - x'

    assert pretty(1 - 2*x, order='lex') == '-2*x + 1'
    assert pretty(1 - 2*x, order='rev-lex') == '1 - 2*x'

    f = 2*x**4 + y**2 - x**2 + y**3
    assert pretty(f, order=None) == \
        """\
   4    2    3    2\n\
2*x  - x  + y  + y \
"""
    assert pretty(f, order='lex') == \
        """\
   4    2    3    2\n\
2*x  - x  + y  + y \
"""
    assert pretty(f, order='rev-lex') == \
        """\
 2    3    2      4\n\
y  + y  - x  + 2*x \
"""

    expr = x - x**3/6 + x**5/120 + O(x**6)
    ascii_str = \
        """\
     3     5        \n\
    x     x     / 6\\\n\
x - -- + --- + O\\x /\n\
    6    120        \
"""
    ucode_str = \
        """\
     3     5        \n\
    x     x     ⎛ 6⎞\n\
x - ── + ─── + O⎝x ⎠\n\
    6    120        \
"""
    assert pretty(expr, order=None) == ascii_str
    assert upretty(expr, order=None) == ucode_str

    assert pretty(expr, order='lex') == ascii_str
    assert upretty(expr, order='lex') == ucode_str

    assert pretty(expr, order='rev-lex') == ascii_str
    assert upretty(expr, order='rev-lex') == ucode_str


def test_EulerGamma():
    assert pretty(EulerGamma) == str(EulerGamma) == 'EulerGamma'
    assert upretty(EulerGamma) == 'γ'


def test_GoldenRatio():
    assert pretty(GoldenRatio) == str(GoldenRatio) == 'GoldenRatio'
    assert upretty(GoldenRatio) == 'φ'


def test_pretty_relational():
    expr = Eq(x, y)
    ascii_str = \
        """\
x = y\
"""
    ucode_str = \
        """\
x = y\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Lt(x, y)
    ascii_str = \
        """\
x < y\
"""
    ucode_str = \
        """\
x < y\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Gt(x, y)
    ascii_str = \
        """\
x > y\
"""
    ucode_str = \
        """\
x > y\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Le(x, y)
    ascii_str = \
        """\
x <= y\
"""
    ucode_str = \
        """\
x ≤ y\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Ge(x, y)
    ascii_str = \
        """\
x >= y\
"""
    ucode_str = \
        """\
x ≥ y\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Ne(x/(y + 1), y**2)
    ascii_str_1 = \
        """\
  x       2\n\
----- != y \n\
1 + y      \
"""
    ascii_str_2 = \
        """\
  x       2\n\
----- != y \n\
y + 1      \
"""
    ucode_str_1 = \
        """\
  x      2\n\
───── ≠ y \n\
1 + y     \
"""
    ucode_str_2 = \
        """\
  x      2\n\
───── ≠ y \n\
y + 1     \
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]


def test_pretty_rational():
    expr = y*x**-2
    ascii_str = \
        """\
y \n\
--\n\
 2\n\
x \
"""
    ucode_str = \
        """\
y \n\
──\n\
 2\n\
x \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = y**Rational(3, 2) * x**Rational(-5, 2)
    ascii_str = \
        """\
 3/2\n\
y   \n\
----\n\
 5/2\n\
x   \
"""
    ucode_str = \
        """\
 3/2\n\
y   \n\
────\n\
 5/2\n\
x   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = sin(x)**3/tan(x)**2
    ascii_str = \
        """\
   3   \n\
sin (x)\n\
-------\n\
   2   \n\
tan (x)\
"""
    ucode_str = \
        """\
   3   \n\
sin (x)\n\
───────\n\
   2   \n\
tan (x)\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_functions():
    """Tests for Abs, conjugate, exp, function braces, and factorial."""
    expr = (2*x + exp(x))
    ascii_str_1 = \
        """\
       x\n\
2*x + E \
"""
    ascii_str_2 = \
        """\
 x      \n\
E  + 2*x\
"""
    ucode_str_1 = \
        """\
       x\n\
2⋅x + ℯ \
"""
    ucode_str_2 = \
        """\
 x      \n\
ℯ  + 2⋅x\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = abs(x)
    ascii_str = \
        """\
|x|\
"""
    ucode_str = \
        """\
│x│\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = abs(x/(x**2 + 1))
    ascii_str_1 = \
        """\
|  x   |\n\
|------|\n\
|     2|\n\
|1 + x |\
"""
    ascii_str_2 = \
        """\
|  x   |\n\
|------|\n\
| 2    |\n\
|x  + 1|\
"""
    ucode_str_1 = \
        """\
│  x   │\n\
│──────│\n\
│     2│\n\
│1 + x │\
"""
    ucode_str_2 = \
        """\
│  x   │\n\
│──────│\n\
│ 2    │\n\
│x  + 1│\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = abs(1 / (y - abs(x)))
    ascii_str = \
        """\
|   1   |\n\
|-------|\n\
|y - |x||\
"""
    ucode_str = \
        """\
│   1   │\n\
│───────│\n\
│y - │x││\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    n = Symbol('n', integer=True)
    expr = factorial(n)
    ascii_str = \
        """\
n!\
"""
    ucode_str = \
        """\
n!\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial(2*n)
    ascii_str = \
        """\
(2*n)!\
"""
    ucode_str = \
        """\
(2⋅n)!\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial(factorial(factorial(n)))
    ascii_str = \
        """\
((n!)!)!\
"""
    ucode_str = \
        """\
((n!)!)!\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial(n + 1)
    ascii_str_1 = \
        """\
(1 + n)!\
"""
    ascii_str_2 = \
        """\
(n + 1)!\
"""
    ucode_str_1 = \
        """\
(1 + n)!\
"""
    ucode_str_2 = \
        """\
(n + 1)!\
"""

    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = subfactorial(n)
    ascii_str = \
        """\
!n\
"""
    ucode_str = \
        """\
!n\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = subfactorial(2*n)
    ascii_str = \
        """\
!(2*n)\
"""
    ucode_str = \
        """\
!(2⋅n)\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    n = Symbol('n', integer=True)
    expr = factorial2(n)
    ascii_str = \
        """\
n!!\
"""
    ucode_str = \
        """\
n!!\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial2(2*n)
    ascii_str = \
        """\
(2*n)!!\
"""
    ucode_str = \
        """\
(2⋅n)!!\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial2(factorial2(factorial2(n)))
    ascii_str = \
        """\
((n!!)!!)!!\
"""
    ucode_str = \
        """\
((n!!)!!)!!\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = factorial2(n + 1)
    ascii_str_1 = \
        """\
(1 + n)!!\
"""
    ascii_str_2 = \
        """\
(n + 1)!!\
"""
    ucode_str_1 = \
        """\
(1 + n)!!\
"""
    ucode_str_2 = \
        """\
(n + 1)!!\
"""

    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = 2*binomial(n, k)
    ascii_str = \
        """\
  /n\\\n\
2*| |\n\
  \\k/\
"""
    ucode_str = \
        """\
  ⎛n⎞\n\
2⋅⎜ ⎟\n\
  ⎝k⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 2*binomial(2*n, k)
    ascii_str = \
        """\
  /2*n\\\n\
2*|   |\n\
  \\ k /\
"""
    ucode_str = \
        """\
  ⎛2⋅n⎞\n\
2⋅⎜   ⎟\n\
  ⎝ k ⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 2*binomial(n**2, k)
    ascii_str = \
        """\
  / 2\\\n\
  |n |\n\
2*|  |\n\
  \\k /\
"""
    ucode_str = \
        """\
  ⎛ 2⎞\n\
  ⎜n ⎟\n\
2⋅⎜  ⎟\n\
  ⎝k ⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = catalan(n)
    ascii_str = \
        """\
C \n\
 n\
"""
    ucode_str = \
        """\
C \n\
 n\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = conjugate(x)
    ascii_str = \
        """\
_\n\
x\
"""
    ucode_str = \
        """\
_\n\
x\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    f = Function('f')
    expr = conjugate(f(x + 1))
    ascii_str_1 = \
        """\
________\n\
f(1 + x)\
"""
    ascii_str_2 = \
        """\
________\n\
f(x + 1)\
"""
    ucode_str_1 = \
        """\
________\n\
f(1 + x)\
"""
    ucode_str_2 = \
        """\
________\n\
f(x + 1)\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = f(x)
    ascii_str = \
        """\
f(x)\
"""
    ucode_str = \
        """\
f(x)\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = f(x, y)
    ascii_str = \
        """\
f(x, y)\
"""
    ucode_str = \
        """\
f(x, y)\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = f(x/(y + 1), y)
    ascii_str_1 = \
        """\
 /  x     \\\n\
f|-----, y|\n\
 \\1 + y   /\
"""
    ascii_str_2 = \
        """\
 /  x     \\\n\
f|-----, y|\n\
 \\y + 1   /\
"""
    ucode_str_1 = \
        """\
 ⎛  x     ⎞\n\
f⎜─────, y⎟\n\
 ⎝1 + y   ⎠\
"""
    ucode_str_2 = \
        """\
 ⎛  x     ⎞\n\
f⎜─────, y⎟\n\
 ⎝y + 1   ⎠\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = f(x**x**x**x**x**x)
    ascii_str = \
        """\
 / / / / / x\\\\\\\\\\
 | | | | \\x /||||
 | | | \\x    /|||
 | | \\x       /||
 | \\x          /|
f\\x             /\
"""
    ucode_str = \
        """\
 ⎛ ⎛ ⎛ ⎛ ⎛ x⎞⎞⎞⎞⎞
 ⎜ ⎜ ⎜ ⎜ ⎝x ⎠⎟⎟⎟⎟
 ⎜ ⎜ ⎜ ⎝x    ⎠⎟⎟⎟
 ⎜ ⎜ ⎝x       ⎠⎟⎟
 ⎜ ⎝x          ⎠⎟
f⎝x             ⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = sin(x)**2
    ascii_str = \
        """\
   2   \n\
sin (x)\
"""
    ucode_str = \
        """\
   2   \n\
sin (x)\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = conjugate(a + b*I)
    ascii_str = \
        """\
_     _\n\
a - I*b\
"""
    ucode_str = \
        """\
_     _\n\
a - ⅈ⋅b\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = conjugate(exp(a + b*I))
    ascii_str = \
        """\
 _     _\n\
 a - I*b\n\
E       \
"""
    ucode_str = \
        """\
 _     _\n\
 a - ⅈ⋅b\n\
ℯ       \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = conjugate( f(1 + conjugate(f(x))) )
    ascii_str_1 = \
        """\
___________\n\
 /    ____\\\n\
f\\1 + f(x)/\
"""
    ascii_str_2 = \
        """\
___________\n\
 /____    \\\n\
f\\f(x) + 1/\
"""
    ucode_str_1 = \
        """\
___________\n\
 ⎛    ____⎞\n\
f⎝1 + f(x)⎠\
"""
    ucode_str_2 = \
        """\
___________\n\
 ⎛____    ⎞\n\
f⎝f(x) + 1⎠\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = f(x/(y + 1), y)
    ascii_str_1 = \
        """\
 /  x     \\\n\
f|-----, y|\n\
 \\1 + y   /\
"""
    ascii_str_2 = \
        """\
 /  x     \\\n\
f|-----, y|\n\
 \\y + 1   /\
"""
    ucode_str_1 = \
        """\
 ⎛  x     ⎞\n\
f⎜─────, y⎟\n\
 ⎝1 + y   ⎠\
"""
    ucode_str_2 = \
        """\
 ⎛  x     ⎞\n\
f⎜─────, y⎟\n\
 ⎝y + 1   ⎠\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = floor(1 / (y - floor(x)))
    ascii_str = \
        """\
     /     1      \\\n\
floor|------------|\n\
     \\y - floor(x)/\
"""
    ucode_str = \
        """\
⎢   1   ⎥\n\
⎢───────⎥\n\
⎣y - ⌊x⌋⎦\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = ceiling(1 / (y - ceiling(x)))
    ascii_str = \
        """\
       /      1       \\\n\
ceiling|--------------|\n\
       \\y - ceiling(x)/\
"""
    ucode_str = \
        """\
⎡   1   ⎤\n\
⎢───────⎥\n\
⎢y - ⌈x⌉⎥\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = euler(n)
    ascii_str = \
        """\
E \n\
 n\
"""
    ucode_str = \
        """\
E \n\
 n\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = euler(1/(1 + 1/(1 + 1/n)))
    ascii_str = \
        """\
E         \n\
     1    \n\
 ---------\n\
       1  \n\
 1 + -----\n\
         1\n\
     1 + -\n\
         n\
"""

    ucode_str = \
        """\
E         \n\
     1    \n\
 ─────────\n\
       1  \n\
 1 + ─────\n\
         1\n\
     1 + ─\n\
         n\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_sqrt():
    expr = sqrt(2)
    ascii_str = \
        """\
  ___\n\
\\/ 2 \
"""
    ucode_str = \
        """\
  ___\n\
╲╱ 2 \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = cbrt(2)
    ascii_str = \
        """\
3 ___\n\
\\/ 2 \
"""
    ucode_str = \
        """\
3 ___\n\
╲╱ 2 \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = root(2, 1000)
    ascii_str = \
        """\
1000___\n\
  \\/ 2 \
"""
    ucode_str = \
        """\
1000___\n\
  ╲╱ 2 \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = sqrt(x**2 + 1)
    ascii_str = \
        """\
   ________\n\
  /  2     \n\
\\/  x  + 1 \
"""
    ucode_str = \
        """\
   ________\n\
  ╱  2     \n\
╲╱  x  + 1 \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = cbrt(1 + sqrt(5))
    ascii_str = \
        """\
   ___________\n\
3 /       ___ \n\
\\/  1 + \\/ 5  \
"""
    ucode_str = \
        """\
   ___________\n\
3 ╱       ___ \n\
╲╱  1 + ╲╱ 5  \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = 2**(1/x)
    ascii_str = \
        """\
x ___\n\
\\/ 2 \
"""
    ucode_str = \
        """\
x ___\n\
╲╱ 2 \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = sqrt(2 + pi)
    ascii_str = \
        """\
  ________\n\
\\/ 2 + pi \
"""
    ucode_str = \
        """\
  _______\n\
╲╱ 2 + π \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = root(2 + (1 + x**2)/(2 + x), 4) + (1 + root(x, 1000))/sqrt(3 + x**2)
    ascii_str = \
        """\
     ____________              \n\
    /      2        1000___    \n\
   /      x  + 1      \\/ x  + 1\n\
4 /   2 + ------  + -----------\n\
\\/        x + 2        ________\n\
                      /  2     \n\
                    \\/  x  + 3 \
"""
    ucode_str = \
        """\
     ____________              \n\
    ╱      2        1000___    \n\
   ╱      x  + 1      ╲╱ x  + 1\n\
4 ╱   2 + ──────  + ───────────\n\
╲╱        x + 2        ________\n\
                      ╱  2     \n\
                    ╲╱  x  + 3 \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_KroneckerDelta():
    expr = KroneckerDelta(x, y)
    ascii_str = \
        """\
d   \n\
 x,y\
"""
    ucode_str = \
        """\
δ   \n\
 x,y\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_product():
    f = symbols('f', cls=Function)
    expr = Product(f((n/3)**2), (n, k**2, l))

    unicode_str = \
        """\
    l           \n\
┬────────┬      \n\
│        │  ⎛ 2⎞\n\
│        │  ⎜n ⎟\n\
│        │ f⎜──⎟\n\
│        │  ⎝9 ⎠\n\
│        │      \n\
       2        \n\
  n = k         """
    ascii_str = \
        """\
    l           \n\
__________      \n\
|        |  / 2\\\n\
|        |  |n |\n\
|        | f|--|\n\
|        |  \\9 /\n\
|        |      \n\
       2        \n\
  n = k         """

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str

    expr = Product(f((n/3)**2), (n, k**2, l), (l, 1, m))

    unicode_str = \
        """\
    m          l           \n\
┬────────┬ ┬────────┬      \n\
│        │ │        │  ⎛ 2⎞\n\
│        │ │        │  ⎜n ⎟\n\
│        │ │        │ f⎜──⎟\n\
│        │ │        │  ⎝9 ⎠\n\
│        │ │        │      \n\
  l = 1           2        \n\
             n = k         """
    ascii_str = \
        """\
    m          l           \n\
__________ __________      \n\
|        | |        |  / 2\\\n\
|        | |        |  |n |\n\
|        | |        | f|--|\n\
|        | |        |  \\9 /\n\
|        | |        |      \n\
  l = 1           2        \n\
             n = k         """

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str


def test_pretty_lambda():
    # S.IdentityFunction is a special case
    expr = Lambda(y, y)
    assert pretty(expr) == 'dummy_for_IdentityFunction -> dummy_for_IdentityFunction'
    assert upretty(expr) == 'dummy_for_IdentityFunction ↦ dummy_for_IdentityFunction'

    expr = Lambda(x, x+1)
    assert pretty(expr) == 'x -> x + 1'
    assert upretty(expr) == 'x ↦ x + 1'

    expr = Lambda(x, x**2)
    ascii_str = \
        """\
      2\n\
x -> x \
"""
    ucode_str = \
        """\
     2\n\
x ↦ x \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Lambda(x, x**2)**2
    ascii_str = \
        """\
         2
/      2\\ \n\
\\x -> x / \
"""
    ucode_str = \
        """\
        2
⎛     2⎞ \n\
⎝x ↦ x ⎠ \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Lambda((x, y), x)
    ascii_str = '(x, y) -> x'
    ucode_str = '(x, y) ↦ x'
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Lambda((x, y), x**2)
    ascii_str = \
        """\
           2\n\
(x, y) -> x \
"""
    ucode_str = \
        """\
          2\n\
(x, y) ↦ x \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_order():
    expr = O(1)
    ascii_str = \
        """\
O(1)\
"""
    ucode_str = \
        """\
O(1)\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(1/x)
    ascii_str = \
        """\
 /1\\\n\
O|-|\n\
 \\x/\
"""
    ucode_str = \
        """\
 ⎛1⎞\n\
O⎜─⎟\n\
 ⎝x⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(x**2 + y**2)
    ascii_str = \
        """\
 / 2    2                  \\\n\
O\\x  + y ; (x, y) -> (0, 0)/\
"""
    ucode_str = \
        """\
 ⎛ 2    2                 ⎞\n\
O⎝x  + y ; (x, y) → (0, 0)⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(1, (x, oo))
    ascii_str = \
        """\
O(1; x -> oo)\
"""
    ucode_str = \
        """\
O(1; x → ∞)\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(1/x, (x, oo))
    ascii_str = \
        """\
 /1         \\\n\
O|-; x -> oo|\n\
 \\x         /\
"""
    ucode_str = \
        """\
 ⎛1       ⎞\n\
O⎜─; x → ∞⎟\n\
 ⎝x       ⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = O(x**2 + y**2, (x, oo), (y, oo))
    ascii_str = \
        """\
 / 2    2                    \\\n\
O\\x  + y ; (x, y) -> (oo, oo)/\
"""
    ucode_str = \
        """\
 ⎛ 2    2                 ⎞\n\
O⎝x  + y ; (x, y) → (∞, ∞)⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_derivatives():
    # Simple
    expr = Derivative(log(x), x, evaluate=False)
    ascii_str = \
        """\
d         \n\
--(log(x))\n\
dx        \
"""
    ucode_str = \
        """\
d         \n\
──(log(x))\n\
dx        \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Derivative(log(x), x, evaluate=False) + x
    ascii_str_1 = \
        """\
    d         \n\
x + --(log(x))\n\
    dx        \
"""
    ascii_str_2 = \
        """\
d             \n\
--(log(x)) + x\n\
dx            \
"""
    ucode_str_1 = \
        """\
    d         \n\
x + ──(log(x))\n\
    dx        \
"""
    ucode_str_2 = \
        """\
d             \n\
──(log(x)) + x\n\
dx            \
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    # basic partial derivatives
    y = Symbol('y')
    expr = Derivative(log(x + y) + x, x)
    ascii_str_1 = \
        """\
d                 \n\
--(log(x + y) + x)\n\
dx                \
"""
    ascii_str_2 = \
        """\
d                 \n\
--(x + log(x + y))\n\
dx                \
"""
    ucode_str_1 = \
        """\
∂                 \n\
──(log(x + y) + x)\n\
∂x                \
"""
    ucode_str_2 = \
        """\
∂                 \n\
──(x + log(x + y))\n\
∂x                \
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2], upretty(expr)

    # Multiple symbols
    expr = Derivative(log(x) + x**2, x, y)
    ascii_str_1 = \
        """\
   2              \n\
  d  /          2\\\n\
-----\\log(x) + x /\n\
dy dx             \
"""
    ascii_str_2 = \
        """\
   2              \n\
  d  / 2         \\\n\
-----\\x  + log(x)/\n\
dy dx             \
"""
    ucode_str_1 = \
        """\
   2              \n\
  d  ⎛          2⎞\n\
─────⎝log(x) + x ⎠\n\
dy dx             \
"""
    ucode_str_2 = \
        """\
   2              \n\
  d  ⎛ 2         ⎞\n\
─────⎝x  + log(x)⎠\n\
dy dx             \
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = Derivative(2*x*y, y, x) + x**2
    ascii_str_1 = \
        """\
   2             \n\
  d             2\n\
-----(2*x*y) + x \n\
dx dy            \
"""
    ascii_str_2 = \
        """\
        2        \n\
 2     d         \n\
x  + -----(2*x*y)\n\
     dx dy       \
"""
    ucode_str_1 = \
        """\
   2             \n\
  ∂             2\n\
─────(2⋅x⋅y) + x \n\
∂x ∂y            \
"""
    ucode_str_2 = \
        """\
        2        \n\
 2     ∂         \n\
x  + ─────(2⋅x⋅y)\n\
     ∂x ∂y       \
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = Derivative(2*x*y, x, x)
    ascii_str = \
        """\
  2       \n\
 d        \n\
---(2*x*y)\n\
  2       \n\
dx        \
"""
    ucode_str = \
        """\
  2       \n\
 ∂        \n\
───(2⋅x⋅y)\n\
  2       \n\
∂x        \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Derivative(2*x*y, (x, 17))
    ascii_str = \
        """\
 17        \n\
d          \n\
----(2*x*y)\n\
  17       \n\
dx         \
"""
    ucode_str = \
        """\
 17        \n\
∂          \n\
────(2⋅x⋅y)\n\
  17       \n\
∂x         \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Derivative(2*x*y, x, x, y)
    ascii_str = \
        """\
   3         \n\
  d          \n\
------(2*x*y)\n\
     2       \n\
dy dx        \
"""
    ucode_str = \
        """\
   3         \n\
  ∂          \n\
──────(2⋅x⋅y)\n\
     2       \n\
∂y ∂x        \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # Greek letters
    alpha = Symbol('alpha')
    beta = Function('beta')
    expr = beta(alpha).diff(alpha)
    ascii_str = \
        """\
  d                \n\
------(beta(alpha))\n\
dalpha             \
"""
    ucode_str = \
        """\
d       \n\
──(β(α))\n\
dα      \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # sympy/sympy#4335
    y = Function('y')
    expr = -y(x).diff(x)
    ucode_str = \
        """\
 d       \n\
-──(y(x))\n\
 dx      \
"""
    ascii_str = \
        """\
  d       \n\
- --(y(x))\n\
  dx      \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_integrals():
    expr = Integral(log(x), x)
    ascii_str = \
        """\
  /         \n\
 |          \n\
 | log(x) dx\n\
 |          \n\
/           \
"""
    ucode_str = \
        """\
⌠          \n\
⎮ log(x) dx\n\
⌡          \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**2, x)
    ascii_str = \
        """\
  /     \n\
 |      \n\
 |  2   \n\
 | x  dx\n\
 |      \n\
/       \
"""
    ucode_str = \
        """\
⌠      \n\
⎮  2   \n\
⎮ x  dx\n\
⌡      \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral((sin(x))**2 / (tan(x))**2)
    ascii_str = \
        """\
  /          \n\
 |           \n\
 |    2      \n\
 | sin (x)   \n\
 | ------- dx\n\
 |    2      \n\
 | tan (x)   \n\
 |           \n\
/            \
"""
    ucode_str = \
        """\
⌠           \n\
⎮    2      \n\
⎮ sin (x)   \n\
⎮ ─────── dx\n\
⎮    2      \n\
⎮ tan (x)   \n\
⌡           \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**(2**x), x)
    ascii_str = \
        """\
  /        \n\
 |         \n\
 |  / x\\   \n\
 |  \\2 /   \n\
 | x     dx\n\
 |         \n\
/          \
"""
    ucode_str = \
        """\
⌠         \n\
⎮  ⎛ x⎞   \n\
⎮  ⎝2 ⎠   \n\
⎮ x     dx\n\
⌡         \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**2, (x, 1, 2))
    ascii_str = \
        """\
  2      \n\
  /      \n\
 |       \n\
 |   2   \n\
 |  x  dx\n\
 |       \n\
/        \n\
1        \
"""
    ucode_str = \
        """\
2      \n\
⌠      \n\
⎮  2   \n\
⎮ x  dx\n\
⌡      \n\
1      \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**2, (x, Rational(1, 2), 10))
    ascii_str = \
        """\
 10      \n\
  /      \n\
 |       \n\
 |   2   \n\
 |  x  dx\n\
 |       \n\
/        \n\
1/2      \
"""
    ucode_str = \
        """\
 10      \n\
 ⌠       \n\
 ⎮   2   \n\
 ⎮  x  dx\n\
 ⌡       \n\
1/2      \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(x**2*y**2, x, y)
    ascii_str = \
        """\
  /  /           \n\
 |  |            \n\
 |  |  2  2      \n\
 |  | x *y  dx dy\n\
 |  |            \n\
/  /             \
"""
    ucode_str = \
        """\
⌠ ⌠            \n\
⎮ ⎮  2  2      \n\
⎮ ⎮ x ⋅y  dx dy\n\
⌡ ⌡            \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(sin(theta)/cos(phi), (theta, 0, pi), (phi, 0, 2*pi))
    ascii_str = \
        """\
 2*pi pi                           \n\
   /   /                           \n\
  |   |                            \n\
  |   |  sin(theta)                \n\
  |   |  ---------- d(theta) d(phi)\n\
  |   |   cos(phi)                 \n\
  |   |                            \n\
 /   /                             \n\
 0   0                             \
"""
    ucode_str = \
        """\
2⋅π π             \n\
 ⌠  ⌠             \n\
 ⎮  ⎮ sin(θ)      \n\
 ⎮  ⎮ ────── dθ dφ\n\
 ⎮  ⎮ cos(φ)      \n\
 ⌡  ⌡             \n\
 0  0             \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(sin(x) + log(x), x)
    ascii_str = \
        """\
  /                    \n\
 |                     \n\
 | (log(x) + sin(x)) dx\n\
 |                     \n\
/                      \
"""
    ucode_str = \
        """\
⌠                     \n\
⎮ (log(x) + sin(x)) dx\n\
⌡                     \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(sin(x), (x, x))
    ascii_str = \
        """\
  x          \n\
  /          \n\
 |           \n\
 |  sin(x) dx\n\
 |           \n\
/            \n\
             \
"""
    ucode_str = \
        """\
x          \n\
⌠          \n\
⎮ sin(x) dx\n\
⌡          \n\
           \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral((x**4 + x**2*exp(x) - x**2 - 2*x*exp(x) - 2*x -
                     exp(x))*exp(x)/((x - 1)**2*(x + 1)**2*(exp(x) + 1)))

    ucode_str = \
        """\
⌠                                            \n\
⎮  x ⎛ x  2      x      x    4    2      ⎞   \n\
⎮ ℯ ⋅⎝ℯ ⋅x  - 2⋅ℯ ⋅x - ℯ  + x  - x  - 2⋅x⎠   \n\
⎮ ──────────────────────────────────────── dx\n\
⎮        ⎛ x    ⎞        2        2          \n\
⎮        ⎝ℯ  + 1⎠⋅(x - 1) ⋅(x + 1)           \n\
⌡                                            \
"""
    assert upretty(expr) == ucode_str


def test_pretty_matrix():
    # Empty Matrix
    expr = Matrix()
    ascii_str = '[]'
    unicode_str = '[]'
    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str
    expr = Matrix(2, 0, lambda i, j: 0)
    ascii_str = '[]'
    unicode_str = '[]'
    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str
    expr = Matrix(0, 2, lambda i, j: 0)
    ascii_str = '[]'
    unicode_str = '[]'
    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str
    expr = Matrix([[x**2 + 1, 1], [y, x + y]])
    ascii_str_1 = \
        """\
[     2       ]
[1 + x     1  ]
[             ]
[  y     x + y]\
"""
    ascii_str_2 = \
        """\
[ 2           ]
[x  + 1    1  ]
[             ]
[  y     x + y]\
"""
    ucode_str_1 = \
        """\
⎡     2       ⎤
⎢1 + x     1  ⎥
⎢             ⎥
⎣  y     x + y⎦\
"""
    ucode_str_2 = \
        """\
⎡ 2           ⎤
⎢x  + 1    1  ⎥
⎢             ⎥
⎣  y     x + y⎦\
"""
    assert pretty(expr) in [ascii_str_1, ascii_str_2]
    assert upretty(expr) in [ucode_str_1, ucode_str_2]

    expr = Matrix([[x/y, y, theta], [0, exp(I*k*phi), 1]])
    ascii_str = \
        """\
[x                 ]
[-     y      theta]
[y                 ]
[                  ]
[    I*k*phi       ]
[0  E           1  ]\
"""
    ucode_str = \
        """\
⎡x           ⎤
⎢─    y     θ⎥
⎢y           ⎥
⎢            ⎥
⎢    ⅈ⋅k⋅φ   ⎥
⎣0  ℯ       1⎦\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    A = MatrixSymbol('A', 3, 3)
    B = MatrixSymbol('B', 3, 3)
    C = MatrixSymbol('C', 3, 3)
    expr = A*(B + C)
    assert upretty(expr) == 'A⋅(B + C)'


def test_pretty_ndim_arrays():
    for ArrayType in (ImmutableDenseNDimArray, ImmutableSparseNDimArray,
                      MutableDenseNDimArray, MutableSparseNDimArray):
        M0 = ArrayType((x,), ())

        assert pretty(M0) == 'x'
        assert upretty(M0) == 'x'

        M = ArrayType([[1/x, y], [z, w]])
        M1 = ArrayType([1/x, y, z])

        M2 = tensorproduct(M1, M)
        M3 = tensorproduct(M, M)

        ascii_str = \
            """\
[1   ]\n\
[-  y]\n\
[x   ]\n\
[    ]\n\
[z  w]\
"""
        ucode_str = \
            """\
⎡1   ⎤\n\
⎢─  y⎥\n\
⎢x   ⎥\n\
⎢    ⎥\n\
⎣z  w⎦\
"""
        assert pretty(M) == ascii_str
        assert upretty(M) == ucode_str

        ascii_str = \
            """\
[1      ]\n\
[-  y  z]\n\
[x      ]\
"""
        ucode_str = \
            """\
⎡1      ⎤\n\
⎢─  y  z⎥\n\
⎣x      ⎦\
"""
        assert pretty(M1) == ascii_str
        assert upretty(M1) == ucode_str

        ascii_str = \
            """\
[[1   y]                       ]\n\
[[--  -]              [z      ]]\n\
[[ 2  x]  [ y    2 ]  [-   y*z]]\n\
[[x    ]  [ -   y  ]  [x      ]]\n\
[[     ]  [ x      ]  [       ]]\n\
[[z   w]  [        ]  [ 2     ]]\n\
[[-   -]  [y*z  w*y]  [z   w*z]]\n\
[[x   x]                       ]\
"""
        ucode_str = \
            """\
⎡⎡1   y⎤                       ⎤\n\
⎢⎢──  ─⎥              ⎡z      ⎤⎥\n\
⎢⎢ 2  x⎥  ⎡ y    2 ⎤  ⎢─   y⋅z⎥⎥\n\
⎢⎢x    ⎥  ⎢ ─   y  ⎥  ⎢x      ⎥⎥\n\
⎢⎢     ⎥  ⎢ x      ⎥  ⎢       ⎥⎥\n\
⎢⎢z   w⎥  ⎢        ⎥  ⎢ 2     ⎥⎥\n\
⎢⎢─   ─⎥  ⎣y⋅z  w⋅y⎦  ⎣z   w⋅z⎦⎥\n\
⎣⎣x   x⎦                       ⎦\
"""
        assert pretty(M2) == ascii_str
        assert upretty(M2) == ucode_str

        ascii_str = \
            """\
[ [1   y]             ]\n\
[ [--  -]             ]\n\
[ [ 2  x]   [ y    2 ]]\n\
[ [x    ]   [ -   y  ]]\n\
[ [     ]   [ x      ]]\n\
[ [z   w]   [        ]]\n\
[ [-   -]   [y*z  w*y]]\n\
[ [x   x]             ]\n\
[                     ]\n\
[[z      ]  [ w      ]]\n\
[[-   y*z]  [ -   w*y]]\n\
[[x      ]  [ x      ]]\n\
[[       ]  [        ]]\n\
[[ 2     ]  [      2 ]]\n\
[[z   w*z]  [w*z  w  ]]\
"""
        ucode_str = \
            """\
⎡ ⎡1   y⎤             ⎤\n\
⎢ ⎢──  ─⎥             ⎥\n\
⎢ ⎢ 2  x⎥   ⎡ y    2 ⎤⎥\n\
⎢ ⎢x    ⎥   ⎢ ─   y  ⎥⎥\n\
⎢ ⎢     ⎥   ⎢ x      ⎥⎥\n\
⎢ ⎢z   w⎥   ⎢        ⎥⎥\n\
⎢ ⎢─   ─⎥   ⎣y⋅z  w⋅y⎦⎥\n\
⎢ ⎣x   x⎦             ⎥\n\
⎢                     ⎥\n\
⎢⎡z      ⎤  ⎡ w      ⎤⎥\n\
⎢⎢─   y⋅z⎥  ⎢ ─   w⋅y⎥⎥\n\
⎢⎢x      ⎥  ⎢ x      ⎥⎥\n\
⎢⎢       ⎥  ⎢        ⎥⎥\n\
⎢⎢ 2     ⎥  ⎢      2 ⎥⎥\n\
⎣⎣z   w⋅z⎦  ⎣w⋅z  w  ⎦⎦\
"""
        assert pretty(M3) == ascii_str
        assert upretty(M3) == ucode_str

        Mrow = ArrayType([[x, y, 1 / z]])
        Mcolumn = ArrayType([[x], [y], [1 / z]])
        Mcol2 = ArrayType([Mcolumn.tolist()])

        ascii_str = \
            """\
[[      1]]\n\
[[x  y  -]]\n\
[[      z]]\
"""
        ucode_str = \
            """\
⎡⎡      1⎤⎤\n\
⎢⎢x  y  ─⎥⎥\n\
⎣⎣      z⎦⎦\
"""
        assert pretty(Mrow) == ascii_str
        assert upretty(Mrow) == ucode_str

        ascii_str = \
            """\
[x]\n\
[ ]\n\
[y]\n\
[ ]\n\
[1]\n\
[-]\n\
[z]\
"""
        ucode_str = \
            """\
⎡x⎤\n\
⎢ ⎥\n\
⎢y⎥\n\
⎢ ⎥\n\
⎢1⎥\n\
⎢─⎥\n\
⎣z⎦\
"""
        assert pretty(Mcolumn) == ascii_str
        assert upretty(Mcolumn) == ucode_str

        ascii_str = \
            """\
[[x]]\n\
[[ ]]\n\
[[y]]\n\
[[ ]]\n\
[[1]]\n\
[[-]]\n\
[[z]]\
"""
        ucode_str = \
            """\
⎡⎡x⎤⎤\n\
⎢⎢ ⎥⎥\n\
⎢⎢y⎥⎥\n\
⎢⎢ ⎥⎥\n\
⎢⎢1⎥⎥\n\
⎢⎢─⎥⎥\n\
⎣⎣z⎦⎦\
"""
        assert pretty(Mcol2) == ascii_str
        assert upretty(Mcol2) == ucode_str


def test_Adjoint():
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    assert pretty(Adjoint(X)) == ' +\nX '
    assert pretty(Adjoint(X + Y)) == '       +\n(X + Y) '
    assert pretty(Adjoint(X) + Adjoint(Y)) == ' +    +\nX  + Y '
    assert pretty(Adjoint(X*Y)) == '     +\n(X*Y) '
    assert pretty(Adjoint(Y)*Adjoint(X)) == ' +  +\nY *X '
    assert pretty(Adjoint(X**2)) == '    +\n/ 2\\ \n\\X / '
    assert pretty(Adjoint(X)**2) == '    2\n/ +\\ \n\\X / '
    assert pretty(Adjoint(Inverse(X))) == '     +\n/ -1\\ \n\\X  / '
    assert pretty(Inverse(Adjoint(X))) == '    -1\n/ +\\  \n\\X /  '
    assert pretty(Adjoint(Transpose(X))) == '    +\n/ T\\ \n\\X / '
    assert pretty(Transpose(Adjoint(X))) == '    T\n/ +\\ \n\\X / '
    assert upretty(Adjoint(X)) == ' †\nX '
    assert upretty(Adjoint(X + Y)) == '       †\n(X + Y) '
    assert upretty(Adjoint(X) + Adjoint(Y)) == ' †    †\nX  + Y '
    assert upretty(Adjoint(X*Y)) == '     †\n(X⋅Y) '
    assert upretty(Adjoint(Y)*Adjoint(X)) == ' †  †\nY ⋅X '
    assert upretty(Adjoint(X**2)) == \
        '    †\n⎛ 2⎞ \n⎝X ⎠ '
    assert upretty(Adjoint(X)**2) == \
        '    2\n⎛ †⎞ \n⎝X ⎠ '
    assert upretty(Adjoint(Inverse(X))) == \
        '     †\n⎛ -1⎞ \n⎝X  ⎠ '
    assert upretty(Inverse(Adjoint(X))) == \
        '    -1\n⎛ †⎞  \n⎝X ⎠  '
    assert upretty(Adjoint(Transpose(X))) == \
        '    †\n⎛ T⎞ \n⎝X ⎠ '
    assert upretty(Transpose(Adjoint(X))) == \
        '    T\n⎛ †⎞ \n⎝X ⎠ '


def test_pretty_Trace():
    # issue sympy/sympy#9044
    X = Matrix([[1, 2], [3, 4]])
    Y = Matrix([[2, 4], [6, 8]])
    ascii_str_1 = \
        """\
  /[1  2]\\
tr|[    ]|
  \\[3  4]/\
"""
    ucode_str_1 = \
        """\
  ⎛⎡1  2⎤⎞
tr⎜⎢    ⎥⎟
  ⎝⎣3  4⎦⎠\
"""
    ascii_str_2 = \
        """\
  /[1  2]\\     /[2  4]\\
tr|[    ]| + tr|[    ]|
  \\[3  4]/     \\[6  8]/\
"""
    ucode_str_2 = \
        """\
  ⎛⎡1  2⎤⎞     ⎛⎡2  4⎤⎞
tr⎜⎢    ⎥⎟ + tr⎜⎢    ⎥⎟
  ⎝⎣3  4⎦⎠     ⎝⎣6  8⎦⎠\
"""
    assert pretty(Trace(X)) == ascii_str_1
    assert upretty(Trace(X)) == ucode_str_1

    assert pretty(Trace(X) + Trace(Y)) == ascii_str_2
    assert upretty(Trace(X) + Trace(Y)) == ucode_str_2


def test_MatrixExpressions():
    n = Symbol('n', integer=True)
    X = MatrixSymbol('X', n, n)

    assert pretty(X) == upretty(X) == 'X'

    Y = X[1:2:3, 4:5:6]

    ascii_str = ucode_str = 'X[1:3, 4:6]'

    assert pretty(Y) == ascii_str
    assert upretty(Y) == ucode_str

    Z = X[1:10:2]

    ascii_str = ucode_str = 'X[1:10:2, :n]'

    assert pretty(Z) == ascii_str
    assert upretty(Z) == ucode_str


def test_pretty_piecewise():
    expr = Piecewise((x, x < 1), (x**2, True))
    ascii_str = \
        """\
/x   for x < 1\n\
|             \n\
< 2           \n\
|x   otherwise\n\
\\             \
"""
    ucode_str = \
        """\
⎧x   for x < 1\n\
⎪             \n\
⎨ 2           \n\
⎪x   otherwise\n\
⎩             \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = -Piecewise((x, x < 1), (x**2, True))
    ascii_str = \
        """\
 //x   for x < 1\\\n\
 ||             |\n\
-|< 2           |\n\
 ||x   otherwise|\n\
 \\\\             /\
"""
    ucode_str = \
        """\
 ⎛⎧x   for x < 1⎞\n\
 ⎜⎪             ⎟\n\
-⎜⎨ 2           ⎟\n\
 ⎜⎪x   otherwise⎟\n\
 ⎝⎩             ⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x + Piecewise((x, x > 0), (y, True)) + Piecewise((x/y, x < 2),
                                                            (y**2, x > 2), (1, True)) + 1
    ascii_str = \
        """\
                      //x            \\    \n\
                      ||-   for x < 2|    \n\
                      ||y            |    \n\
    //x  for x > 0\\   ||             |    \n\
x + |<            | + |< 2           | + 1\n\
    \\\\y  otherwise/   ||y   for x > 2|    \n\
                      ||             |    \n\
                      ||1   otherwise|    \n\
                      \\\\             /    \
"""
    ucode_str = \
        """\
                      ⎛⎧x            ⎞    \n\
                      ⎜⎪─   for x < 2⎟    \n\
                      ⎜⎪y            ⎟    \n\
    ⎛⎧x  for x > 0⎞   ⎜⎪             ⎟    \n\
x + ⎜⎨            ⎟ + ⎜⎨ 2           ⎟ + 1\n\
    ⎝⎩y  otherwise⎠   ⎜⎪y   for x > 2⎟    \n\
                      ⎜⎪             ⎟    \n\
                      ⎜⎪1   otherwise⎟    \n\
                      ⎝⎩             ⎠    \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x - Piecewise((x, x > 0), (y, True)) + Piecewise((x/y, x < 2),
                                                            (y**2, x > 2), (1, True)) + 1
    ascii_str = \
        """\
                      //x            \\    \n\
                      ||-   for x < 2|    \n\
                      ||y            |    \n\
    //x  for x > 0\\   ||             |    \n\
x - |<            | + |< 2           | + 1\n\
    \\\\y  otherwise/   ||y   for x > 2|    \n\
                      ||             |    \n\
                      ||1   otherwise|    \n\
                      \\\\             /    \
"""
    ucode_str = \
        """\
                      ⎛⎧x            ⎞    \n\
                      ⎜⎪─   for x < 2⎟    \n\
                      ⎜⎪y            ⎟    \n\
    ⎛⎧x  for x > 0⎞   ⎜⎪             ⎟    \n\
x - ⎜⎨            ⎟ + ⎜⎨ 2           ⎟ + 1\n\
    ⎝⎩y  otherwise⎠   ⎜⎪y   for x > 2⎟    \n\
                      ⎜⎪             ⎟    \n\
                      ⎜⎪1   otherwise⎟    \n\
                      ⎝⎩             ⎠    \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x*Piecewise((x, x > 0), (y, True))
    ascii_str = \
        """\
  //x  for x > 0\\\n\
x*|<            |\n\
  \\\\y  otherwise/\
"""
    ucode_str = \
        """\
  ⎛⎧x  for x > 0⎞\n\
x⋅⎜⎨            ⎟\n\
  ⎝⎩y  otherwise⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Piecewise((x, x > 0), (y, True))*Piecewise((x/y, x < 2), (y**2, x >
                                                                     2), (1, True))
    ascii_str = \
        """\
                //x            \\\n\
                ||-   for x < 2|\n\
                ||y            |\n\
//x  for x > 0\\ ||             |\n\
|<            |*|< 2           |\n\
\\\\y  otherwise/ ||y   for x > 2|\n\
                ||             |\n\
                ||1   otherwise|\n\
                \\\\             /\
"""
    ucode_str = \
        """\
                ⎛⎧x            ⎞\n\
                ⎜⎪─   for x < 2⎟\n\
                ⎜⎪y            ⎟\n\
⎛⎧x  for x > 0⎞ ⎜⎪             ⎟\n\
⎜⎨            ⎟⋅⎜⎨ 2           ⎟\n\
⎝⎩y  otherwise⎠ ⎜⎪y   for x > 2⎟\n\
                ⎜⎪             ⎟\n\
                ⎜⎪1   otherwise⎟\n\
                ⎝⎩             ⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = -Piecewise((x, x > 0), (y, True))*Piecewise((x/y, x < 2), (y**2, x
                                                                      > 2), (1, True))
    ascii_str = \
        """\
                 //x            \\\n\
                 ||-   for x < 2|\n\
                 ||y            |\n\
 //x  for x > 0\\ ||             |\n\
-|<            |*|< 2           |\n\
 \\\\y  otherwise/ ||y   for x > 2|\n\
                 ||             |\n\
                 ||1   otherwise|\n\
                 \\\\             /\
"""
    ucode_str = \
        """\
                 ⎛⎧x            ⎞\n\
                 ⎜⎪─   for x < 2⎟\n\
                 ⎜⎪y            ⎟\n\
 ⎛⎧x  for x > 0⎞ ⎜⎪             ⎟\n\
-⎜⎨            ⎟⋅⎜⎨ 2           ⎟\n\
 ⎝⎩y  otherwise⎠ ⎜⎪y   for x > 2⎟\n\
                 ⎜⎪             ⎟\n\
                 ⎜⎪1   otherwise⎟\n\
                 ⎝⎩             ⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Piecewise((0, abs(1/y) < 1), (1, abs(y) < 1), (y*meijerg(((2, 1),
                                                                     ()), ((), (1, 0)), 1/y), True))
    ascii_str = \
        """\
/                                |1|    \n\
|            0               for |-| < 1\n\
|                                |y|    \n\
|                                       \n\
<            1               for |y| < 1\n\
|                                       \n\
|   __0, 2 /2, 1       | 1\\             \n\
|y*/__     |           | -|   otherwise \n\
\\  \\_|2, 2 \\      1, 0 | y/             \
"""
    ucode_str = \
        """\
⎧                                │1│    \n\
⎪            0               for │─│ < 1\n\
⎪                                │y│    \n\
⎪                                       \n\
⎨            1               for │y│ < 1\n\
⎪                                       \n\
⎪  ╭─╮0, 2 ⎛2, 1       │ 1⎞             \n\
⎪y⋅│╶┐     ⎜           │ ─⎟   otherwise \n\
⎩  ╰─╯2, 2 ⎝      1, 0 │ y⎠             \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # XXX: We have to use evaluate=False here because Piecewise._eval_power
    # denests the power.
    expr = Pow(Piecewise((x, x > 0), (y, True)), 2, evaluate=False)
    ascii_str = \
        """\
               2\n\
//x  for x > 0\\ \n\
|<            | \n\
\\\\y  otherwise/ \
"""
    ucode_str = \
        """\
               2\n\
⎛⎧x  for x > 0⎞ \n\
⎜⎨            ⎟ \n\
⎝⎩y  otherwise⎠ \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_seq():
    expr = ()
    ascii_str = \
        """\
()\
"""
    ucode_str = \
        """\
()\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = []
    ascii_str = \
        """\
[]\
"""
    ucode_str = \
        """\
[]\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = {}
    expr_2 = {}
    ascii_str = \
        """\
{}\
"""
    ucode_str = \
        """\
{}\
"""
    assert pretty(expr) == ascii_str
    assert pretty(expr_2) == ascii_str
    assert upretty(expr) == ucode_str
    assert upretty(expr_2) == ucode_str

    expr = 1/x,
    ascii_str = \
        """\
 1  \n\
(-,)\n\
 x  \
"""
    ucode_str = \
        """\
⎛1 ⎞\n\
⎜─,⎟\n\
⎝x ⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = [x**2, 1/x, x, y, sin(theta)**2/cos(phi)**2]
    ascii_str = \
        """\
                 2        \n\
  2  1        sin (theta) \n\
[x , -, x, y, -----------]\n\
     x            2       \n\
               cos (phi)  \
"""
    ucode_str = \
        """\
⎡                2   ⎤\n\
⎢ 2  1        sin (θ)⎥\n\
⎢x , ─, x, y, ───────⎥\n\
⎢    x           2   ⎥\n\
⎣             cos (φ)⎦\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = (x**2, 1/x, x, y, sin(theta)**2/cos(phi)**2)
    ascii_str = \
        """\
                 2        \n\
  2  1        sin (theta) \n\
(x , -, x, y, -----------)\n\
     x            2       \n\
               cos (phi)  \
"""
    ucode_str = \
        """\
⎛                2   ⎞\n\
⎜ 2  1        sin (θ)⎟\n\
⎜x , ─, x, y, ───────⎟\n\
⎜    x           2   ⎟\n\
⎝             cos (φ)⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Tuple(x**2, 1/x, x, y, sin(theta)**2/cos(phi)**2)
    ascii_str = \
        """\
                 2        \n\
  2  1        sin (theta) \n\
(x , -, x, y, -----------)\n\
     x            2       \n\
               cos (phi)  \
"""
    ucode_str = \
        """\
⎛                2   ⎞\n\
⎜ 2  1        sin (θ)⎟\n\
⎜x , ─, x, y, ───────⎟\n\
⎜    x           2   ⎟\n\
⎝             cos (φ)⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = {x: sin(x)}
    expr_2 = Dict({x: sin(x)})
    ascii_str = \
        """\
{x: sin(x)}\
"""
    ucode_str = \
        """\
{x: sin(x)}\
"""
    assert pretty(expr) == ascii_str
    assert pretty(expr_2) == ascii_str
    assert upretty(expr) == ucode_str
    assert upretty(expr_2) == ucode_str

    expr = {1/x: 1/y, x: sin(x)**2}
    expr_2 = Dict({1/x: 1/y, x: sin(x)**2})
    ascii_str = \
        """\
 1  1        2    \n\
{-: -, x: sin (x)}\n\
 x  y             \
"""
    ucode_str = \
        """\
⎧1  1        2   ⎫\n\
⎨─: ─, x: sin (x)⎬\n\
⎩x  y            ⎭\
"""
    assert pretty(expr) == ascii_str
    assert pretty(expr_2) == ascii_str
    assert upretty(expr) == ucode_str
    assert upretty(expr_2) == ucode_str

    # There used to be a bug with pretty-printing sequences of even height.
    expr = [x**2]
    ascii_str = \
        """\
  2 \n\
[x ]\
"""
    ucode_str = \
        """\
⎡ 2⎤\n\
⎣x ⎦\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = x**2,
    ascii_str = \
        """\
  2  \n\
(x ,)\
"""
    ucode_str = \
        """\
⎛ 2 ⎞\n\
⎝x ,⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Tuple(x**2)
    ascii_str = \
        """\
  2  \n\
(x ,)\
"""
    ucode_str = \
        """\
⎛ 2 ⎞\n\
⎝x ,⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = {x**2: 1}
    expr_2 = Dict({x**2: 1})
    ascii_str = \
        """\
  2    \n\
{x : 1}\
"""
    ucode_str = \
        """\
⎧ 2   ⎫\n\
⎨x : 1⎬\n\
⎩     ⎭\
"""
    assert pretty(expr) == ascii_str
    assert pretty(expr_2) == ascii_str
    assert upretty(expr) == ucode_str
    assert upretty(expr_2) == ucode_str


def test_any_object_in_sequence():
    # Cf. issue sympy/sympy#5306
    b1 = Basic()
    b2 = Basic(Basic())

    expr = [b2, b1]
    assert pretty(expr) == '[Basic(Basic()), Basic()]'
    assert upretty(expr) == '[Basic(Basic()), Basic()]'

    expr = {b2, b1}
    assert pretty(expr) == '{Basic(), Basic(Basic())}'
    assert upretty(expr) == '{Basic(), Basic(Basic())}'

    expr = {b2: b1, b1: b2}
    expr2 = Dict({b2: b1, b1: b2})
    assert pretty(expr) == '{Basic(): Basic(Basic()), Basic(Basic()): Basic()}'
    assert pretty(
        expr2) == '{Basic(): Basic(Basic()), Basic(Basic()): Basic()}'
    assert upretty(
        expr) == '{Basic(): Basic(Basic()), Basic(Basic()): Basic()}'
    assert upretty(
        expr2) == '{Basic(): Basic(Basic()), Basic(Basic()): Basic()}'


def test_pretty_sets():
    assert pretty(FiniteSet(x*y, x**2)) == \
        """\
  2      \n\
{x , x*y}\
"""
    assert pretty(FiniteSet(*range(1, 6))) == '{1, 2, 3, 4, 5}'
    assert pretty(FiniteSet(*range(1, 13))) == '{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}'

    assert pretty({x*y, x**2}) == \
        """\
  2      \n\
{x , x*y}\
"""
    assert pretty(set(range(1, 6))) == '{1, 2, 3, 4, 5}'
    assert pretty(set(range(1, 13))) == \
        '{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}'

    assert pretty(frozenset({x*y, x**2})) == \
        """\
            2       \n\
frozenset({x , x*y})\
"""
    assert pretty(frozenset(range(1, 6))) == 'frozenset({1, 2, 3, 4, 5})'
    assert pretty(frozenset(range(1, 13))) == \
        'frozenset({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12})'

    assert pretty(Range(0, 3, 1)) == '{0, 1, 2}'

    ascii_str = '{0, 1, ..., 29}'
    ucode_str = '{0, 1, …, 29}'
    assert pretty(Range(0, 30, 1)) == ascii_str
    assert upretty(Range(0, 30, 1)) == ucode_str

    ascii_str = '{0, 2, ..., oo}'
    ucode_str = '{0, 2, …, ∞}'
    assert pretty(Range(0, oo, 2)) == ascii_str
    assert upretty(Range(0, oo, 2)) == ucode_str

    ascii_str = '{-oo, ..., -3, -2}'
    ucode_str = '{-∞, …, -3, -2}'
    assert pretty(Range(-2, -oo, -1)) == ascii_str
    assert upretty(Range(-2, -oo, -1)) == ucode_str


def test_pretty_Union():
    a, b = Interval(2, 3), Interval(4, 7)
    ucode_str = '[2, 3] ∪ [4, 7]'
    ascii_str = '[2, 3] U [4, 7]'
    assert upretty(Union(a, b)) == ucode_str
    assert pretty(Union(a, b)) == ascii_str

    # issue sympy/sympy#9877
    ucode_str1 = '(2, 3) ∪ ([1, 2] \\ {x})'
    a, b, c = Interval(2, 3, True, True), Interval(1, 2), FiniteSet(x)
    assert upretty(Union(a, Complement(b, c))) == ucode_str1


def test_pretty_Intersection():
    a, b = Interval(x, y), Interval(z, w)
    ucode_str = '[x, y] ∩ [z, w]'
    ascii_str = '[x, y] n [z, w]'
    assert upretty(Intersection(a, b)) == ucode_str
    assert pretty(Intersection(a, b)) == ascii_str

    # issue sympy/sympy#9877
    ucode_str2 = '{x} ∩ {y} ∩ ({z} \\ [1, 2])'
    d, e, f, g = FiniteSet(x), FiniteSet(y), FiniteSet(z), Interval(1, 2)
    assert upretty(Intersection(d, e, Complement(f, g))) == ucode_str2


def test_ProductSet_paranthesis():
    ucode_str = '([4, 7] × {1, 2}) ∪ ([2, 3] × [4, 7])'

    a, b = Interval(2, 3), Interval(4, 7)
    assert upretty(Union(a*b, b*FiniteSet(1, 2))) == ucode_str


def test_sympyissue_10413():
    ascii_str = '[2, 3] x [4, 7]'
    ucode_str = '[2, 3] × [4, 7]'

    a, b = Interval(2, 3), Interval(4, 7)
    assert pretty(a*b) == ascii_str
    assert upretty(a*b) == ucode_str


def test_pretty_limits():
    expr = Limit(x, x, oo)
    ascii_str = \
        """\
 lim x\n\
x->oo \
"""
    ucode_str = \
        """\
lim x\n\
x─→∞ \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(x**2, x, 0)
    ascii_str = \
        """\
      2\n\
 lim x \n\
x->0+  \
"""
    ucode_str = \
        """\
      2\n\
 lim x \n\
x─→0⁺  \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(1/x, x, 0)
    ascii_str = \
        """\
     1\n\
 lim -\n\
x->0+x\
"""
    ucode_str = \
        """\
     1\n\
 lim ─\n\
x─→0⁺x\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(sin(x)/x, x, 0)
    ascii_str = \
        """\
     sin(x)\n\
 lim ------\n\
x->0+  x   \
"""
    ucode_str = \
        """\
     sin(x)\n\
 lim ──────\n\
x─→0⁺  x   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(sin(x)/x, x, 0, '-')
    ascii_str = \
        """\
     sin(x)\n\
 lim ------\n\
x->0-  x   \
"""
    ucode_str = \
        """\
     sin(x)\n\
 lim ──────\n\
x─→0⁻  x   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(sin(x)/x, x, 0, 'real')
    ascii_str = \
        """\
    sin(x)\n\
lim ------\n\
x->0  x   \
"""
    ucode_str = \
        """\
    sin(x)\n\
lim ──────\n\
x─→0  x   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Limit(x - 1/y, x, 1)
    ascii_str = \
        """\
     /    1\\\n\
 lim |x - -|\n\
x->1+\\    y/\
"""
    ucode_str = \
        """\
     ⎛    1⎞\n\
 lim ⎜x - ─⎟\n\
x─→1⁺⎝    y⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    # issue sympy/sympy#10803
    expr = 1/Limit(x, x, oo)
    ascii_str = \
        """\
        -1\n\
/ lim x\\  \n\
\\x->oo /  \
"""
    ucode_str = \
        """\
       -1\n\
⎛lim x⎞  \n\
⎝x─→∞ ⎠  \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_RootOf():
    expr = RootOf(x**5 + 11*x - 2, 0)
    ascii_str = \
        """\
      / 5              \\\n\
RootOf\\x  + 11*x - 2, 0/\
"""
    ucode_str = \
        """\
      ⎛ 5              ⎞\n\
RootOf⎝x  + 11⋅x - 2, 0⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = RootOf(x**5 + y*x - 2, x, 0)
    ascii_str = \
        """\
      / 5                \\\n\
RootOf\\x  + x*y - 2, x, 0/\
"""
    ucode_str = \
        """\
      ⎛ 5                ⎞\n\
RootOf⎝x  + x⋅y - 2, x, 0⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_RootSum():
    expr = RootSum(x**5 + 11*x - 2, auto=False)
    ascii_str = \
        """\
       / 5           \\\n\
RootSum\\x  + 11*x - 2/\
"""
    ucode_str = \
        """\
       ⎛ 5           ⎞\n\
RootSum⎝x  + 11⋅x - 2⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = RootSum(x**5 + 11*x - 2, Lambda(z, exp(z)))
    ascii_str = \
        """\
       / 5                   z\\\n\
RootSum\\x  + 11*x - 2, z -> E /\
"""
    ucode_str = \
        """\
       ⎛ 5                  z⎞\n\
RootSum⎝x  + 11⋅x - 2, z ↦ ℯ ⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_GroebnerBasis():
    expr = groebner([], x, y)

    ascii_str = \
        """\
GroebnerBasis([], x, y, domain=ZZ, order=lex)\
"""
    ucode_str = \
        """\
GroebnerBasis([], x, y, domain=ℤ, order=lex)\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    F = [x**2 - 3*y - x + 1, y**2 - 2*x + y - 1]
    expr = groebner(F, x, y, order='grlex')

    ascii_str = \
        """\
             /[ 2                 2              ]                              \\\n\
GroebnerBasis\\[x  - x - 3*y + 1, y  - 2*x + y - 1], x, y, domain=ZZ, order=grlex/\
"""
    ucode_str = \
        """\
             ⎛⎡ 2                 2              ⎤                             ⎞\n\
GroebnerBasis⎝⎣x  - x - 3⋅y + 1, y  - 2⋅x + y - 1⎦, x, y, domain=ℤ, order=grlex⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = expr.set_order('lex')

    ascii_str = \
        """\
             /[       2           4      3      2           ]                            \\\n\
GroebnerBasis\\[2*x - y  - y + 1, y  + 2*y  - 3*y  - 16*y + 7], x, y, domain=ZZ, order=lex/\
"""
    ucode_str = \
        """\
             ⎛⎡       2           4      3      2           ⎤                           ⎞\n\
GroebnerBasis⎝⎣2⋅x - y  - y + 1, y  + 2⋅y  - 3⋅y  - 16⋅y + 7⎦, x, y, domain=ℤ, order=lex⎠\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_Boolean():
    expr = Not(x, evaluate=False)

    assert pretty(expr) == 'Not(x)'
    assert upretty(expr) == '¬x'

    expr = And(x, y)

    assert pretty(expr) == 'And(x, y)'
    assert upretty(expr) == 'x ∧ y'

    expr = Or(x, y)

    assert pretty(expr) == 'Or(x, y)'
    assert upretty(expr) == 'x ∨ y'

    expr = And(*[a, b, c, d, e, f])

    assert pretty(expr) == 'And(a, b, c, d, e, f)'
    assert upretty(expr) == 'a ∧ b ∧ c ∧ d ∧ e ∧ f'

    expr = Or(*[a, b, c, d, e, f])

    assert pretty(expr) == 'Or(a, b, c, d, e, f)'
    assert upretty(expr) == 'a ∨ b ∨ c ∨ d ∨ e ∨ f'

    expr = Xor(x, y, evaluate=False)

    assert pretty(expr) == 'Xor(x, y)'
    assert upretty(expr) == 'x ⊻ y'

    expr = Nand(x, y, evaluate=False)

    assert pretty(expr) == 'Nand(x, y)'
    assert upretty(expr) == 'x ⊼ y'

    expr = Nor(x, y, evaluate=False)

    assert pretty(expr) == 'Nor(x, y)'
    assert upretty(expr) == 'x ⊽ y'

    expr = Implies(x, y, evaluate=False)

    assert pretty(expr) == 'Implies(x, y)'
    assert upretty(expr) == 'x → y'

    # don't sort args
    expr = Implies(y, x, evaluate=False)

    assert pretty(expr) == 'Implies(y, x)'
    assert upretty(expr) == 'y → x'

    expr = Equivalent(x, y, evaluate=False)

    assert pretty(expr) == 'Equivalent(x, y)'
    assert upretty(expr) == 'x ≡ y'

    expr = Equivalent(y, x, evaluate=False)

    assert pretty(expr) == 'Equivalent(x, y)'
    assert upretty(expr) == 'x ≡ y'

    expr = ~(x & y)
    assert upretty(expr) == '¬(x ∧ y)'

    expr = x & (y | z)
    assert upretty(expr) == 'x ∧ (y ∨ z)'

    expr = (y | z) & (x | z)
    assert upretty(expr) == '(x ∨ z) ∧ (y ∨ z)'

    # issue sympy/sympy#7179
    assert upretty(Not(Equivalent(x, y))) == 'x ≢ y'
    assert upretty(Not(Implies(x, y))) == 'x ↛ y'


def test_pretty_Domain():
    expr = FF(23)

    assert pretty(expr) == 'GF(23)'
    assert upretty(expr) == '𝔽₂₃'

    expr = FF(2, [1, 1, 1])

    assert pretty(expr) == 'GF(4)'
    assert upretty(expr) == '𝔽₄'

    expr = ZZ

    assert pretty(expr) == 'ZZ'
    assert upretty(expr) == 'ℤ'

    expr = QQ

    assert pretty(expr) == 'QQ'
    assert upretty(expr) == 'ℚ'

    expr = RR

    assert pretty(expr) == 'RR'
    assert upretty(expr) == 'ℝ'

    assert upretty(RealField(prec=100)) == 'ℝ₁₀₀'

    expr = QQ.inject(x)

    assert pretty(expr) == 'QQ[x]'
    assert upretty(expr) == 'ℚ[x]'

    expr = QQ.inject(x, y)

    assert pretty(expr) == 'QQ[x, y]'
    assert upretty(expr) == 'ℚ[x, y]'

    expr = ZZ.inject(x).field

    assert pretty(expr) == 'ZZ(x)'
    assert upretty(expr) == 'ℤ(x)'

    expr = ZZ.inject(x, y).field

    assert pretty(expr) == 'ZZ(x, y)'
    assert upretty(expr) == 'ℤ(x, y)'

    expr = QQ.poly_ring(x, y, order=grlex)

    assert pretty(expr) == 'QQ[x, y, order=grlex]'
    assert upretty(expr) == 'ℚ[x, y, order=grlex]'

    expr = QQ.poly_ring(x, y, order=ilex)

    assert pretty(expr) == 'QQ[x, y, order=ilex]'
    assert upretty(expr) == 'ℚ[x, y, order=ilex]'

    expr = ZZ.frac_field(x, y, order=grlex)
    assert upretty(expr) == 'ℤ(x, y, order=grlex)'


def test_pretty_prec():
    assert xpretty(Float(0.3), full_prec=True, wrap_line=False) == '0.300000000000000'
    assert xpretty(Float(0.3), full_prec='auto', wrap_line=False) == '0.300000000000000'
    assert xpretty(Float(0.3), full_prec=False, wrap_line=False) == '0.3'
    assert xpretty(Float(0.3)*x, full_prec=True, use_unicode=False, wrap_line=False) in [
        '0.300000000000000*x',
        'x*0.300000000000000'
    ]
    assert xpretty(Float('0.3')*x, full_prec='auto', use_unicode=False, wrap_line=False) in [
        '0.3*x',
        'x*0.3'
    ]
    assert xpretty(Float('0.3')*x, full_prec=False, use_unicode=False, wrap_line=False) in [
        '0.3*x',
        'x*0.3'
    ]


def test_pprint():
    fd = StringIO()
    sso = sys.stdout
    sys.stdout = fd
    try:
        pprint(pi, use_unicode=False, wrap_line=False)
    finally:
        sys.stdout = sso
    assert fd.getvalue() == 'pi\n'


def test_pretty_class():
    """Test that the printer dispatcher correctly handles classes."""
    class C:
        pass   # C has no .__class__ and this was causing problems

    class D:
        pass

    assert pretty( C ) == str( C )
    assert pretty( D ) == str( D )


def test_pretty_no_wrap_line():
    huge_expr = 0
    for i in range(20):
        huge_expr += i*sin(i + x)
    assert xpretty(huge_expr            ).find('\n') != -1
    assert xpretty(huge_expr, wrap_line=False).find('\n') == -1


def test_settings():
    pytest.raises(TypeError, lambda: pretty(Integer(4), method='garbage'))


def test_pretty_sum():
    n = Symbol('n')
    expr = Sum(k**k, (k, 0, n))
    ascii_str = \
        """\
  n     \n\
 ___    \n\
 \\  `   \n\
  \\    k\n\
  /   k \n\
 /__,   \n\
k = 0   \
"""
    ucode_str = \
        """\
  n     \n\
 ___    \n\
 ╲      \n\
  ╲    k\n\
  ╱   k \n\
 ╱      \n\
 ‾‾‾    \n\
k = 0   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(k**(Integral(x**n, (x, -oo, oo))), (k, 0, n**n))
    ascii_str = \
        """\
    n             \n\
   n              \n\
______            \n\
\\     `           \n\
 \\        oo      \n\
  \\        /      \n\
   \\      |       \n\
    \\     |   n   \n\
     )    |  x  dx\n\
    /     |       \n\
   /     /        \n\
  /      -oo      \n\
 /      k         \n\
/_____,           \n\
 k = 0            \
"""
    ucode_str = \
        """\
   n            \n\
  n             \n\
______          \n\
╲               \n\
 ╲      ∞       \n\
  ╲     ⌠       \n\
   ╲    ⎮   n   \n\
    ╲   ⎮  x  dx\n\
    ╱   ⌡       \n\
   ╱    -∞      \n\
  ╱    k        \n\
 ╱              \n\
╱               \n\
‾‾‾‾‾‾          \n\
k = 0           \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(k**(
        Integral(x**n, (x, -oo, oo))), (k, 0, Integral(x**x, (x, -oo, oo))))
    ascii_str = \
        """\
 oo                 \n\
  /                 \n\
 |                  \n\
 |   x              \n\
 |  x  dx           \n\
 |                  \n\
/                   \n\
-oo                 \n\
 ______             \n\
 \\     `            \n\
  \\         oo      \n\
   \\         /      \n\
    \\       |       \n\
     \\      |   n   \n\
      )     |  x  dx\n\
     /      |       \n\
    /      /        \n\
   /       -oo      \n\
  /       k         \n\
 /_____,            \n\
  k = 0             \
"""
    ucode_str = \
        """\
∞                 \n\
⌠                 \n\
⎮   x             \n\
⎮  x  dx          \n\
⌡                 \n\
-∞                \n\
 ______           \n\
 ╲                \n\
  ╲       ∞       \n\
   ╲      ⌠       \n\
    ╲     ⎮   n   \n\
     ╲    ⎮  x  dx\n\
     ╱    ⌡       \n\
    ╱     -∞      \n\
   ╱     k        \n\
  ╱               \n\
 ╱                \n\
 ‾‾‾‾‾‾           \n\
 k = 0            \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(k**(Integral(x**n, (x, -oo, oo))), (
        k, x + n + x**2 + n**2 + (x/n) + (1/x), Integral(x**x, (x, -oo, oo))))
    ascii_str = \
        """\
          oo                          \n\
           /                          \n\
          |                           \n\
          |   x                       \n\
          |  x  dx                    \n\
          |                           \n\
         /                            \n\
         -oo                          \n\
          ______                      \n\
          \\     `                     \n\
           \\                  oo      \n\
            \\                  /      \n\
             \\                |       \n\
              \\               |   n   \n\
               )              |  x  dx\n\
              /               |       \n\
             /               /        \n\
            /                -oo      \n\
           /                k         \n\
          /_____,                     \n\
     2        2       1   x           \n\
k = n  + n + x  + x + - + -           \n\
                      x   n           \
"""
    ucode_str = \
        """\
          ∞                          \n\
          ⌠                          \n\
          ⎮   x                      \n\
          ⎮  x  dx                   \n\
          ⌡                          \n\
          -∞                         \n\
           ______                    \n\
           ╲                         \n\
            ╲                ∞       \n\
             ╲               ⌠       \n\
              ╲              ⎮   n   \n\
               ╲             ⎮  x  dx\n\
               ╱             ⌡       \n\
              ╱              -∞      \n\
             ╱              k        \n\
            ╱                        \n\
           ╱                         \n\
           ‾‾‾‾‾‾                    \n\
     2        2       1   x          \n\
k = n  + n + x  + x + ─ + ─          \n\
                      x   n          \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(k**(
        Integral(x**n, (x, -oo, oo))), (k, 0, x + n + x**2 + n**2 + (x/n) + (1/x)))
    ascii_str = \
        """\
 2        2       1   x           \n\
n  + n + x  + x + - + -           \n\
                  x   n           \n\
        ______                    \n\
        \\     `                   \n\
         \\                oo      \n\
          \\                /      \n\
           \\              |       \n\
            \\             |   n   \n\
             )            |  x  dx\n\
            /             |       \n\
           /             /        \n\
          /              -oo      \n\
         /              k         \n\
        /_____,                   \n\
         k = 0                    \
"""
    ucode_str = \
        """\
 2        2       1   x          \n\
n  + n + x  + x + ─ + ─          \n\
                  x   n          \n\
         ______                  \n\
         ╲                       \n\
          ╲              ∞       \n\
           ╲             ⌠       \n\
            ╲            ⎮   n   \n\
             ╲           ⎮  x  dx\n\
             ╱           ⌡       \n\
            ╱            -∞      \n\
           ╱            k        \n\
          ╱                      \n\
         ╱                       \n\
         ‾‾‾‾‾‾                  \n\
         k = 0                   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(x, (x, 0, oo))
    ascii_str = \
        """\
  oo   \n\
 __    \n\
 \\ `   \n\
  )   x\n\
 /_,   \n\
x = 0  \
"""
    ucode_str = \
        """\
  ∞    \n\
 ___   \n\
 ╲     \n\
  ╲   x\n\
  ╱    \n\
 ╱     \n\
 ‾‾‾   \n\
x = 0  \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(x**2, (x, 0, oo))
    ascii_str = \
        """\
  oo    \n\
 ___    \n\
 \\  `   \n\
  \\    2\n\
  /   x \n\
 /__,   \n\
x = 0   \
"""
    ucode_str = \
        """\
  ∞     \n\
 ___    \n\
 ╲      \n\
  ╲    2\n\
  ╱   x \n\
 ╱      \n\
 ‾‾‾    \n\
x = 0   \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(x/2, (x, 0, oo))
    ascii_str = \
        """\
  oo   \n\
 ___   \n\
 \\  `  \n\
  \\   x\n\
   )  -\n\
  /   2\n\
 /__,  \n\
x = 0  \
"""
    ucode_str = \
        """\
  ∞    \n\
 ____  \n\
 ╲     \n\
  ╲   x\n\
   ╲  ─\n\
   ╱  2\n\
  ╱    \n\
 ╱     \n\
 ‾‾‾‾  \n\
x = 0  \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(x**3/2, (x, 0, oo))
    ascii_str = \
        """\
  oo    \n\
____    \n\
\\   `   \n\
 \\     3\n\
  \\   x \n\
  /   --\n\
 /    2 \n\
/___,   \n\
x = 0   \
"""
    ucode_str = \
        """\
  ∞     \n\
 ____   \n\
 ╲      \n\
  ╲    3\n\
   ╲  x \n\
   ╱  ──\n\
  ╱   2 \n\
 ╱      \n\
 ‾‾‾‾   \n\
x = 0   \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum((x**3*y**(x/2))**n, (x, 0, oo))
    ascii_str = \
        """\
  oo          \n\
____          \n\
\\   `         \n\
 \\           n\n\
  \\   /    x\\ \n\
   )  |    -| \n\
  /   | 3  2| \n\
 /    \\x *y / \n\
/___,         \n\
x = 0         \
"""
    ucode_str = \
        """\
  ∞           \n\
_____         \n\
╲             \n\
 ╲           n\n\
  ╲   ⎛    x⎞ \n\
   ╲  ⎜    ─⎟ \n\
   ╱  ⎜ 3  2⎟ \n\
  ╱   ⎝x ⋅y ⎠ \n\
 ╱            \n\
╱             \n\
‾‾‾‾‾         \n\
x = 0         \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(1/x**2, (x, 0, oo))
    ascii_str = \
        """\
  oo    \n\
____    \n\
\\   `   \n\
 \\    1 \n\
  \\   --\n\
  /    2\n\
 /    x \n\
/___,   \n\
x = 0   \
"""
    ucode_str = \
        """\
  ∞     \n\
 ____   \n\
 ╲      \n\
  ╲   1 \n\
   ╲  ──\n\
   ╱   2\n\
  ╱   x \n\
 ╱      \n\
 ‾‾‾‾   \n\
x = 0   \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(1/y**(a/b), (x, 0, oo))
    ascii_str = \
        """\
  oo      \n\
____      \n\
\\   `     \n\
 \\     -a \n\
  \\    ---\n\
  /     b \n\
 /    y   \n\
/___,     \n\
x = 0     \
"""
    ucode_str = \
        """\
  ∞       \n\
 ____     \n\
 ╲        \n\
  ╲    -a \n\
   ╲   ───\n\
   ╱    b \n\
  ╱   y   \n\
 ╱        \n\
 ‾‾‾‾     \n\
x = 0     \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Sum(1/y**(a/b), (x, 0, oo), (y, 1, 2))
    ascii_str = \
        """\
  2     oo     \n\
____  ____     \n\
\\   ` \\   `    \n\
 \\     \\     -a\n\
  \\     \\    --\n\
  /     /    b \n\
 /     /    y  \n\
/___, /___,    \n\
y = 1 x = 0    \
"""
    ucode_str = \
        """\
  2     ∞      \n\
____  ____     \n\
╲     ╲        \n\
 ╲     ╲     -a\n\
  ╲     ╲    ──\n\
  ╱     ╱    b \n\
 ╱     ╱    y  \n\
╱     ╱        \n\
‾‾‾‾  ‾‾‾‾     \n\
y = 1 x = 0    \
"""
    expr = Sum(1/(1 + 1/(
        1 + 1/k)) + 1, (k, 111, 1 + 1/n), (k, 1/(1 + m), oo)) + 1/(1 + 1/k)
    ascii_str = \
        """\
               1                         \n\
           1 + -                         \n\
    oo         n                         \n\
  _____    _____                         \n\
  \\    `   \\    `                        \n\
   \\        \\     /        1    \\        \n\
    \\        \\    |1 + ---------|        \n\
     \\        \\   |          1  |     1  \n\
      )        )  |    1 + -----| + -----\n\
     /        /   |            1|       1\n\
    /        /    |        1 + -|   1 + -\n\
   /        /     \\            k/       k\n\
  /____,   /____,                        \n\
      1   k = 111                        \n\
k = -----                                \n\
    m + 1                                \
"""
    ucode_str = \
        """\
               1                         \n\
           1 + ─                         \n\
    ∞          n                         \n\
  ______   ______                        \n\
  ╲        ╲                             \n\
   ╲        ╲     ⎛        1    ⎞        \n\
    ╲        ╲    ⎜1 + ─────────⎟        \n\
     ╲        ╲   ⎜          1  ⎟        \n\
      ╲        ╲  ⎜    1 + ─────⎟     1  \n\
      ╱        ╱  ⎜            1⎟ + ─────\n\
     ╱        ╱   ⎜        1 + ─⎟       1\n\
    ╱        ╱    ⎝            k⎠   1 + ─\n\
   ╱        ╱                           k\n\
  ╱        ╱                             \n\
  ‾‾‾‾‾‾   ‾‾‾‾‾‾                        \n\
      1   k = 111                        \n\
k = ─────                                \n\
    m + 1                                \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    n = Symbol('n', integer=True)
    expr = Sum(-1/(n**3 - 1), (n, -oo, -3))
    ascii_str = """\
   -3         \n\
 ____         \n\
 \\   `        \n\
  \\      -1   \n\
   \\    ------\n\
   /     3    \n\
  /     n  - 1\n\
 /___,        \n\
n = -oo       \
"""
    ucode_str = """\
  -3         \n\
 ____        \n\
 ╲           \n\
  ╲     -1   \n\
   ╲   ──────\n\
   ╱    3    \n\
  ╱    n  - 1\n\
 ╱           \n\
 ‾‾‾‾        \n\
n = -∞       \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_Subs():
    f = Function('f')
    expr = Subs(f(x), (x, phi**2))
    ascii_str = \
        """\
(f(x))|     2\n\
      |x=phi \
"""
    unicode_str = \
        """\
(f(x))│   2\n\
      │x=φ \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str

    expr = Subs(f(x).diff(x), (x, 0))
    ascii_str = \
        """\
/d       \\|   \n\
|--(f(x))||   \n\
\\dx      /|x=0\
"""
    unicode_str = \
        """\
⎛d       ⎞│   \n\
⎜──(f(x))⎟│   \n\
⎝dx      ⎠│x=0\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str

    expr = Subs(f(x).diff(x)/y, (x, 0), (y, Rational(1, 2)))
    ascii_str = \
        """\
/d       \\|          \n\
|--(f(x))||          \n\
|dx      ||          \n\
|--------||          \n\
\\   y    /|x=0, y=1/2\
"""
    unicode_str = \
        """\
⎛d       ⎞│          \n\
⎜──(f(x))⎟│          \n\
⎜dx      ⎟│          \n\
⎜────────⎟│          \n\
⎝   y    ⎠│x=0, y=1/2\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == unicode_str


def test_gammas():
    assert upretty(lowergamma(x, y)) == 'γ(x, y)'
    assert pretty(lowergamma(x, y)) == 'lowergamma(x, y)'
    assert upretty(uppergamma(x, y)) == 'Γ(x, y)'
    assert pretty(uppergamma(x, y)) == 'uppergamma(x, y)'
    assert upretty(gamma(x)) == 'Γ(x)'
    assert pretty(gamma(x)) == 'gamma(x)'
    # issue sympy/sympy#11841
    assert upretty(Function('gamma')(x)) == 'γ(x)'


def test_deltas():
    assert upretty(DiracDelta(x)) == 'δ(x)'
    assert pretty(DiracDelta(x)) == 'DiracDelta(x)'

    ucode_str = \
        """\
 ⎛x + 1⎞\n\
δ⎜─────⎟\n\
 ⎝x + 2⎠\
"""
    assert upretty(DiracDelta((x + 1)/(x + 2))) == ucode_str

    ucode_str = \
        """\
 (3)   \n\
δ   (x)\
"""
    assert upretty(DiracDelta(x, 3)) == ucode_str


def test_hyper():
    expr = hyper((), (), z)
    ucode_str = \
        """\
 ┌─  ⎛  │  ⎞\n\
 ├─  ⎜  │ z⎟\n\
0╵ 0 ⎝  │  ⎠\
"""
    ascii_str = \
        """\
  _         \n\
 |_  /  |  \\\n\
 |   |  | z|\n\
0  0 \\  |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper((), (1,), x)
    ucode_str = \
        """\
 ┌─  ⎛  │  ⎞\n\
 ├─  ⎜  │ x⎟\n\
0╵ 1 ⎝1 │  ⎠\
"""
    ascii_str = \
        """\
  _         \n\
 |_  /  |  \\\n\
 |   |  | x|\n\
0  1 \\1 |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper([2], [1], x)
    ucode_str = \
        """\
 ┌─  ⎛2 │  ⎞\n\
 ├─  ⎜  │ x⎟\n\
1╵ 1 ⎝1 │  ⎠\
"""
    ascii_str = \
        """\
  _         \n\
 |_  /2 |  \\\n\
 |   |  | x|\n\
1  1 \\1 |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper((pi/3, -2*k), (3, 4, 5, -3), x)
    ucode_str = \
        """\
     ⎛  π         │  ⎞\n\
 ┌─  ⎜  ─, -2⋅k   │  ⎟\n\
 ├─  ⎜  3         │ x⎟\n\
2╵ 4 ⎜            │  ⎟\n\
     ⎝3, 4, 5, -3 │  ⎠\
"""
    ascii_str = \
        """\
                      \n\
  _  /  pi        |  \\\n\
 |_  |  --, -2*k  |  |\n\
 |   |  3         | x|\n\
2  4 |            |  |\n\
     \\3, 4, 5, -3 |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper((pi, Rational(2, 3), -2*k), (3, 4, 5, -3), x**2)
    ucode_str = \
        """\
 ┌─  ⎛π, 2/3, -2⋅k │  2⎞\n\
 ├─  ⎜             │ x ⎟\n\
3╵ 4 ⎝3, 4, 5, -3  │   ⎠\
"""
    ascii_str = \
        """\
  _                      \n\
 |_  /pi, 2/3, -2*k |  2\\\n\
 |   |              | x |\n\
3  4 \\ 3, 4, 5, -3  |   /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = hyper([1, 2], [3, 4], 1/(1/(1/(1/x + 1) + 1) + 1))
    ucode_str = \
        """\
     ⎛     │       1      ⎞\n\
     ⎜     │ ─────────────⎟\n\
     ⎜     │         1    ⎟\n\
 ┌─  ⎜1, 2 │ 1 + ─────────⎟\n\
 ├─  ⎜     │           1  ⎟\n\
2╵ 2 ⎜3, 4 │     1 + ─────⎟\n\
     ⎜     │             1⎟\n\
     ⎜     │         1 + ─⎟\n\
     ⎝     │             x⎠\
"""

    ascii_str = \
        """\
                           \n\
     /     |       1      \\\n\
     |     | -------------|\n\
  _  |     |         1    |\n\
 |_  |1, 2 | 1 + ---------|\n\
 |   |     |           1  |\n\
2  2 |3, 4 |     1 + -----|\n\
     |     |             1|\n\
     |     |         1 + -|\n\
     \\     |             x/\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_meijerg():
    expr = meijerg([pi, pi, x], [1], [0, 1], [1, 2, 3], z)
    ucode_str = \
        """\
╭─╮2, 3 ⎛π, π, x     1    │  ⎞\n\
│╶┐     ⎜                 │ z⎟\n\
╰─╯4, 5 ⎝ 0, 1    1, 2, 3 │  ⎠\
"""
    ascii_str = \
        """\
 __2, 3 /pi, pi, x     1    |  \\\n\
/__     |                   | z|\n\
\\_|4, 5 \\  0, 1     1, 2, 3 |  /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = meijerg([1, pi/7], [2, pi, 5], [], [], z**2)
    ucode_str = \
        """\
        ⎛   π          │   ⎞\n\
╭─╮0, 2 ⎜1, ─  2, π, 5 │  2⎟\n\
│╶┐     ⎜   7          │ z ⎟\n\
╰─╯5, 0 ⎜              │   ⎟\n\
        ⎝              │   ⎠\
"""
    ascii_str = \
        """\
        /   pi           |   \\\n\
 __0, 2 |1, --  2, pi, 5 |  2|\n\
/__     |   7            | z |\n\
\\_|5, 0 |                |   |\n\
        \\                |   /\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ucode_str = \
        """\
╭─╮ 1, 10 ⎛1, 1, 1, 1, 1, 1, 1, 1, 1, 1  1 │  ⎞\n\
│╶┐       ⎜                                │ z⎟\n\
╰─╯11,  2 ⎝             1                1 │  ⎠\
"""
    ascii_str = \
        """\
 __ 1, 10 /1, 1, 1, 1, 1, 1, 1, 1, 1, 1  1 |  \\\n\
/__       |                                | z|\n\
\\_|11,  2 \\             1                1 |  /\
"""

    expr = meijerg([1]*10, [1], [1], [1], z)
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = meijerg([1, 2, ], [4, 3], [3], [4, 5], 1/(1/(1/(1/x + 1) + 1) + 1))

    ucode_str = \
        """\
        ⎛           │       1      ⎞\n\
        ⎜           │ ─────────────⎟\n\
        ⎜           │         1    ⎟\n\
╭─╮1, 2 ⎜1, 2  4, 3 │ 1 + ─────────⎟\n\
│╶┐     ⎜           │           1  ⎟\n\
╰─╯4, 3 ⎜ 3    4, 5 │     1 + ─────⎟\n\
        ⎜           │             1⎟\n\
        ⎜           │         1 + ─⎟\n\
        ⎝           │             x⎠\
"""

    ascii_str = \
        """\
        /           |       1      \\\n\
        |           | -------------|\n\
        |           |         1    |\n\
 __1, 2 |1, 2  4, 3 | 1 + ---------|\n\
/__     |           |           1  |\n\
\\_|4, 3 | 3    4, 5 |     1 + -----|\n\
        |           |             1|\n\
        |           |         1 + -|\n\
        \\           |             x/\
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = Integral(expr, x)

    ucode_str = \
        """\
⌠                                        \n\
⎮         ⎛           │       1      ⎞   \n\
⎮         ⎜           │ ─────────────⎟   \n\
⎮         ⎜           │         1    ⎟   \n\
⎮ ╭─╮1, 2 ⎜1, 2  4, 3 │ 1 + ─────────⎟   \n\
⎮ │╶┐     ⎜           │           1  ⎟ dx\n\
⎮ ╰─╯4, 3 ⎜ 3    4, 5 │     1 + ─────⎟   \n\
⎮         ⎜           │             1⎟   \n\
⎮         ⎜           │         1 + ─⎟   \n\
⎮         ⎝           │             x⎠   \n\
⌡                                        \
"""

    ascii_str = \
        """\
  /                                       \n\
 |                                        \n\
 |         /           |       1      \\   \n\
 |         |           | -------------|   \n\
 |         |           |         1    |   \n\
 |  __1, 2 |1, 2  4, 3 | 1 + ---------|   \n\
 | /__     |           |           1  | dx\n\
 | \\_|4, 3 | 3    4, 5 |     1 + -----|   \n\
 |         |           |             1|   \n\
 |         |           |         1 + -|   \n\
 |         \\           |             x/   \n\
 |                                        \n\
/                                         \
"""

    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_noncommutative():
    A, B, C = symbols('A,B,C', commutative=False)

    expr = A*B*C**-1
    ascii_str = \
        """\
     -1\n\
A*B*C  \
"""
    ucode_str = \
        """\
     -1\n\
A⋅B⋅C  \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = C**-1*A*B
    ascii_str = \
        """\
 -1    \n\
C  *A*B\
"""
    ucode_str = \
        """\
 -1    \n\
C  ⋅A⋅B\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = A*C**-1*B
    ascii_str = \
        """\
   -1  \n\
A*C  *B\
"""
    ucode_str = \
        """\
   -1  \n\
A⋅C  ⋅B\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    expr = A*C**-1*B/x
    ascii_str = \
        """\
   -1  \n\
A*C  *B\n\
-------\n\
   x   \
"""
    ucode_str = \
        """\
   -1  \n\
A⋅C  ⋅B\n\
───────\n\
   x   \
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_special_functions():
    # atan2
    expr = atan2(y/sqrt(200), sqrt(x))
    ascii_str = \
        """\
     /  ___         \\\n\
     |\\/ 2 *y    ___|\n\
atan2|-------, \\/ x |\n\
     \\   20         /\
"""
    ucode_str = \
        """\
     ⎛  ___         ⎞\n\
     ⎜╲╱ 2 ⋅y    ___⎟\n\
atan2⎜───────, ╲╱ x ⎟\n\
     ⎝   20         ⎠\
"""
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_pretty_geometry():
    e = Segment((0, 1), (0, 2))
    assert pretty(e) == 'Segment(Point(0, 1), Point(0, 2))'
    e = Ray((1, 1), angle=4.02*pi)
    assert pretty(e) == 'Ray(Point(1, 1), Point(2, tan(pi/50) + 1))'


def test_expint():
    expr = Ei(x)
    string = 'Ei(x)'
    assert pretty(expr) == string
    assert upretty(expr) == string

    expr = expint(1, z)
    ucode_str = 'E₁(z)'
    ascii_str = 'expint(1, z)'
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    assert pretty(Shi(x)) == 'Shi(x)'
    assert pretty(Si(x)) == 'Si(x)'
    assert pretty(Ci(x)) == 'Ci(x)'
    assert pretty(Chi(x)) == 'Chi(x)'
    assert upretty(Shi(x)) == 'Shi(x)'
    assert upretty(Si(x)) == 'Si(x)'
    assert upretty(Ci(x)) == 'Ci(x)'
    assert upretty(Chi(x)) == 'Chi(x)'


def test_elliptic_functions():
    ascii_str = \
        """\
 /  1  \\\n\
K|-----|\n\
 \\z + 1/\
"""
    ucode_str = \
        """\
 ⎛  1  ⎞\n\
K⎜─────⎟\n\
 ⎝z + 1⎠\
"""
    expr = elliptic_k(1/(z + 1))
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
        """\
 / |  1  \\\n\
F|1|-----|\n\
 \\ |z + 1/\
"""
    ucode_str = \
        """\
 ⎛ │  1  ⎞\n\
F⎜1│─────⎟\n\
 ⎝ │z + 1⎠\
"""
    expr = elliptic_f(1, 1/(1 + z))
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
        """\
 /  1  \\\n\
E|-----|\n\
 \\z + 1/\
"""
    ucode_str = \
        """\
 ⎛  1  ⎞\n\
E⎜─────⎟\n\
 ⎝z + 1⎠\
"""
    expr = elliptic_e(1/(z + 1))
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
        """\
 / |  1  \\\n\
E|1|-----|\n\
 \\ |z + 1/\
"""
    ucode_str = \
        """\
 ⎛ │  1  ⎞\n\
E⎜1│─────⎟\n\
 ⎝ │z + 1⎠\
"""
    expr = elliptic_e(1, 1/(1 + z))
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
        """\
  / |4\\\n\
Pi|3|-|\n\
  \\ |x/\
"""
    ucode_str = \
        """\
 ⎛ │4⎞\n\
Π⎜3│─⎟\n\
 ⎝ │x⎠\
"""
    expr = elliptic_pi(3, 4/x)
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str

    ascii_str = \
        """\
  /   4| \\\n\
Pi|3; -|6|\n\
  \\   x| /\
"""
    ucode_str = \
        """\
 ⎛   4│ ⎞\n\
Π⎜3; ─│6⎟\n\
 ⎝   x│ ⎠\
"""
    expr = elliptic_pi(3, 4/x, 6)
    assert pretty(expr) == ascii_str
    assert upretty(expr) == ucode_str


def test_RandomDomain():
    X = Normal('x1', 0, 1)
    assert upretty(where(X > 0)) == 'Domain: 0 < x₁ ∧ x₁ < ∞'

    D = Die('d1', 6)
    assert upretty(where(D > 4)) == 'Domain: d₁ = 5 ∨ d₁ = 6'

    A = Exponential('a', 1)
    B = Exponential('b', 1)
    assert upretty(pspace(Tuple(A, B)).domain) == \
        'Domain: 0 ≤ a ∧ 0 ≤ b ∧ a < ∞ ∧ b < ∞'


def test_PrettyPoly():
    R = QQ.inject(x, y)
    F = R.field

    expr = F.convert(x/(x + y))
    assert pretty(expr) == 'x/(x + y)'
    assert upretty(expr) == 'x/(x + y)'

    expr = R.convert(x + y)
    assert pretty(expr) == 'x + y'
    assert upretty(expr) == 'x + y'


def test_sympyissue_6285():
    assert pretty(Pow(2, -5, evaluate=False)) == '1 \n--\n 5\n2 '
    assert pretty(Pow(x, (1/pi))) == 'pi___\n\\/ x '


def test_sympyissue_6359():
    assert pretty(Integral(x**2, x)**2) == \
        """\
          2
/  /     \\ \n\
| |      | \n\
| |  2   | \n\
| | x  dx| \n\
| |      | \n\
\\/       / \
"""
    assert upretty(Integral(x**2, x)**2) == \
        """\
         2
⎛⌠      ⎞ \n\
⎜⎮  2   ⎟ \n\
⎜⎮ x  dx⎟ \n\
⎝⌡      ⎠ \
"""

    assert pretty(Sum(x**2, (x, 0, 1))**2) == \
        """\
          2
/  1     \\ \n\
| ___    | \n\
| \\  `   | \n\
|  \\    2| \n\
|  /   x | \n\
| /__,   | \n\
\\x = 0   / \
"""
    assert upretty(Sum(x**2, (x, 0, 1))**2) == \
        """\
          2
⎛  1     ⎞ \n\
⎜ ___    ⎟ \n\
⎜ ╲      ⎟ \n\
⎜  ╲    2⎟ \n\
⎜  ╱   x ⎟ \n\
⎜ ╱      ⎟ \n\
⎜ ‾‾‾    ⎟ \n\
⎝x = 0   ⎠ \
"""

    assert pretty(Product(x**2, (x, 1, 2))**2) == \
        """\
           2
/  2      \\ \n\
|______   | \n\
||    |  2| \n\
||    | x | \n\
||    |   | \n\
\\x = 1    / \
"""
    assert upretty(Product(x**2, (x, 1, 2))**2) == \
        """\
           2
⎛  2      ⎞ \n\
⎜┬────┬   ⎟ \n\
⎜│    │  2⎟ \n\
⎜│    │ x ⎟ \n\
⎜│    │   ⎟ \n\
⎝x = 1    ⎠ \
"""

    f = Function('f')
    assert pretty(Derivative(f(x), x)**2) == \
        """\
          2
/d       \\ \n\
|--(f(x))| \n\
\\dx      / \
"""
    assert upretty(Derivative(f(x), x)**2) == \
        """\
          2
⎛d       ⎞ \n\
⎜──(f(x))⎟ \n\
⎝dx      ⎠ \
"""


def test_sympyissue_6739():
    ascii_str = \
        """\
  1  \n\
-----\n\
  ___\n\
\\/ x \
"""
    ucode_str = \
        """\
  1  \n\
─────\n\
  ___\n\
╲╱ x \
"""
    assert pretty(1/sqrt(x)) == ascii_str
    assert upretty(1/sqrt(x)) == ucode_str


def test_complicated_symbol_unchanged():
    for symb_name in ['dexpr2_d1tau', 'dexpr2^d1tau']:
        assert pretty(Symbol(symb_name)) == symb_name


def test_Tr():
    A, B = symbols('A B', commutative=False)
    t = Tr(A*B)
    assert pretty(t) == r'Tr(A*B)'
    assert upretty(t) == 'Tr(A⋅B)'


def test_pretty_Add():
    eq = Mul(-2, x - 2, evaluate=False) + 5
    assert pretty(eq) == '-2*(x - 2) + 5'
    eq = Add(0, 0, evaluate=False)
    assert pretty(eq) == '0 + 0'
    assert upretty(eq) == '0 + 0'

    eq = Add(y, x, evaluate=False)
    assert upretty(eq, order='none') == 'y + x'


def test_pretty_Mul():
    eq = Mul(1, 1, evaluate=False)
    assert pretty(eq) == '1*1'
    assert upretty(eq) == '1⋅1'

    eq = Mul(y, x, evaluate=False)
    assert upretty(eq, order='none') == 'y⋅x'


def test_sympyissue_7180():
    assert upretty(Equivalent(x, y)) == 'x ≡ y'


def test_pretty_Complement():
    assert pretty(S.Reals - S.Naturals) == r'(-oo, oo) \ Naturals()'
    assert upretty(S.Reals - S.Naturals) == r'ℝ \ ℕ'
    assert pretty(S.Reals - S.Naturals0) == r'(-oo, oo) \ Naturals0()'
    assert upretty(S.Reals - S.Naturals0) == r'ℝ \ ℕ₀'


def test_pretty_SymmetricDifference():
    assert upretty(SymmetricDifference(Interval(2, 3), Interval(3, 5),
                                       evaluate=False)) == '[2, 3] ∆ [3, 5]'
    with pytest.raises(NotImplementedError):
        pretty(SymmetricDifference(Interval(2, 3), Interval(3, 5), evaluate=False))


def test_pretty_Contains():
    assert pretty(Contains(x, S.Integers)) == 'Contains(x, Integers())'
    assert upretty(Contains(x, S.Integers)) == 'x ∈ ℤ'


def test_sympyissue_8292():
    e = Add(Mul(Add(x**4, x),
                Pow(Add(x, -1, evaluate=False), -1,
                    evaluate=False),
                evaluate=False),
            Mul(-2, Pow(Add(x, -1, evaluate=False), -4, evaluate=False),
                Pow(Add(x, -1, evaluate=False), 4, evaluate=False),
                evaluate=False),
            evaluate=False)
    ucode_str = \
        """\
           4    4    \n\
  2⋅(x - 1)    x  + x\n\
- ────────── + ──────\n\
          4    x - 1 \n\
   (x - 1)           \
"""
    ascii_str = \
        """\
           4    4    \n\
  2*(x - 1)    x  + x\n\
- ---------- + ------\n\
          4    x - 1 \n\
   (x - 1)           \
"""
    assert pretty(e) == ascii_str
    assert upretty(e) == ucode_str


def test_sympyissue_8344():
    e = Add(Mul(2, x, y**2,
                Pow(Pow(1, 2, evaluate=False),
                    -1, evaluate=False), evaluate=False),
            1, evaluate=False)
    ucode_str = \
        """\
     2    \n\
2⋅x⋅y     \n\
────── + 1\n\
   2      \n\
  1       \
"""
    assert upretty(e) == ucode_str


def test_sympyissue_6324():
    x = Pow(2, 3, evaluate=False)
    y = Pow(10, -2, evaluate=False)
    e = Mul(x, y, evaluate=False)
    ucode_str = \
        """\
  3\n\
 2 \n\
───\n\
  2\n\
10 \
"""
    assert upretty(e) == ucode_str


def test_sympyissue_6134():
    phi = Function('phi')

    e = lamda*x*Integral(phi(t)*pi*sin(pi*t), (t, 0, 1)) + lamda*x**2*Integral(phi(t)*2*pi*sin(2*pi*t), (t, 0, 1))
    ucode_str = \
        """\
     1                              1                   \n\
   2 ⌠                              ⌠                   \n\
λ⋅x ⋅⎮ 2⋅π⋅φ(t)⋅sin(2⋅π⋅t) dt + λ⋅x⋅⎮ π⋅φ(t)⋅sin(π⋅t) dt\n\
     ⌡                              ⌡                   \n\
     0                              0                   \
"""
    assert upretty(e) == ucode_str


def test_BaseVectorField():
    assert upretty(BaseVectorField(R2_r, 1)) == '∂_y'


def test_BaseScalarField():
    v = BaseVectorField(R2_r, 1)
    g = Function('g')
    s_field = g(R2.x, R2.y)
    assert pretty(v(s_field)) == \
        """\
/  d              \\|      \n\
|-----(g(x, xi_2))||      \n\
\\dxi_2            /|xi_2=y\
"""


def test_MatrixElement():
    X = MatrixSymbol('X', 2, 2)
    n = Symbol('n', integer=True)
    assert pretty(X[n, 0]) == 'X[n, 0]'

    ucode_str = \
        """\
⎡X₀₀  1⎤\n\
⎢      ⎥\n\
⎣ 2   3⎦\
"""
    assert upretty(Matrix([[X[0, 0], 1], [2, 3]])) == ucode_str


def test_AlgebraicElement():
    K = QQ.algebraic_field(sqrt(2))
    a = K([1, 1])
    ucode_str = \
        """\
      ___\n\
1 + ╲╱ 2 \
"""
    assert upretty(a) == ucode_str


def test_sympyissue_11801():
    assert pretty(Symbol('')) == ''
    assert upretty(Symbol('')) == ''


def test_Order():
    assert upretty(O(1)) == 'O(1)'


def test_pretty_Mod():
    ascii_str1 = 'x mod 7'
    ucode_str1 = 'x mod 7'

    ascii_str2 = '(x + 1) mod 7'
    ucode_str2 = '(x + 1) mod 7'

    ascii_str3 = '2*x mod 7'
    ucode_str3 = '2⋅x mod 7'

    ascii_str4 = '(x mod 7) + 1'
    ucode_str4 = '(x mod 7) + 1'

    ascii_str5 = '2*(x mod 7)'
    ucode_str5 = '2⋅(x mod 7)'

    x = symbols('x', integer=True)
    assert pretty(Mod(x, 7)) == ascii_str1
    assert upretty(Mod(x, 7)) == ucode_str1
    assert pretty(Mod(x + 1, 7)) == ascii_str2
    assert upretty(Mod(x + 1, 7)) == ucode_str2
    assert pretty(Mod(2 * x, 7)) == ascii_str3
    assert upretty(Mod(2 * x, 7)) == ucode_str3
    assert pretty(Mod(x, 7) + 1) == ascii_str4
    assert upretty(Mod(x, 7) + 1) == ucode_str4
    assert pretty(2 * Mod(x, 7)) == ascii_str5
    assert upretty(2 * Mod(x, 7)) == ucode_str5


def test_Pow_roots():
    expr = root(pi, E)

    assert pretty(expr) == \
        """\
E ____\n\
\\/ pi \
"""
    assert upretty(expr) == \
        """\
ℯ ___\n\
╲╱ π \
"""

    expr = root(pi, pi)

    assert pretty(expr) == \
        """\
pi____\n\
\\/ pi \
"""
    assert upretty(expr) == \
        """\
π ___\n\
╲╱ π \
"""

    expr = root(7, pi)  # issue diofant/diofant#888

    assert pretty(expr) == \
        """\
pi___\n\
\\/ 7 \
"""
    assert upretty(expr) == \
        """\
π ___\n\
╲╱ 7 \
"""

    expr = root(pi, EulerGamma)

    assert upretty(expr) == \
        """\
γ ___\n\
╲╱ π \
"""
    assert pretty(expr) == \
        """\
EulerGamma____\n\
        \\/ pi \
"""

    expr = root(7, Symbol('x_17'))

    assert upretty(expr) == \
        """\
x₁₇___\n\
 ╲╱ 7 \
"""
    assert pretty(expr) == \
        """\
x_17___\n\
  \\/ 7 \
"""
