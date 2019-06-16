import pytest

from diofant import (CC, FF, QQ, ZZ, Abs, Add, BlockMatrix, Chi, Ci,
                     Complement, Contains, CosineTransform, Dict, Dummy, Ei,
                     Eq, Equivalent, FallingFactorial, FiniteSet, Float,
                     FourierTransform, Function, Identity, ImageSet,
                     IndexedBase, Integer, Integral, Intersection, Interval,
                     InverseCosineTransform, InverseFourierTransform,
                     InverseLaplaceTransform, InverseMellinTransform,
                     InverseSineTransform, Lambda, LaplaceTransform, Limit,
                     Matrix, MatrixSymbol, Max, MellinTransform, Min, Mod, Mul,
                     Not, Order, Piecewise, Poly, Pow, Product, Range,
                     Rational, RisingFactorial, RootOf, RootSum, S, Shi, Si,
                     SineTransform, Subs, Sum, Symbol, SymmetricDifference,
                     Tuple, Union, Wild, Ynm, Znm, acot, airyai, airyaiprime,
                     airybi, airybiprime, arg, asin, assoc_laguerre,
                     assoc_legendre, binomial, catalan, cbrt, ceiling,
                     chebyshevt, chebyshevu, conjugate, cos, cot, coth, diff,
                     dirichlet_eta, divisor_sigma, elliptic_e, elliptic_f,
                     elliptic_k, elliptic_pi, exp, expint, factorial,
                     factorial2, false, field, floor, gamma, gegenbauer,
                     hermite, hyper, im, jacobi, laguerre, legendre, lerchphi,
                     log, lowergamma, meijerg, oo, pi, polar_lift, polylog, re,
                     ring, root, sin, sqrt, subfactorial, symbols, tan,
                     totient, true, uppergamma, zeta)
from diofant.abc import a, b, mu, t, tau, w, x, y, z
from diofant.combinatorics.permutations import Cycle, Permutation
from diofant.core.trace import Tr
from diofant.diffgeom import (CovarDerivativeOp, Differential, Manifold, Patch,
                              TensorProduct, metric_to_Christoffel_2nd)
from diofant.diffgeom.rn import R2, R2_r
from diofant.functions import (DiracDelta, Heaviside, KroneckerDelta,
                               LeviCivita, besseli, besselj, besselk, bessely,
                               euler, fresnelc, fresnels, hankel1, hankel2, jn,
                               yn)
from diofant.logic import Implies
from diofant.logic.boolalg import And, Or, Xor
from diofant.matrices import (Adjoint, HadamardProduct, Inverse, Transpose,
                              ZeroMatrix)
from diofant.parsing.sympy_parser import parse_expr
from diofant.printing.latex import (LatexPrinter, greek_letters_set, latex,
                                    other_symbols, tex_greek_dictionary,
                                    translate)
from diofant.stats import Die, Exponential, Normal, pspace, where
from diofant.tensor import (ImmutableDenseNDimArray, ImmutableSparseNDimArray,
                            MutableDenseNDimArray, MutableSparseNDimArray,
                            tensorproduct)


__all__ = ()

k, m, n = symbols('k m n', integer=True)


def test_printmethod():
    class R(Abs):
        def _latex(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert latex(R(x)) == "foo(x)"

    class R(Abs):
        def _latex(self, printer):
            return "foo"
    assert latex(R(x)) == "foo"


def test_latex_basic():
    assert latex(1 + x) == "x + 1"
    assert latex(x**2) == "x^{2}"
    assert latex(x**(1 + x)) == "x^{x + 1}"
    assert latex(x**3 + x + 1 + x**2) == "x^{3} + x^{2} + x + 1"

    assert latex(2*x*y) == "2 x y"
    assert latex(2*x*y, mul_symbol='dot') == r"2 \cdot x \cdot y"

    assert latex(1/x) == r"\frac{1}{x}"
    assert latex(1/x, fold_short_frac=True) == "1 / x"
    assert latex(-Rational(3, 2)) == r"- \frac{3}{2}"
    assert latex(-Rational(3, 2), fold_short_frac=True) == r"- 3 / 2"
    assert latex(1/x**2) == r"\frac{1}{x^{2}}"
    assert latex(x/2) == r"\frac{x}{2}"
    assert latex(x/2, fold_short_frac=True) == "x / 2"
    assert latex((x + y)/(2*x)) == r"\frac{x + y}{2 x}"
    assert latex((x + y)/(2*x), fold_short_frac=True) == \
        r"\left(x + y\right) / 2 x"
    assert latex((x + y)/(2*x), long_frac_ratio=0) == \
        r"\frac{1}{2 x} \left(x + y\right)"
    assert latex((x + y)/x) == r"\frac{1}{x} \left(x + y\right)"
    assert latex((x + y)/x, long_frac_ratio=3) == r"\frac{x + y}{x}"

    assert latex(2*Integral(x, x)/3) == r"\frac{2}{3} \int x\, dx"
    assert latex(2*Integral(x, x)/3, fold_short_frac=True) == \
        r"\left(2 \int x\, dx\right) / 3"
    assert latex(Integral(x, x)/2/pi, long_frac_ratio=0) == \
        r"\frac{1}{2 \pi} \int x\, dx"

    assert latex(sqrt(x)) == r"\sqrt{x}"
    assert latex(cbrt(x)) == r"\sqrt[3]{x}"
    assert latex(sqrt(x)**3) == r"x^{\frac{3}{2}}"
    assert latex(sqrt(x), itex=True) == r"\sqrt{x}"
    assert latex(cbrt(x), itex=True) == r"\root{3}{x}"
    assert latex(sqrt(x)**3, itex=True) == r"x^{\frac{3}{2}}"
    assert latex(x**Rational(3, 4)) == r"x^{\frac{3}{4}}"
    assert latex(x**Rational(3, 4), fold_frac_powers=True) == "x^{3/4}"
    assert latex((x + 1)**Rational(3, 4)) == \
        r"\left(x + 1\right)^{\frac{3}{4}}"
    assert latex((x + 1)**Rational(3, 4), fold_frac_powers=True) == \
        r"\left(x + 1\right)^{3/4}"

    assert latex(1.5e20*x) == r"1.5 \cdot 10^{20} x"
    assert latex(1.5e20*x, mul_symbol='dot') == r"1.5 \cdot 10^{20} \cdot x"
    assert latex(1.5e20*x, mul_symbol='times') == r"1.5 \times 10^{20} \times x"

    assert latex(x*y, order='none') == 'x y'

    assert latex(1/sin(x)) == r"\frac{1}{\sin{\left (x \right )}}"
    assert latex(sin(x)**-1) == r"\frac{1}{\sin{\left (x \right )}}"
    assert latex(sin(x)**Rational(3, 2)) == \
        r"\sin^{\frac{3}{2}}{\left (x \right )}"
    assert latex(sin(x)**Rational(3, 2), fold_frac_powers=True) == \
        r"\sin^{3/2}{\left (x \right )}"

    assert latex(~x) == r"\neg x"
    assert latex(x & y) == r"x \wedge y"
    assert latex(x & y & z) == r"x \wedge y \wedge z"
    assert latex(x | y) == r"x \vee y"
    assert latex(x | y | z) == r"x \vee y \vee z"
    assert latex((x & y) | z) == r"z \vee \left(x \wedge y\right)"
    assert latex(Implies(x, y)) == r"x \Rightarrow y"
    assert latex(~(x >> ~y)) == r"x \not\Rightarrow \neg y"
    assert latex(Implies(Or(x, y), z)) == r"\left(x \vee y\right) \Rightarrow z"
    assert latex(Implies(z, Or(x, y))) == r"z \Rightarrow \left(x \vee y\right)"

    assert latex(~x, symbol_names={x: "x_i"}) == r"\neg x_i"
    assert latex(x & y, symbol_names={x: "x_i", y: "y_i"}) == \
        r"x_i \wedge y_i"
    assert latex(x & y & z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) == \
        r"x_i \wedge y_i \wedge z_i"
    assert latex(x | y, symbol_names={x: "x_i", y: "y_i"}) == r"x_i \vee y_i"
    assert latex(x | y | z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) == \
        r"x_i \vee y_i \vee z_i"
    assert latex((x & y) | z, symbol_names={x: "x_i", y: "y_i", z: "z_i"}) == \
        r"z_i \vee \left(x_i \wedge y_i\right)"
    assert latex(Implies(x, y), symbol_names={x: "x_i", y: "y_i"}) == \
        r"x_i \Rightarrow y_i"
    assert latex(Tuple(x, y)) == r'\left ( x, \quad y\right )'


def test_latex_builtins():
    assert latex(True) == r"\mathrm{True}"
    assert latex(False) == r"\mathrm{False}"
    assert latex(None) == r"\mathrm{None}"
    assert latex(true) == r"\mathrm{true}"
    assert latex(false) == r'\mathrm{false}'


def test_latex_cycle():
    assert latex(Cycle(1, 2, 4)) == r"\left( 1\; 2\; 4\right)"
    assert latex(Cycle(1, 2)(4, 5, 6)) == r"\left( 1\; 2\right)\left( 4\; 5\; 6\right)"
    assert latex(Cycle()) == r"\left( \right)"


def test_latex_permutation():
    assert latex(Permutation(1, 2, 4)) == r"\left( 1\; 2\; 4\right)"
    assert latex(Permutation(1, 2)(4, 5, 6)) == r"\left( 1\; 2\right)\left( 4\; 5\; 6\right)"
    assert latex(Permutation()) == r"\left( \right)"
    assert latex(Permutation(2, 4)*Permutation(5)) == r"\left( 2\; 4\right)\left( 5\right)"
    assert latex(Permutation(5)) == r"\left( 5\right)"


def test_latex_Float():
    assert latex(Float(1.0e100)) == r"1.0 \cdot 10^{100}"
    assert latex(Float(1.0e-100)) == r"1.0 \cdot 10^{-100}"
    assert latex(Float(1.0e-100), mul_symbol="times") == r"1.0 \times 10^{-100}"
    assert latex(1.0*oo) == r"\infty"
    assert latex(-1.0*oo) == r"- \infty"


def test_latex_symbols():
    Gamma, lmbda, rho = symbols('Gamma, lambda, rho')
    tau, Tau, TAU, taU = symbols('tau, Tau, TAU, taU')
    assert latex(tau) == r"\tau"
    assert latex(Tau) == "T"
    assert latex(TAU) == r"\tau"
    assert latex(taU) == r"\tau"
    # Check that all capitalized greek letters are handled explicitly
    capitalized_letters = {l.capitalize() for l in greek_letters_set}
    assert len(capitalized_letters - set(tex_greek_dictionary.keys())) == 0
    assert latex(Gamma + lmbda) == r"\Gamma + \lambda"
    assert latex(Gamma * lmbda) == r"\Gamma \lambda"
    assert latex(Symbol('q1')) == r"q_{1}"
    assert latex(Symbol('q21')) == r"q_{21}"
    assert latex(Symbol('epsilon0')) == r"\epsilon_{0}"
    assert latex(Symbol('omega1')) == r"\omega_{1}"
    assert latex(Symbol('91')) == r"91"
    assert latex(Symbol('alpha_new')) == r"\alpha_{new}"
    assert latex(Symbol('C^orig')) == r"C^{orig}"
    assert latex(Symbol('x^alpha')) == r"x^{\alpha}"
    assert latex(Symbol('beta^alpha')) == r"\beta^{\alpha}"
    assert latex(Symbol('e^Alpha')) == r"e^{A}"
    assert latex(Symbol('omega_alpha^beta')) == r"\omega^{\beta}_{\alpha}"
    assert latex(Symbol('omega') ** Symbol('beta')) == r"\omega^{\beta}"

    assert latex(Dummy('x')) == 'x'


def test_latex_functions():
    assert latex(exp(x)) == "e^{x}"
    assert latex(exp(1) + exp(2)) == "e + e^{2}"

    f = Function('f')
    assert latex(f(x)) == r'f{\left (x \right )}'
    assert latex(f) == r'f'

    g = Function('g')
    assert latex(g(x, y)) == r'g{\left (x,y \right )}'
    assert latex(g) == r'g'

    h = Function('h')
    assert latex(h(x, y, z)) == r'h{\left (x,y,z \right )}'
    assert latex(h) == r'h'

    Li = Function('Li')
    assert latex(Li) == r'\operatorname{Li}'
    assert latex(Li(x)) == r'\operatorname{Li}{\left (x \right )}'

    beta = Function('beta')

    # not to be confused with the beta function
    assert latex(beta(x)) == r"\beta{\left (x \right )}"
    assert latex(beta) == r"\beta"

    a1 = Function('a_1')

    assert latex(a1) == r"\operatorname{a_{1}}"
    assert latex(a1(x)) == r"\operatorname{a_{1}}{\left (x \right )}"

    # issue sympy/sympy#5868
    omega1 = Function('omega1')
    assert latex(omega1) == r"\omega_{1}"
    assert latex(omega1(x)) == r"\omega_{1}{\left (x \right )}"

    assert latex(sin(x)) == r"\sin{\left (x \right )}"
    assert latex(sin(x), fold_func_brackets=True) == r"\sin {x}"
    assert latex(sin(2*x**2), fold_func_brackets=True) == \
        r"\sin {2 x^{2}}"
    assert latex(sin(x**2), fold_func_brackets=True) == \
        r"\sin {x^{2}}"

    assert latex(asin(x)**2) == r"\operatorname{asin}^{2}{\left (x \right )}"
    assert latex(asin(x)**2, inv_trig_style="full") == \
        r"\arcsin^{2}{\left (x \right )}"
    assert latex(asin(x)**2, inv_trig_style="power") == \
        r"\sin^{-1}{\left (x \right )}^{2}"
    assert latex(asin(x**2), inv_trig_style="power",
                 fold_func_brackets=True) == \
        r"\sin^{-1} {x^{2}}"
    pytest.raises(ValueError, lambda: latex(asin(x)**2,
                                            inv_trig_style="spam"))
    assert (latex(acot(x)**2, fold_func_brackets=True) ==
            r'\operatorname{acot}^{2}x')

    assert latex(factorial(k)) == r"k!"
    assert latex(factorial(k)**2) == r'k!^{2}'
    assert latex(factorial(-k)) == r"\left(- k\right)!"

    assert latex(subfactorial(k)) == r"!k"
    assert latex(subfactorial(k)**2) == r'!k^{2}'
    assert latex(subfactorial(-k)) == r"!\left(- k\right)"

    assert latex(factorial2(k)) == r"k!!"
    assert latex(factorial2(-k)) == r"\left(- k\right)!!"
    assert latex(factorial2(k)**2) == r'k!!^{2}'

    assert latex(binomial(2, k)) == r"{\binom{2}{k}}"
    assert latex(binomial(2, k)**3) == r'{\binom{2}{k}}^{3}'

    assert latex(FallingFactorial(3, k)) == r"{\left(3\right)}_{k}"
    assert (latex(FallingFactorial(3, k)**2) ==
            r'\left({\left(3\right)}_{k}\right)^{2}')
    assert (latex(FallingFactorial(x, 1/y)) ==
            r'{\left(x\right)}_{\left(\frac{1}{y}\right)}')
    assert latex(RisingFactorial(3, k)) == r"{3}^{\left(k\right)}"
    assert (latex(RisingFactorial(3, k)**2) ==
            r'\left({3}^{\left(k\right)}\right)^{2}')
    assert (latex(RisingFactorial(1/x, y)) ==
            r'{\left(\frac{1}{x}\right)}^{\left(y\right)}')

    assert latex(floor(x)) == r"\lfloor{x}\rfloor"
    assert latex(floor(x)**2) == r"\lfloor{x}\rfloor^{2}"
    assert latex(ceiling(x)) == r"\lceil{x}\rceil"
    assert latex(ceiling(x)**2) == r"\lceil{x}\rceil^{2}"
    assert latex(Min(x, 2, x**3)) == r"\min\left(2, x, x^{3}\right)"
    assert latex(Min(x, y)**2) == r"\min\left(x, y\right)^{2}"
    assert latex(Max(x, 2, x**3)) == r"\max\left(2, x, x^{3}\right)"
    assert latex(Max(x, y)**2) == r"\max\left(x, y\right)^{2}"
    assert latex(Abs(x)) == r"\left|{x}\right|"
    assert latex(Abs(x)**3) == r'\left|{x}\right|^{3}'
    assert latex(re(x)) == r"\Re{x}"
    assert latex(re(x)**3) == r"\left(\Re{x}\right)^{3}"
    assert latex(re(x + y)) == r"\Re{x} + \Re{y}"
    assert latex(re(1/x, evaluate=False)) == r'\Re{\left(\frac{1}{x}\right)}'
    assert latex(im(x)) == r"\Im{x}"
    assert latex(im(x)**3) == r"\left(\Im{x}\right)^{3}"
    assert latex(im(1/x, evaluate=False)) == r'\Im{\left(\frac{1}{x}\right)}'
    assert latex(conjugate(x)) == r"\overline{x}"
    assert latex(conjugate(x)**3) == r'\overline{x}^{3}'
    assert latex(gamma(x)) == r"\Gamma{\left(x \right)}"
    assert latex(gamma(x)**3) == r"\Gamma^{3}{\left(x \right)}"
    w = Wild('w')
    assert latex(gamma(w)) == r"\Gamma{\left(w \right)}"
    assert latex(Order(x)) == r"\mathcal{O}\left(x\right)"
    assert latex(Order(x, x)) == r"\mathcal{O}\left(x\right)"
    assert latex(Order(x, (x, 0))) == r"\mathcal{O}\left(x\right)"
    assert latex(Order(x, (x, oo))) == r"\mathcal{O}\left(x; x\rightarrow{}\infty\right)"
    assert latex(Order(x, x, y)) == r"\mathcal{O}\left(x; \left ( x, \quad y\right )\rightarrow{}\left ( 0, \quad 0\right )\right)"
    assert latex(Order(x, x, y)) == r"\mathcal{O}\left(x; \left ( x, \quad y\right )\rightarrow{}\left ( 0, \quad 0\right )\right)"
    assert latex(Order(x, (x, oo), (y, oo))) == r"\mathcal{O}\left(x; \left ( x, \quad y\right )\rightarrow{}\left ( \infty, \quad \infty\right )\right)"
    assert latex(lowergamma(x, y)) == r'\gamma\left(x, y\right)'
    assert latex(lowergamma(x, y)**3) == r'\gamma^{3}\left(x, y\right)'
    assert latex(uppergamma(x, y)) == r'\Gamma\left(x, y\right)'
    assert latex(uppergamma(x, y)**3) == r'\Gamma^{3}\left(x, y\right)'

    assert latex(cot(x)) == r'\cot{\left (x \right )}'
    assert latex(coth(x)) == r'\coth{\left (x \right )}'
    assert latex(re(x)) == r'\Re{x}'
    assert latex(im(x)) == r'\Im{x}'
    assert latex(root(x, y)) == r'x^{\frac{1}{y}}'
    assert latex(arg(x)) == r'\arg{\left (x \right )}'
    assert latex(zeta(x)) == r'\zeta\left(x\right)'

    assert latex(zeta(x)) == r"\zeta\left(x\right)"
    assert latex(zeta(x)**2) == r"\zeta^{2}\left(x\right)"
    assert latex(zeta(x, y)) == r"\zeta\left(x, y\right)"
    assert latex(zeta(x, y)**2) == r"\zeta^{2}\left(x, y\right)"
    assert latex(dirichlet_eta(x)) == r"\eta\left(x\right)"
    assert latex(dirichlet_eta(x)**2) == r"\eta^{2}\left(x\right)"
    assert latex(polylog(x, y)) == r"\operatorname{Li}_{x}\left(y\right)"
    assert latex(
        polylog(x, y)**2) == r"\operatorname{Li}_{x}^{2}\left(y\right)"
    assert latex(lerchphi(x, y, n)) == r"\Phi\left(x, y, n\right)"
    assert latex(lerchphi(x, y, n)**2) == r"\Phi^{2}\left(x, y, n\right)"

    assert latex(elliptic_k(z)) == r"K\left(z\right)"
    assert latex(elliptic_k(z)**2) == r"K^{2}\left(z\right)"
    assert latex(elliptic_f(x, y)) == r"F\left(x\middle| y\right)"
    assert latex(elliptic_f(x, y)**2) == r"F^{2}\left(x\middle| y\right)"
    assert latex(elliptic_e(x, y)) == r"E\left(x\middle| y\right)"
    assert latex(elliptic_e(x, y)**2) == r"E^{2}\left(x\middle| y\right)"
    assert latex(elliptic_e(z)) == r"E\left(z\right)"
    assert latex(elliptic_e(z)**2) == r"E^{2}\left(z\right)"
    assert latex(elliptic_pi(x, y, z)) == r"\Pi\left(x; y\middle| z\right)"
    assert latex(elliptic_pi(x, y, z)**2) == \
        r"\Pi^{2}\left(x; y\middle| z\right)"
    assert latex(elliptic_pi(x, y)) == r"\Pi\left(x\middle| y\right)"
    assert latex(elliptic_pi(x, y)**2) == r"\Pi^{2}\left(x\middle| y\right)"

    assert latex(Ei(x)) == r'\operatorname{Ei}{\left (x \right )}'
    assert latex(Ei(x)**2) == r'\operatorname{Ei}^{2}{\left (x \right )}'
    assert latex(expint(x, y)**2) == r'\operatorname{E}_{x}^{2}\left(y\right)'
    assert latex(expint(x, y)) == r'\operatorname{E}_{x}\left(y\right)'
    assert latex(Shi(x)**2) == r'\operatorname{Shi}^{2}{\left (x \right )}'
    assert latex(Si(x)**2) == r'\operatorname{Si}^{2}{\left (x \right )}'
    assert latex(Ci(x)**2) == r'\operatorname{Ci}^{2}{\left (x \right )}'
    assert latex(Chi(x)**2) == r'\operatorname{Chi}^{2}{\left (x \right )}'
    assert latex(Chi(x)) == r'\operatorname{Chi}{\left (x \right )}'

    assert latex(
        jacobi(n, a, b, x)) == r'P_{n}^{\left(a,b\right)}\left(x\right)'
    assert latex(jacobi(n, a, b, x)**2) == r'\left(P_{n}^{\left(a,b\right)}\left(x\right)\right)^{2}'
    assert latex(
        gegenbauer(n, a, x)) == r'C_{n}^{\left(a\right)}\left(x\right)'
    assert latex(gegenbauer(n, a, x)**2) == r'\left(C_{n}^{\left(a\right)}\left(x\right)\right)^{2}'
    assert latex(chebyshevt(n, x)) == r'T_{n}\left(x\right)'
    assert latex(
        chebyshevt(n, x)**2) == r'\left(T_{n}\left(x\right)\right)^{2}'
    assert latex(chebyshevu(n, x)) == r'U_{n}\left(x\right)'
    assert latex(
        chebyshevu(n, x)**2) == r'\left(U_{n}\left(x\right)\right)^{2}'
    assert latex(legendre(n, x)) == r'P_{n}\left(x\right)'
    assert latex(legendre(n, x)**2) == r'\left(P_{n}\left(x\right)\right)^{2}'
    assert latex(
        assoc_legendre(n, a, x)) == r'P_{n}^{\left(a\right)}\left(x\right)'
    assert latex(assoc_legendre(n, a, x)**2) == r'\left(P_{n}^{\left(a\right)}\left(x\right)\right)^{2}'
    assert latex(laguerre(n, x)) == r'L_{n}\left(x\right)'
    assert latex(laguerre(n, x)**2) == r'\left(L_{n}\left(x\right)\right)^{2}'
    assert latex(
        assoc_laguerre(n, a, x)) == r'L_{n}^{\left(a\right)}\left(x\right)'
    assert latex(assoc_laguerre(n, a, x)**2) == r'\left(L_{n}^{\left(a\right)}\left(x\right)\right)^{2}'
    assert latex(hermite(n, x)) == r'H_{n}\left(x\right)'
    assert latex(hermite(n, x)**2) == r'\left(H_{n}\left(x\right)\right)^{2}'

    theta = Symbol("theta", extended_real=True)
    phi = Symbol("phi", extended_real=True)
    assert latex(Ynm(n, m, theta, phi)) == r'Y_{n}^{m}\left(\theta,\phi\right)'
    assert latex(Ynm(n, m, theta, phi)**3) == r'\left(Y_{n}^{m}\left(\theta,\phi\right)\right)^{3}'
    assert latex(Znm(n, m, theta, phi)) == r'Z_{n}^{m}\left(\theta,\phi\right)'
    assert latex(Znm(n, m, theta, phi)**3) == r'\left(Z_{n}^{m}\left(\theta,\phi\right)\right)^{3}'

    # Test latex printing of function names with "_"
    assert (latex(polar_lift(0, evaluate=False)) ==
            r"\operatorname{polar\_lift}{\left (0 \right )}")
    assert (latex(polar_lift(0, evaluate=False)**3) ==
            r"\operatorname{polar\_lift}^{3}{\left (0 \right )}")

    assert latex(totient(n)) == r'\phi\left(n\right)'
    assert latex(totient(n)) == r'\phi\left(n\right)'
    assert latex(totient(n)**2) == r'\left(\phi\left(n\right)\right)^{2}'

    assert latex(divisor_sigma(x)) == r"\sigma\left(x\right)"
    assert latex(divisor_sigma(x)**2) == r"\sigma^{2}\left(x\right)"
    assert latex(divisor_sigma(x, y)) == r"\sigma_y\left(x\right)"
    assert latex(divisor_sigma(x, y)**2) == r"\sigma^{2}_y\left(x\right)"

    assert latex(Mod(x, 7)) == r'x\bmod{7}'
    assert latex(Mod(x + 1, 7)) == r'\left(x + 1\right)\bmod{7}'
    assert latex(Mod(2 * x, 7)) == r'\left(2 x\right)\bmod{7}'
    assert latex(Mod(x, 7) + 1) == r'\left(x\bmod{7}\right) + 1'
    assert latex(2 * Mod(x, 7)) == r'2 \left(x\bmod{7}\right)'
    assert latex(Mod(x, 7)**2) == r'\left(x\bmod{7}\right)^{2}'

    # some unknown function name should get rendered with \operatorname
    fjlkd = Function('fjlkd')
    assert latex(fjlkd(x)) == r'\operatorname{fjlkd}{\left (x \right )}'
    # even when it is referred to without an argument
    assert latex(fjlkd) == r'\operatorname{fjlkd}'

    assert latex(airyai(x)) == r'Ai\left(x\right)'
    assert latex(airyai(x)**3) == r'Ai^{3}\left(x\right)'
    assert latex(airybi(x)) == r'Bi\left(x\right)'
    assert latex(airybi(x)**3) == r'Bi^{3}\left(x\right)'
    assert latex(airyaiprime(x)) == r'Ai^\prime\left(x\right)'
    assert latex(airyaiprime(x)**2) == r'{Ai^\prime}^{2}\left(x\right)'
    assert latex(airybiprime(x)) == r'Bi^\prime\left(x\right)'
    assert latex(airybiprime(x)**2) == r'{Bi^\prime}^{2}\left(x\right)'

    assert latex(euler(x)) == r'E_{x}'


def test_hyper_printing():
    assert latex(meijerg(Tuple(pi, pi, x), Tuple(1),
                         (0, 1), Tuple(1, 2, 3/pi), z)) == \
        r'{G_{4, 5}^{2, 3}\left(\begin{matrix} \pi, \pi, x & 1 \\0, 1 & 1, 2, \frac{3}{\pi} \end{matrix} \middle| {z} \right)}'
    assert latex(meijerg(Tuple(), Tuple(1), (0,), Tuple(), z)) == \
        r'{G_{1, 1}^{1, 0}\left(\begin{matrix}  & 1 \\0 &  \end{matrix} \middle| {z} \right)}'
    assert latex(latex(meijerg(Tuple(), Tuple(1), (0,), Tuple(), z)**2)) == \
        r'{{G_{1, 1}^{1, 0}\left(\begin{matrix}  & 1 \\0 &  ' \
        r'\end{matrix} \middle| {z} \right)}}^{2}'
    assert latex(hyper((x, 2), (3,), z)) == \
        r'{{}_{2}F_{1}\left(\begin{matrix} x, 2 ' \
        r'\\ 3 \end{matrix}\middle| {z} \right)}'
    assert latex(hyper((x, 2), (3,), z)**3) == \
        r'{{{}_{2}F_{1}\left(\begin{matrix}' \
        r' x, 2 \\ 3 \end{matrix}\middle| {z} \right)}}^{3}'
    assert latex(hyper(Tuple(), Tuple(1), z)) == \
        r'{{}_{0}F_{1}\left(\begin{matrix}  ' \
        r'\\ 1 \end{matrix}\middle| {z} \right)}'


def test_latex_bessel():
    assert latex(besselj(n, z**2)**k) == r'J^{k}_{n}\left(z^{2}\right)'
    assert latex(bessely(n, z)) == r'Y_{n}\left(z\right)'
    assert latex(besseli(n, z)) == r'I_{n}\left(z\right)'
    assert latex(besselk(n, z)) == r'K_{n}\left(z\right)'
    assert latex(hankel1(n, z**2)**2) == \
        r'\left(H^{(1)}_{n}\left(z^{2}\right)\right)^{2}'
    assert latex(hankel2(n, z)) == r'H^{(2)}_{n}\left(z\right)'
    assert latex(jn(n, z)) == r'j_{n}\left(z\right)'
    assert latex(yn(n, z)) == r'y_{n}\left(z\right)'


def test_latex_fresnel():
    assert latex(fresnels(z)) == r'S\left(z\right)'
    assert latex(fresnelc(z)) == r'C\left(z\right)'
    assert latex(fresnels(z)**2) == r'S^{2}\left(z\right)'
    assert latex(fresnelc(z)**2) == r'C^{2}\left(z\right)'


def test_latex_brackets():
    assert latex((-1)**x) == r"\left(-1\right)^{x}"


def test_latex_indexed():
    Psi_symbol = Symbol('Psi_0', complex=True, extended_real=False)
    Psi_indexed = IndexedBase(Symbol('Psi', complex=True, extended_real=False))
    symbol_latex = latex(Psi_symbol * conjugate(Psi_symbol))
    indexed_latex = latex(Psi_indexed[0] * conjugate(Psi_indexed[0]))
    # \\overline{\\Psi_{0}} \\Psi_{0}   vs.   \\Psi_{0} \\overline{\\Psi_{0}}
    assert symbol_latex.split() == indexed_latex.split() \
        or symbol_latex.split() == indexed_latex.split()[::-1]

    # Symbol('gamma') gives r'\gamma'
    assert latex(IndexedBase('gamma')) == r'\gamma'
    assert latex(IndexedBase('a b')) == 'a b'
    assert latex(IndexedBase('a_b')) == 'a_{b}'


def test_latex_derivatives():
    # regular "d" for ordinary derivatives
    assert latex(diff(x**3, x, evaluate=False)) == \
        r"\frac{d}{d x} x^{3}"
    assert latex(diff(sin(x) + x**2, x, evaluate=False)) == \
        r"\frac{d}{d x}\left(x^{2} + \sin{\left (x \right )}\right)"
    assert latex(diff(diff(sin(x) + x**2, x, evaluate=False), evaluate=False)) == \
        r"\frac{d^{2}}{d x^{2}} \left(x^{2} + \sin{\left (x \right )}\right)"
    assert latex(diff(diff(diff(sin(x) + x**2, x, evaluate=False), evaluate=False), evaluate=False)) == \
        r"\frac{d^{3}}{d x^{3}} \left(x^{2} + \sin{\left (x \right )}\right)"

    # \partial for partial derivatives
    assert latex(diff(sin(x * y), x, evaluate=False)) == \
        r"\frac{\partial}{\partial x} \sin{\left (x y \right )}"
    assert latex(diff(sin(x * y) + x**2, x, evaluate=False)) == \
        r"\frac{\partial}{\partial x}\left(x^{2} + \sin{\left (x y \right )}\right)"
    assert latex(diff(diff(sin(x*y) + x**2, x, evaluate=False), x, evaluate=False)) == \
        r"\frac{\partial^{2}}{\partial x^{2}} \left(x^{2} + \sin{\left (x y \right )}\right)"
    assert latex(diff(diff(diff(sin(x*y) + x**2, x, evaluate=False), x, evaluate=False), x, evaluate=False)) == \
        r"\frac{\partial^{3}}{\partial x^{3}} \left(x^{2} + \sin{\left (x y \right )}\right)"

    # mixed partial derivatives
    f = Function("f")
    assert latex(diff(diff(f(x, y), x, evaluate=False), y, evaluate=False)) == \
        r"\frac{\partial^{2}}{\partial x\partial y}  " + latex(f(x, y))

    assert latex(diff(diff(diff(f(x, y), x, evaluate=False), x, evaluate=False), y, evaluate=False)) == \
        r"\frac{\partial^{3}}{\partial x^{2}\partial y}  " + latex(f(x, y))

    # use ordinary d when one of the variables has been integrated out
    assert latex(diff(Integral(exp(-x * y), (x, 0, oo)), y, evaluate=False)) == \
        r"\frac{d}{d y} \int_{0}^{\infty} e^{- x y}\, dx"


def test_latex_subs():
    assert latex(Subs(x, (x, 1))) == r'\left. x \right|_{\substack{ x=1 }}'
    assert latex(Subs(x*y, (x, 1),
                      (y, 2))) == r'\left. x y \right|_{\substack{ x=1\\ y=2 }}'


def test_latex_integrals():
    assert latex(Integral(log(x), x)) == r"\int \log{\left (x \right )}\, dx"
    assert latex(Integral(x**2, (x, 0, 1))) == r"\int_{0}^{1} x^{2}\, dx"
    assert latex(Integral(x**2, (x, 10, 20))) == r"\int_{10}^{20} x^{2}\, dx"
    assert latex(Integral(
        y*x**2, (x, 0, 1), y)) == r"\int\int_{0}^{1} x^{2} y\, dx\, dy"
    assert latex(Integral(y*x**2, (x, 0, 1), y), mode='equation*') \
        == r"\begin{equation*}\int\int\limits_{0}^{1} x^{2} y\, dx\, dy\end{equation*}"
    assert latex(Integral(y*x**2, (x, 0, 1), y), mode='equation*', itex=True) \
        == r"$$\int\int_{0}^{1} x^{2} y\, dx\, dy$$"
    assert latex(Integral(x, (x, 0))) == r"\int^{0} x\, dx"
    assert latex(Integral(x*y, x, y)) == r"\iint x y\, dx\, dy"
    assert latex(Integral(x*y*z, x, y, z)) == r"\iiint x y z\, dx\, dy\, dz"
    assert latex(Integral(x*y*z*t, x, y, z, t)) == \
        r"\iiiint t x y z\, dx\, dy\, dz\, dt"
    assert latex(Integral(x, x, x, x, x, x, x)) == \
        r"\int\int\int\int\int\int x\, dx\, dx\, dx\, dx\, dx\, dx"
    assert latex(Integral(x, x, y, (z, 0, 1))) == \
        r"\int_{0}^{1}\int\int x\, dx\, dy\, dz"

    # issue sympy/sympy#10806
    assert latex(Integral(x, x)**2) == r"\left(\int x\, dx\right)^{2}"
    assert latex(Integral(z + x, z)) == r"\int \left(x + z\right)\, dz"


def test_latex_sets():
    for s in (frozenset, set):
        assert latex(s([x*y, x**2])) == r"\left\{x^{2}, x y\right\}"
        assert latex(s(range(1, 6))) == r"\left\{1, 2, 3, 4, 5\right\}"
        assert latex(s(range(1, 13))) == \
            r"\left\{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12\right\}"

    s = FiniteSet
    assert latex(s(*[x*y, x**2])) == r"\left\{x^{2}, x y\right\}"
    assert latex(s(*range(1, 6))) == r"\left\{1, 2, 3, 4, 5\right\}"
    assert latex(s(*range(1, 13))) == \
        r"\left\{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12\right\}"

    assert (latex(Intersection(s(x), Interval(y, 2))) ==
            r'\left\{x\right\} \cap \left[y, 2\right]')


def test_latex_Range():
    assert latex(Range(1, 51)) == \
        r'\left\{1, 2, \ldots, 50\right\}'
    assert latex(Range(1, 4)) == r'\left\{1, 2, 3\right\}'


def test_latex_intervals():
    a = Symbol('a', extended_real=True)
    assert latex(Interval(0, 0)) == r"\left\{0\right\}"
    assert latex(Interval(0, a)) == r"\left[0, a\right]"
    assert latex(Interval(0, a, False, False)) == r"\left[0, a\right]"
    assert latex(Interval(0, a, True, False)) == r"\left(0, a\right]"
    assert latex(Interval(0, a, False, True)) == r"\left[0, a\right)"
    assert latex(Interval(0, a, True, True)) == r"\left(0, a\right)"


def test_latex_emptyset():
    assert latex(S.EmptySet) == r"\emptyset"


def test_latex_union():
    assert latex(Union(Interval(0, 1), Interval(2, 3))) == \
        r"\left[0, 1\right] \cup \left[2, 3\right]"
    assert latex(Union(Interval(1, 1), Interval(2, 2), Interval(3, 4))) == \
        r"\left\{1, 2\right\} \cup \left[3, 4\right]"


def test_latex_symmetric_difference():
    assert latex(SymmetricDifference(Interval(2, 5), Interval(4, 7),
                                     evaluate=False)) == r'\left[2, 5\right] \triangle \left[4, 7\right]'


def test_latex_Complement():
    assert latex(Complement(S.Reals, S.Naturals)) == r"\mathbb{R} \setminus \mathbb{N}"


def test_latex_productset():
    line = Interval(0, 1)
    bigline = Interval(0, 10)
    fset = FiniteSet(1, 2, 3)
    assert latex(line**2) == r"%s^2" % latex(line)
    assert latex(line * bigline * fset) == r"%s \times %s \times %s" % (
        latex(line), latex(bigline), latex(fset))


def test_latex_Naturals():
    assert latex(S.Naturals) == r"\mathbb{N}"


def test_latex_Naturals0():
    assert latex(S.Naturals0) == r"\mathbb{N}_0"


def test_latex_Integers():
    assert latex(S.Integers) == r"\mathbb{Z}"


def test_latex_Rationals():
    assert latex(S.Rationals) == r"\mathbb{Q}"


def test_latex_ImageSet():
    x = Symbol('x')
    assert latex(ImageSet(Lambda(x, x**2), S.Naturals)) == \
        r"\left\{x^{2}\; |\; x \in \mathbb{N}\right\}"


def test_latex_Contains():
    x = Symbol('x')
    assert latex(Contains(x, S.Naturals)) == r"x \in \mathbb{N}"


def test_latex_sum():
    assert latex(Sum(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"\sum_{\substack{-2 \leq x \leq 2\\-5 \leq y \leq 5}} x y^{2}"
    assert latex(Sum(x**2, (x, -2, 2))) == \
        r"\sum_{x=-2}^{2} x^{2}"
    assert latex(Sum(x**2 + y, (x, -2, 2))) == \
        r"\sum_{x=-2}^{2} \left(x^{2} + y\right)"


def test_latex_product():
    assert latex(Product(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"\prod_{\substack{-2 \leq x \leq 2\\-5 \leq y \leq 5}} x y^{2}"
    assert latex(Product(x**2, (x, -2, 2))) == \
        r"\prod_{x=-2}^{2} x^{2}"
    assert latex(Product(x**2 + y, (x, -2, 2))) == \
        r"\prod_{x=-2}^{2} \left(x^{2} + y\right)"


def test_latex_limits():
    assert latex(Limit(x, x, oo)) == r"\lim_{x \to \infty} x"

    # issue sympy/sympy#8175
    f = Function('f')
    assert latex(Limit(f(x), x, 0)) == r"\lim_{x \to 0^+} f{\left (x \right )}"
    assert latex(Limit(f(x), x, 0, "-")) == r"\lim_{x \to 0^-} f{\left (x \right )}"
    assert latex(Limit(f(x), x, 0, "real")) == r"\lim_{x \to 0} f{\left (x \right )}"

    assert latex(Limit(x + x**2, x, 0)) == r"\lim_{x \to 0^+}\left(x^{2} + x\right)"
    assert latex(Limit(x > 0, x, 0)) == r"\lim_{x \to 0^+}\left(x > 0\right)"


def test_sympyissue_3568():
    beta = Symbol(r'\beta')
    y = beta + x
    assert latex(y) in [r'\beta + x', r'x + \beta']

    beta = Symbol(r'beta')
    y = beta + x
    assert latex(y) in [r'\beta + x', r'x + \beta']


def test_latex():
    assert latex((2*tau)**Rational(7, 2)) == "8 \\sqrt{2} \\tau^{\\frac{7}{2}}"
    assert latex((2*mu)**Rational(7, 2), mode='equation*') == \
        "\\begin{equation*}8 \\sqrt{2} \\mu^{\\frac{7}{2}}\\end{equation*}"
    assert latex((2*mu)**Rational(7, 2), mode='equation', itex=True) == \
        "$$8 \\sqrt{2} \\mu^{\\frac{7}{2}}$$"
    assert latex([2/x, y]) == r"\left [ \frac{2}{x}, \quad y\right ]"


def test_latex_dict():
    d = {Integer(1): 1, x**2: 2, x: 3, x**3: 4}
    assert latex(d) == r'\left \{ 1 : 1, \quad x : 3, \quad x^{2} : 2, \quad x^{3} : 4\right \}'
    D = Dict(d)
    assert latex(D) == r'\left \{ 1 : 1, \quad x : 3, \quad x^{2} : 2, \quad x^{3} : 4\right \}'


def test_latex_list():
    l = [Symbol('omega1'), Symbol('a'), Symbol('alpha')]
    assert latex(l) == r'\left [ \omega_{1}, \quad a, \quad \alpha\right ]'


def test_latex_rational():
    # tests issue sympy/sympy#3973
    assert latex(-Rational(1, 2)) == "- \\frac{1}{2}"
    assert latex(Rational(-1, 2)) == "- \\frac{1}{2}"
    assert latex(Rational(1, -2)) == "- \\frac{1}{2}"
    assert latex(-Rational(-1, 2)) == "\\frac{1}{2}"
    assert latex(-Rational(1, 2)*x) == "- \\frac{x}{2}"
    assert latex(-Rational(1, 2)*x + Rational(-2, 3)*y) == \
        "- \\frac{x}{2} - \\frac{2 y}{3}"


def test_latex_inverse():
    # tests issue sympy/sympy#4129
    assert latex(1/x) == "\\frac{1}{x}"
    assert latex(1/(x + y)) == "\\frac{1}{x + y}"


def test_latex_DiracDelta():
    assert latex(DiracDelta(x)) == r"\delta\left(x\right)"
    assert latex(DiracDelta(x)**2) == r"\left(\delta\left(x\right)\right)^{2}"
    assert latex(DiracDelta(x, 0)) == r"\delta\left(x\right)"
    assert latex(DiracDelta(x, 5)) == \
        r"\delta^{\left( 5 \right)}\left( x \right)"
    assert latex(DiracDelta(x, 5)**2) == \
        r"\left(\delta^{\left( 5 \right)}\left( x \right)\right)^{2}"


def test_latex_Heaviside():
    assert latex(Heaviside(x)) == r"\theta\left(x\right)"
    assert latex(Heaviside(x)**2) == r"\left(\theta\left(x\right)\right)^{2}"


def test_latex_KroneckerDelta():
    assert latex(KroneckerDelta(x, y)) == r"\delta_{x y}"
    assert latex(KroneckerDelta(x, y + 1)) == r"\delta_{x, y + 1}"
    # issue sympy/sympy#6578
    assert latex(KroneckerDelta(x + 1, y)) == r"\delta_{y, x + 1}"

    assert (latex(Pow(KroneckerDelta(x, y), 2, evaluate=False)) ==
            r'\left(\delta_{x y}\right)^{2}')


def test_latex_LeviCivita():
    assert latex(LeviCivita(x, y, z)) == r"\varepsilon_{x y z}"
    assert latex(LeviCivita(x, y, z)**2) == r"\left(\varepsilon_{x y z}\right)^{2}"
    assert latex(LeviCivita(x, y, z + 1)) == r"\varepsilon_{x, y, z + 1}"
    assert latex(LeviCivita(x, y + 1, z)) == r"\varepsilon_{x, y + 1, z}"
    assert latex(LeviCivita(x + 1, y, z)) == r"\varepsilon_{x + 1, y, z}"


def test_mode():
    expr = x + y
    assert latex(expr) == 'x + y'
    assert latex(expr, mode='plain') == 'x + y'
    assert latex(expr, mode='inline') == '$x + y$'
    assert latex(
        expr, mode='equation*') == '\\begin{equation*}x + y\\end{equation*}'
    assert latex(
        expr, mode='equation') == '\\begin{equation}x + y\\end{equation}'

    pytest.raises(ValueError, lambda: latex(expr, mode='spam'))


def test_latex_Piecewise():
    p = Piecewise((x, x < 1), (x**2, True))
    assert latex(p) == r"\begin{cases} x & \text{for}\: x < 1 \\x^{2} &" \
                       r" \text{otherwise} \end{cases}"
    assert latex(p, itex=True) == r"\begin{cases} x & \text{for}\: x \lt 1 \\x^{2} &" \
                                  r" \text{otherwise} \end{cases}"
    p = Piecewise((x, x < 0), (0, x >= 0))
    assert latex(p) == r"\begin{cases} x & \text{for}\: x < 0 \\0 &" \
                       r" \text{for}\: x \geq 0 \end{cases}"
    A, B = symbols("A B", commutative=False)
    p = Piecewise((A**2, Eq(A, B)), (A*B, True))
    s = r"\begin{cases} A^{2} & \text{for}\: A = B \\A B & \text{otherwise} \end{cases}"
    assert latex(p) == s
    assert latex(A*p) == r"A \left(%s\right)" % s
    assert latex(p*A) == r"\left(%s\right) A" % s


def test_latex_Matrix():
    M = Matrix([[1 + x, y], [y, x - 1]])
    assert latex(M) == \
        r'\left[\begin{matrix}x + 1 & y\\y & x - 1\end{matrix}\right]'
    assert latex(M, mode='inline') == \
        r'$\left[\begin{smallmatrix}x + 1 & y\\' \
        r'y & x - 1\end{smallmatrix}\right]$'
    assert latex(M, mat_str='array') == \
        r'\left[\begin{array}{cc}x + 1 & y\\y & x - 1\end{array}\right]'
    assert latex(M, mat_str='bmatrix') == \
        r'\left[\begin{bmatrix}x + 1 & y\\y & x - 1\end{bmatrix}\right]'
    assert latex(M, mat_delim=None, mat_str='bmatrix') == \
        r'\begin{bmatrix}x + 1 & y\\y & x - 1\end{bmatrix}'
    M2 = Matrix(1, 11, range(11))
    assert latex(M2) == \
        r'\left[\begin{array}{ccccccccccc}' \
        r'0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10\end{array}\right]'

    n, m, l = symbols('n,m,l')
    X = MatrixSymbol('X', n, n)
    Y = MatrixSymbol('Y', m, m)
    Z = MatrixSymbol('Z', n, m)
    B = BlockMatrix([[X, Z], [ZeroMatrix(m, n), Y]])
    assert latex(B) == r'\left[\begin{matrix}X & Z' \
                       r'\\\mathbb{0} & Y\end{matrix}\right]'

    assert latex(Identity(3)) == r'\mathbb{I}'


def test_latex_matrix_with_functions():
    t = symbols('t')
    theta1 = symbols('theta1', cls=Function)

    M = Matrix([[sin(theta1(t)), cos(theta1(t))],
                [cos(theta1(t).diff(t)), sin(theta1(t).diff(t))]])

    expected = (r'\left[\begin{matrix}\sin{\left '
                r'(\theta_{1}{\left (t \right )} \right )} & '
                r'\cos{\left (\theta_{1}{\left (t \right )} \right '
                r')}\\\cos{\left (\frac{d}{d t} \theta_{1}{\left (t '
                r'\right )} \right )} & \sin{\left (\frac{d}{d t} '
                r'\theta_{1}{\left (t \right )} \right '
                r')}\end{matrix}\right]')

    assert latex(M) == expected


def test_latex_NDimArray():
    for ArrayType in (ImmutableDenseNDimArray, ImmutableSparseNDimArray,
                      MutableDenseNDimArray, MutableSparseNDimArray):
        M = ArrayType([[1 / x, y], [z, w]])
        M1 = ArrayType([1 / x, y, z])

        M2 = tensorproduct(M1, M)
        M3 = tensorproduct(M, M)

        assert latex(M) == latex(M, mat_str='matrix') == r"\left[\begin{matrix}\frac{1}{x} & y\\z & w\end{matrix}\right]"
        assert latex(M, mat_delim=None) == r"\begin{matrix}\frac{1}{x} & y\\z & w\end{matrix}"
        assert latex(M, mode='inline') == r"$\left[\begin{smallmatrix}1 / x & y\\z & w\end{smallmatrix}\right]$"
        assert latex(M1) == r"\left[\begin{matrix}\frac{1}{x} & y & z\end{matrix}\right]"
        assert latex(M2) == r"\left[\begin{matrix}" \
            r"\left[\begin{matrix}\frac{1}{x^{2}} & \frac{y}{x}\\\frac{z}{x} & \frac{w}{x}\end{matrix}\right] & " \
            r"\left[\begin{matrix}\frac{y}{x} & y^{2}\\y z & w y\end{matrix}\right] & " \
            r"\left[\begin{matrix}\frac{z}{x} & y z\\z^{2} & w z\end{matrix}\right]" \
            r"\end{matrix}\right]"
        assert latex(M3) == r"""\left[\begin{matrix}"""\
            r"""\left[\begin{matrix}\frac{1}{x^{2}} & \frac{y}{x}\\\frac{z}{x} & \frac{w}{x}\end{matrix}\right] & """\
            r"""\left[\begin{matrix}\frac{y}{x} & y^{2}\\y z & w y\end{matrix}\right]\\"""\
            r"""\left[\begin{matrix}\frac{z}{x} & y z\\z^{2} & w z\end{matrix}\right] & """\
            r"""\left[\begin{matrix}\frac{w}{x} & w y\\w z & w^{2}\end{matrix}\right]"""\
            r"""\end{matrix}\right]"""
        assert latex(ArrayType()) == r"\left[\begin{matrix}\end{matrix}\right]"

        Mrow = ArrayType([[x, y, 1/z]])
        Mcolumn = ArrayType([[x], [y], [1/z]])
        Mcol2 = ArrayType([Mcolumn.tolist()])

        assert latex(Mrow) == r"\left[\left[\begin{matrix}x & y & \frac{1}{z}\end{matrix}\right]\right]"
        assert latex(Mcolumn) == r"\left[\begin{matrix}x\\y\\\frac{1}{z}\end{matrix}\right]"
        assert latex(Mcol2) == r'\left[\begin{matrix}\left[\begin{matrix}x\\y\\\frac{1}{z}\end{matrix}\right]\end{matrix}\right]'

        Arow = ArrayType([x]*11)

        assert latex(Arow) == r"\left[\begin{array}x & x & x & x & x & x & x & x & x & x & x\end{array}\right]"


def test_latex_mul_symbol():
    assert latex(4*4**x, mul_symbol='times') == "4 \\times 4^{x}"
    assert latex(4*4**x, mul_symbol='dot') == "4 \\cdot 4^{x}"
    assert latex(4*4**x, mul_symbol='ldot') == r"4 \,.\, 4^{x}"

    assert latex(4*x, mul_symbol='times') == "4 \\times x"
    assert latex(4*x, mul_symbol='dot') == "4 \\cdot x"
    assert latex(4*x, mul_symbol='ldot') == r"4 \,.\, x"


def test_latex_sympyissue_4381():
    y = 4*4**log(2)
    assert latex(y) == r'4 \cdot 4^{\log{\left (2 \right )}}'
    assert latex(1/y) == r'\frac{1}{4 \cdot 4^{\log{\left (2 \right )}}}'


def test_latex_sympyissue_4576():
    assert latex(Symbol("beta_13_2")) == r"\beta_{13 2}"
    assert latex(Symbol("beta_132_20")) == r"\beta_{132 20}"
    assert latex(Symbol("beta_13")) == r"\beta_{13}"
    assert latex(Symbol("x_a_b")) == r"x_{a b}"
    assert latex(Symbol("x_1_2_3")) == r"x_{1 2 3}"
    assert latex(Symbol("x_a_b1")) == r"x_{a b1}"
    assert latex(Symbol("x_a_1")) == r"x_{a 1}"
    assert latex(Symbol("x_1_a")) == r"x_{1 a}"
    assert latex(Symbol("x_1^aa")) == r"x^{aa}_{1}"
    assert latex(Symbol("x_1__aa")) == r"x^{aa}_{1}"
    assert latex(Symbol("x_11^a")) == r"x^{a}_{11}"
    assert latex(Symbol("x_11__a")) == r"x^{a}_{11}"
    assert latex(Symbol("x_a_a_a_a")) == r"x_{a a a a}"
    assert latex(Symbol("x_a_a^a^a")) == r"x^{a a}_{a a}"
    assert latex(Symbol("x_a_a__a__a")) == r"x^{a a}_{a a}"
    assert latex(Symbol("alpha_11")) == r"\alpha_{11}"
    assert latex(Symbol("alpha_11_11")) == r"\alpha_{11 11}"
    assert latex(Symbol("alpha_alpha")) == r"\alpha_{\alpha}"
    assert latex(Symbol("alpha^aleph")) == r"\alpha^{\aleph}"
    assert latex(Symbol("alpha__aleph")) == r"\alpha^{\aleph}"


def test_latex_pow_fraction():
    x = Symbol('x')
    # Testing exp
    assert 'e^{-x}' in latex(exp(-x)/2).replace(' ', '')  # Remove Whitespace

    # Testing just e^{-x} in case future changes alter behavior of muls or fracs
    # In particular current output is \frac{1}{2}e^{- x} but perhaps this will
    # change to \frac{e^{-x}}{2}

    # Testing general, non-exp, power
    assert '3^{-x}' in latex(3**-x/2).replace(' ', '')


def test_noncommutative():
    A, B, C = symbols('A,B,C', commutative=False)

    assert latex(A*B*C**-1) == "A B C^{-1}"
    assert latex(C**-1*A*B) == "C^{-1} A B"
    assert latex(A*C**-1*B) == "A C^{-1} B"


def test_latex_order():
    expr = x**3 + x**2*y + 3*x*y**3 + y**4

    assert latex(expr, order='lex') == "x^{3} + x^{2} y + 3 x y^{3} + y^{4}"
    assert latex(
        expr, order='rev-lex') == "y^{4} + 3 x y^{3} + x^{2} y + x^{3}"


def test_latex_Lambda():
    assert latex(Lambda(x, x + 1)) == \
        r"\left( x \mapsto x + 1 \right)"
    assert latex(Lambda((x, y), x + 1)) == \
        r"\left( \left ( x, \quad y\right ) \mapsto x + 1 \right)"


def test_latex_PolyElement():
    Ruv,  u, v = ring("u,v", ZZ)
    Rxyz,  x, y, z = ring("x,y,z", Ruv)

    assert latex(x - x) == r"0"
    assert latex(x - 1) == r"x - 1"
    assert latex(x + 1) == r"x + 1"

    assert latex((u**2 + 3*u*v + 1)*x**2*y + u + 1) == r"\left({u}^{2} + 3 u v + 1\right) {x}^{2} y + u + 1"
    assert latex((u**2 + 3*u*v + 1)*x**2*y + (u + 1)*x) == r"\left({u}^{2} + 3 u v + 1\right) {x}^{2} y + \left(u + 1\right) x"
    assert latex((u**2 + 3*u*v + 1)*x**2*y + (u + 1)*x + 1) == r"\left({u}^{2} + 3 u v + 1\right) {x}^{2} y + \left(u + 1\right) x + 1"
    assert latex((-u**2 + 3*u*v - 1)*x**2*y - (u + 1)*x - 1) == r"-\left({u}^{2} - 3 u v + 1\right) {x}^{2} y - \left(u + 1\right) x - 1"

    assert latex(-(v**2 + v + 1)*x + 3*u*v + 1) == r"-\left({v}^{2} + v + 1\right) x + 3 u v + 1"
    assert latex(-(v**2 + v + 1)*x - 3*u*v + 1) == r"-\left({v}^{2} + v + 1\right) x - 3 u v + 1"


def test_latex_FracElement():
    Fuv,  u, v = field("u,v", ZZ)
    Fxyzt,  x, y, z, t = field("x,y,z,t", Fuv)

    assert latex(x - x) == r"0"
    assert latex(x - 1) == r"x - 1"
    assert latex(x + 1) == r"x + 1"

    assert latex(x/3) == r"\frac{x}{3}"
    assert latex(x/z) == r"\frac{x}{z}"
    assert latex(x*y/z) == r"\frac{x y}{z}"
    assert latex(x/(z*t)) == r"\frac{x}{z t}"
    assert latex(x*y/(z*t)) == r"\frac{x y}{z t}"

    assert latex((x - 1)/y) == r"\frac{x - 1}{y}"
    assert latex((x + 1)/y) == r"\frac{x + 1}{y}"
    assert latex((-x - 1)/y) == r"\frac{-x - 1}{y}"
    assert latex((x + 1)/(y*z)) == r"\frac{x + 1}{y z}"
    assert latex(-y/(x + 1)) == r"\frac{-y}{x + 1}"
    assert latex(y*z/(x + 1)) == r"\frac{y z}{x + 1}"

    assert latex(((u + 1)*x*y + 1)/((v - 1)*z - 1)) == r"\frac{\left(u + 1\right) x y + 1}{\left(v - 1\right) z - 1}"
    assert latex(((u + 1)*x*y + 1)/((v - 1)*z - t*u*v - 1)) == r"\frac{\left(u + 1\right) x y + 1}{\left(v - 1\right) z - u v t - 1}"


def test_latex_Poly():
    assert latex(Poly(x**2 + 2 * x, x)) == \
        r"\operatorname{Poly}{\left( x^{2} + 2 x, x, domain=\mathbb{Z} \right)}"
    assert latex(Poly(x/y, x)) == \
        r"\operatorname{Poly}{\left( \frac{x}{y}, x, domain=\mathbb{Z}\left(y\right) \right)}"
    assert latex(Poly(2.0*x + y)) == \
        r"\operatorname{Poly}{\left( 2.0 x + 1.0 y, x, y, domain=\mathbb{R} \right)}"


def test_latex_domains():
    assert latex(FF(2)) == r'\mathbb{F}_{2}'
    assert latex(CC) == r'\mathbb{C}'


def test_latex_RootOf():
    assert latex(RootOf(x**5 + x + 3, 0)) == \
        r"\operatorname{RootOf} {\left(x^{5} + x + 3, 0\right)}"
    assert latex(RootOf(x**5 + x + y, x, 0)) == \
        r"\operatorname{RootOf} {\left(x^{5} + x + y, x, 0\right)}"


def test_latex_RootSum():
    assert latex(RootSum(x**5 + x + 3, sin)) == \
        r"\operatorname{RootSum} {\left(x^{5} + x + 3, \left( x \mapsto \sin{\left (x \right )} \right)\right)}"


def test_settings():
    pytest.raises(TypeError, lambda: latex(x*y, method="garbage"))


def test_latex_numbers():
    assert latex(catalan(n)) == r"C_{n}"


def test_lamda():
    assert latex(Symbol('lamda')) == r"\lambda"
    assert latex(Symbol('Lamda')) == r"\Lambda"


def test_custom_symbol_names():
    x = Symbol('x')
    y = Symbol('y')
    assert latex(x) == "x"
    assert latex(x, symbol_names={x: "x_i"}) == "x_i"
    assert latex(x + y, symbol_names={x: "x_i"}) == "x_i + y"
    assert latex(x**2, symbol_names={x: "x_i"}) == "x_i^{2}"
    assert latex(x + y, symbol_names={x: "x_i", y: "y_j"}) == "x_i + y_j"


def test_matAdd():
    C = MatrixSymbol('C', 5, 5)
    B = MatrixSymbol('B', 5, 5)
    l = LatexPrinter()
    assert l._print_MatAdd(C - 2*B) == 'C - 2 B'
    assert l._print_MatAdd(C + 2*B) in ['2 B + C', 'C + 2 B']
    assert l._print_MatAdd(B - 2*C) == 'B - 2 C'
    assert l._print_MatAdd(B + 2*C) in ['B + 2 C', '2 C + B']


def test_matMul():
    A = MatrixSymbol('A', 5, 5)
    B = MatrixSymbol('B', 5, 5)
    x = Symbol('x')
    l = LatexPrinter()
    assert l._print_MatMul(2*A) == '2 A'
    assert l._print_MatMul(2*x*A) == '2 x A'
    assert l._print_MatMul(-2*A) == '-2 A'
    assert l._print_MatMul(1.5*A) == '1.5 A'
    assert l._print_MatMul(sqrt(2)*A) == r'\sqrt{2} A'
    assert l._print_MatMul(-sqrt(2)*A) == r'- \sqrt{2} A'
    assert l._print_MatMul(2*sqrt(2)*x*A) == r'2 \sqrt{2} x A'
    assert l._print_MatMul(-2*A*(A + 2*B)) in [r'-2 A \left(A + 2 B\right)',
                                               r'-2 A \left(2 B + A\right)']


def test_latex_MatrixSlice():
    assert latex(MatrixSymbol('X', 10, 10)[:5, 1:9:2]) == \
        r'X\left[:5, 1:9:2\right]'
    assert latex(MatrixSymbol('X', 10, 10)[5, :5:2]) == \
        r'X\left[5, :5:2\right]'


def test_latex_RandomDomain():
    X = Normal('x1', 0, 1)
    assert latex(where(X > 0)) == r"Domain: 0 < x_{1} \wedge x_{1} < \infty"

    D = Die('d1', 6)
    assert latex(where(D > 4)) == r"Domain: d_{1} = 5 \vee d_{1} = 6"

    A = Exponential('a', 1)
    B = Exponential('b', 1)
    assert latex(
        pspace(Tuple(A, B)).domain) == \
        r"Domain: 0 \leq a \wedge 0 \leq b \wedge a < \infty \wedge b < \infty"


def test_PrettyPoly():
    F = QQ.frac_field(x, y)
    R = QQ.poly_ring(x, y)

    assert latex(F.convert(x/(x + y))) == latex(x/(x + y))
    assert latex(R.convert(x + y)) == latex(x + y)


def test_integral_transforms():
    x = Symbol("x")
    k = Symbol("k")
    f = Function("f")
    a = Symbol("a")
    b = Symbol("b")

    assert latex(MellinTransform(f(x), x, k)) == r"\mathcal{M}_{x}\left[f{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseMellinTransform(f(k), k, x, a, b)) == r"\mathcal{M}^{-1}_{k}\left[f{\left (k \right )}\right]\left(x\right)"

    assert latex(LaplaceTransform(f(x), x, k)) == r"\mathcal{L}_{x}\left[f{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseLaplaceTransform(f(k), k, x, (a, b))) == r"\mathcal{L}^{-1}_{k}\left[f{\left (k \right )}\right]\left(x\right)"

    assert latex(FourierTransform(f(x), x, k)) == r"\mathcal{F}_{x}\left[f{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseFourierTransform(f(k), k, x)) == r"\mathcal{F}^{-1}_{k}\left[f{\left (k \right )}\right]\left(x\right)"

    assert latex(CosineTransform(f(x), x, k)) == r"\mathcal{COS}_{x}\left[f{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseCosineTransform(f(k), k, x)) == r"\mathcal{COS}^{-1}_{k}\left[f{\left (k \right )}\right]\left(x\right)"

    assert latex(SineTransform(f(x), x, k)) == r"\mathcal{SIN}_{x}\left[f{\left (x \right )}\right]\left(k\right)"
    assert latex(InverseSineTransform(f(k), k, x)) == r"\mathcal{SIN}^{-1}_{k}\left[f{\left (k \right )}\right]\left(x\right)"


def test_PolynomialRing():
    assert latex(QQ.poly_ring(x, y)) == r"\mathbb{Q}\left[x, y\right]"


def test_Tr():
    # TODO: Handle indices
    A, B = symbols('A B', commutative=False)
    t = Tr(A*B)
    assert latex(t) == r'\mbox{Tr}\left(A B\right)'


def test_Adjoint():
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    assert latex(Adjoint(X)) == r'X^\dag'
    assert latex(Adjoint(X + Y)) == r'\left(X + Y\right)^\dag'
    assert latex(Adjoint(X) + Adjoint(Y)) == r'X^\dag + Y^\dag'
    assert latex(Adjoint(X*Y)) == r'\left(X Y\right)^\dag'
    assert latex(Adjoint(Y)*Adjoint(X)) == r'Y^\dag X^\dag'
    assert latex(Adjoint(X**2)) == r'\left(X^{2}\right)^\dag'
    assert latex(Adjoint(X)**2) == r'\left(X^\dag\right)^{2}'
    assert latex(Adjoint(Inverse(X))) == r'\left(X^{-1}\right)^\dag'
    assert latex(Inverse(Adjoint(X))) == r'\left(X^\dag\right)^{-1}'
    assert latex(Adjoint(Transpose(X))) == r'\left(X^T\right)^\dag'
    assert latex(Transpose(Adjoint(X))) == r'\left(X^\dag\right)^T'


def test_Hadamard():
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    assert latex(HadamardProduct(X, Y*Y)) == r'X \circ \left(Y Y\right)'
    assert latex(HadamardProduct(X, Y)*Y) == r'\left(X \circ Y\right) Y'


def test_ZeroMatrix():
    assert latex(ZeroMatrix(1, 1)) == r"\mathbb{0}"


def test_boolean_args_order():
    syms = symbols('a:f')

    expr = And(*syms)
    assert latex(expr) == 'a \\wedge b \\wedge c \\wedge d \\wedge e \\wedge f'

    expr = Or(*syms)
    assert latex(expr) == 'a \\vee b \\vee c \\vee d \\vee e \\vee f'

    expr = Equivalent(*syms)
    assert latex(expr) == 'a \\equiv b \\equiv c \\equiv d \\equiv e \\equiv f'

    expr = Xor(*syms)
    assert latex(expr) == 'a \\veebar b \\veebar c \\veebar d \\veebar e \\veebar f'


def test_booleans():
    assert latex(Not(And(x, y))) == r'\neg (x \wedge y)'


def test_imaginary():
    i = sqrt(-1)
    assert latex(i) == r'i'


def test_builtins_without_args():
    assert latex(sin) == r'\sin'
    assert latex(cos) == r'\cos'
    assert latex(tan) == r'\tan'
    assert latex(log) == r'\log'
    assert latex(Ei) == r'\operatorname{Ei}'
    assert latex(zeta) == r'\zeta'


def test_latex_greek_functions():
    # bug because capital greeks that have roman equivalents should not use
    # \Alpha, \Beta, \Eta, etc.
    s = Function('Alpha')
    assert latex(s) == r'A'
    assert latex(s(x)) == r'A{\left (x \right )}'
    s = Function('Beta')
    assert latex(s) == r'B'
    s = Function('Eta')
    assert latex(s) == r'H'
    assert latex(s(x)) == r'H{\left (x \right )}'

    # bug because diofant.core.numbers.Pi is special
    p = Function('Pi')
    # assert latex(p(x)) == r'\Pi{\left (x \right )}'
    assert latex(p) == r'\Pi'

    # bug because not all greeks are included
    c = Function('chi')
    assert latex(c(x)) == r'\chi{\left (x \right )}'
    assert latex(c) == r'\chi'


def test_translate():
    s = 'Alpha'
    assert translate(s) == 'A'
    s = 'Beta'
    assert translate(s) == 'B'
    s = 'Eta'
    assert translate(s) == 'H'
    s = 'omicron'
    assert translate(s) == 'o'
    s = 'Pi'
    assert translate(s) == r'\Pi'
    s = 'pi'
    assert translate(s) == r'\pi'
    s = 'LamdaHatDOT'
    assert translate(s) == r'\dot{\hat{\Lambda}}'


def test_other_symbols():
    for s in other_symbols:
        assert latex(symbols(s)) == "\\"+s


def test_modifiers():
    # Test each modifier individually in the simplest case (with funny capitalizations)
    assert latex(symbols("xMathring")) == r"\mathring{x}"
    assert latex(symbols("xCheck")) == r"\check{x}"
    assert latex(symbols("xBreve")) == r"\breve{x}"
    assert latex(symbols("xAcute")) == r"\acute{x}"
    assert latex(symbols("xGrave")) == r"\grave{x}"
    assert latex(symbols("xTilde")) == r"\tilde{x}"
    assert latex(symbols("xPrime")) == r"{x}'"
    assert latex(symbols("xddDDot")) == r"\ddddot{x}"
    assert latex(symbols("xDdDot")) == r"\dddot{x}"
    assert latex(symbols("xDDot")) == r"\ddot{x}"
    assert latex(symbols("xBold")) == r"\boldsymbol{x}"
    assert latex(symbols("xnOrM")) == r"\left\|{x}\right\|"
    assert latex(symbols("xAVG")) == r"\left\langle{x}\right\rangle"
    assert latex(symbols("xHat")) == r"\hat{x}"
    assert latex(symbols("xDot")) == r"\dot{x}"
    assert latex(symbols("xBar")) == r"\bar{x}"
    assert latex(symbols("xVec")) == r"\vec{x}"
    assert latex(symbols("xAbs")) == r"\left|{x}\right|"
    assert latex(symbols("xMag")) == r"\left|{x}\right|"
    assert latex(symbols("xPrM")) == r"{x}'"
    assert latex(symbols("xBM")) == r"\boldsymbol{x}"
    # Test strings that are *only* the names of modifiers
    assert latex(symbols("Mathring")) == r"Mathring"
    assert latex(symbols("Check")) == r"Check"
    assert latex(symbols("Breve")) == r"Breve"
    assert latex(symbols("Acute")) == r"Acute"
    assert latex(symbols("Grave")) == r"Grave"
    assert latex(symbols("Tilde")) == r"Tilde"
    assert latex(symbols("Prime")) == r"Prime"
    assert latex(symbols("DDot")) == r"\dot{D}"
    assert latex(symbols("Bold")) == r"Bold"
    assert latex(symbols("NORm")) == r"NORm"
    assert latex(symbols("AVG")) == r"AVG"
    assert latex(symbols("Hat")) == r"Hat"
    assert latex(symbols("Dot")) == r"Dot"
    assert latex(symbols("Bar")) == r"Bar"
    assert latex(symbols("Vec")) == r"Vec"
    assert latex(symbols("Abs")) == r"Abs"
    assert latex(symbols("Mag")) == r"Mag"
    assert latex(symbols("PrM")) == r"PrM"
    assert latex(symbols("BM")) == r"BM"
    assert latex(symbols("hbar")) == r"\hbar"
    # Check a few combinations
    assert latex(symbols("xvecdot")) == r"\dot{\vec{x}}"
    assert latex(symbols("xDotVec")) == r"\vec{\dot{x}}"
    assert latex(symbols("xHATNorm")) == r"\left\|{\hat{x}}\right\|"
    # Check a couple big, ugly combinations
    assert latex(symbols('xMathringBm_yCheckPRM__zbreveAbs')) == r"\boldsymbol{\mathring{x}}^{\left|{\breve{z}}\right|}_{{\check{y}}'}"
    assert latex(symbols('alphadothat_nVECDOT__tTildePrime')) == r"\hat{\dot{\alpha}}^{{\tilde{t}}'}_{\dot{\vec{n}}}"


def test_greek_symbols():
    assert latex(Symbol('alpha')) == r'\alpha'
    assert latex(Symbol('beta')) == r'\beta'
    assert latex(Symbol('gamma')) == r'\gamma'
    assert latex(Symbol('delta')) == r'\delta'
    assert latex(Symbol('epsilon')) == r'\epsilon'
    assert latex(Symbol('zeta')) == r'\zeta'
    assert latex(Symbol('eta')) == r'\eta'
    assert latex(Symbol('theta')) == r'\theta'
    assert latex(Symbol('iota')) == r'\iota'
    assert latex(Symbol('kappa')) == r'\kappa'
    assert latex(Symbol('lambda')) == r'\lambda'
    assert latex(Symbol('mu')) == r'\mu'
    assert latex(Symbol('nu')) == r'\nu'
    assert latex(Symbol('xi')) == r'\xi'
    assert latex(Symbol('omicron')) == r'o'
    assert latex(Symbol('pi')) == r'\pi'
    assert latex(Symbol('rho')) == r'\rho'
    assert latex(Symbol('sigma')) == r'\sigma'
    assert latex(Symbol('tau')) == r'\tau'
    assert latex(Symbol('upsilon')) == r'\upsilon'
    assert latex(Symbol('phi')) == r'\phi'
    assert latex(Symbol('chi')) == r'\chi'
    assert latex(Symbol('psi')) == r'\psi'
    assert latex(Symbol('omega')) == r'\omega'

    assert latex(Symbol('Alpha')) == r'A'
    assert latex(Symbol('Beta')) == r'B'
    assert latex(Symbol('Gamma')) == r'\Gamma'
    assert latex(Symbol('Delta')) == r'\Delta'
    assert latex(Symbol('Epsilon')) == r'E'
    assert latex(Symbol('Zeta')) == r'Z'
    assert latex(Symbol('Eta')) == r'H'
    assert latex(Symbol('Theta')) == r'\Theta'
    assert latex(Symbol('Iota')) == r'I'
    assert latex(Symbol('Kappa')) == r'K'
    assert latex(Symbol('Lambda')) == r'\Lambda'
    assert latex(Symbol('Mu')) == r'M'
    assert latex(Symbol('Nu')) == r'N'
    assert latex(Symbol('Xi')) == r'\Xi'
    assert latex(Symbol('Omicron')) == r'O'
    assert latex(Symbol('Pi')) == r'\Pi'
    assert latex(Symbol('Rho')) == r'P'
    assert latex(Symbol('Sigma')) == r'\Sigma'
    assert latex(Symbol('Tau')) == r'T'
    assert latex(Symbol('Upsilon')) == r'\Upsilon'
    assert latex(Symbol('Phi')) == r'\Phi'
    assert latex(Symbol('Chi')) == r'X'
    assert latex(Symbol('Psi')) == r'\Psi'
    assert latex(Symbol('Omega')) == r'\Omega'

    assert latex(Symbol('varepsilon')) == r'\varepsilon'
    assert latex(Symbol('varkappa')) == r'\varkappa'
    assert latex(Symbol('varphi')) == r'\varphi'
    assert latex(Symbol('varpi')) == r'\varpi'
    assert latex(Symbol('varrho')) == r'\varrho'
    assert latex(Symbol('varsigma')) == r'\varsigma'
    assert latex(Symbol('vartheta')) == r'\vartheta'


def test_builtin_without_args_mismatched_names():
    assert latex(CosineTransform) == r'\operatorname{CosineTransform}'


def test_builtin_no_args():
    assert latex(Chi) == r'\operatorname{Chi}'
    assert latex(gamma) == r'\Gamma'
    assert latex(KroneckerDelta) == r'\delta'
    assert latex(DiracDelta) == r'\delta'
    assert latex(lowergamma) == r'\gamma'


def test_sympyissue_6853():
    p = Function('Pi')
    assert latex(p(x)) == r"\Pi{\left (x \right )}"


def test_Mul():
    e = Mul(-2, x + 1, evaluate=False)
    assert latex(e) == r'- 2 \left(x + 1\right)'
    e = Mul(2, x + 1, evaluate=False)
    assert latex(e) == r'2 \left(x + 1\right)'
    e = Mul(Rational(1, 2), x + 1, evaluate=False)
    assert latex(e) == r'\frac{1}{2} \left(x + 1\right)'
    e = Mul(y, x + 1, evaluate=False)
    assert latex(e) == r'y \left(x + 1\right)'
    e = Mul(-y, x + 1, evaluate=False)
    assert latex(e) == r'- y \left(x + 1\right)'
    e = Mul(-2, x + 1)
    assert latex(e) == r'- 2 x - 2'
    e = Mul(2, x + 1)
    assert latex(e) == r'2 x + 2'
    e = Mul(1, 1, evaluate=False)
    assert latex(e) == r'1 \cdot 1'

    e = x*(z - 1 + x + 1/(x - y))/(x - y)
    assert latex(e) == r'\frac{x}{x - y} \left(x + z - 1 + \frac{1}{x - y}\right)'


def test_Add():
    e = Add(0, 0, evaluate=False)
    assert latex(e) == '0 + 0'

    assert latex(x + y + 1, order='none') == '1 + x + y'


def test_Pow():
    e = Pow(2, 2, evaluate=False)
    assert latex(e) == r'2^{2}'

    # issue sympy/sympy#11011
    assert (latex(Float('1.453e4500')**x) ==
            r"\left(1.453 \cdot 10^{4500}\right)^{x}")

    assert latex(1/sqrt(x)) == r'\frac{1}{\sqrt{x}}'


def test_sympyissue_7180():
    assert latex(Equivalent(x, y)) == r"x \equiv y"
    assert latex(Not(Equivalent(x, y))) == r"x \not\equiv y"


def test_sympyissue_8409():
    assert latex(Rational(1, 2)**n) == r"\left(\frac{1}{2}\right)^{n}"


def test_sympyissue_8470():
    e = parse_expr("-B*A", evaluate=False)
    assert latex(e) == r"A \left(- B\right)"


def test_sympyissue_7117():
    # See also issue sympy/sympy#5031 (hence the evaluate=False in these).
    e = Eq(x + 1, 2*x)
    q = Mul(2, e, evaluate=False)
    assert latex(q) == r"2 \left(x + 1 = 2 x\right)"
    q = Add(6, e, evaluate=False)
    assert latex(q) == r"6 + \left(x + 1 = 2 x\right)"
    q = Pow(e, 2, evaluate=False)
    assert latex(q) == r"\left(x + 1 = 2 x\right)^{2}"


def test_sympyissue_2934():
    assert latex(Symbol(r'\frac{a_1}{b_1}')) == '\\frac{a_1}{b_1}'


def test_diffgeom():
    assert latex(R2) == r'\mathbb{R}^{2}'
    r2 = Manifold('r', 2)
    assert latex(r2) == r'\mathrm{r}'
    foo = Manifold('Foo', 2)
    assert latex(foo) == r'\mathrm{Foo}'
    R2_origin = Patch('origin', R2)
    assert latex(R2_origin) == r'\mathrm{origin}_{\mathbb{R}^{2}}'
    assert latex(R2_r) == r'\mathrm{rectangular}^{\mathrm{origin}}_{\mathbb{R}^{2}}'
    tp1 = TensorProduct(R2.dx, R2.dy)
    assert latex(tp1) == r'\mathrm{d}x\otimes\mathrm{d}y'
    ch = metric_to_Christoffel_2nd(TensorProduct(R2.dx, R2.dx) +
                                   TensorProduct(R2.dy, R2.dy))
    cvd = CovarDerivativeOp(R2.x*R2.e_x, ch)
    assert latex(cvd) == r'\mathbb{\nabla}_{\boldsymbol{\mathrm{x}} \partial_{x}}'

    g = Function('g')
    s_field = g(R2.x, R2.y)
    dg = Differential(s_field)
    assert latex(dg) == r'd(g{\left (\boldsymbol{\mathrm{x}},' \
                        r'\boldsymbol{\mathrm{y}} \right )})'


def test_sympyissue_10489():
    latexSymbolWithBrace = 'C_{x_{0}}'
    s = Symbol(latexSymbolWithBrace)
    assert latex(s) == latexSymbolWithBrace
    assert latex(cos(s)) == r'\cos{\left (C_{x_{0}} \right )}'


def test_sympyissue_10889():
    A, B = symbols('A B', commutative=False)
    e = Mul(-1, A*B - B*A)
    assert latex(e) == "- (A B - B A)"
