import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.printing.latex import latex

x = Symbol('x')

def test_latex_basic():
    assert latex(1+x) == "$1+x$"
    assert latex(x**2) == "${x}^{2}$"
    assert latex(x**(1+x)) == "${x}^{(1+x)}$"

def test_latex_symbols():
    Gamma, lmbda, rho = Symbol('Gamma'), Symbol('lambda'), Symbol('rho')
    mass, volume = Symbol('mass'), Symbol('volume')
    assert latex(Gamma + lmbda) == "$\Gamma+\lambda$"
    assert latex(Gamma * lmbda) == "$\Gamma \lambda$"
    assert latex(volume * rho == mass) == r"$\mathrm{volume} \cdot \rho = \mathrm{mass}$"
    assert latex(volume / mass * rho == 1) == r"$\mathrm{volume} \cdot {\mathrm{mass}}^{(-1)} \cdot \rho = 1$"
    assert latex(mass**3 * volume**3) == r"${\mathrm{mass}}^{3} \cdot {\mathrm{volume}}^{3}$"

def test_latex_functions():
    assert latex(exp(x)) == "${e}^{x}$"
    assert latex(exp(I*pi)) == "${e}^{\mathrm{i} \pi}$"

def test_latex_integrals():
    assert latex(integrate(log(x), x, evaluate=False)) == "$\int log(x)\,dx$"
    assert latex(integrate(x**2, (x,0,1), evaluate=False)) == "$\int^0_1 {x}^{2}\,dx$"
    
def test_latex_limits():
    assert latex(limit(x, x, oo, evaluate=False)) == "$\lim_{x \to \infty}x$"
