import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.printing import latex

x = Symbol('x')

def test_latex_basic():
    assert latex(1+x) == "$1+x$"
    assert latex(x**2) == "${x}^{2}$"
    assert latex(x**(1+x)) == "${x}^{(1+x)}$"
    
def test_latex_functions():
    assert latex(exp(x)) == "${e}^{x}$"
    
def test_latex_integrals():
    assert latex(integrate(log(x), x, evaluate=False)) == "$\int log(x)\,dx$"
    assert latex(integrate(x**2, (x,0,1), evaluate=False)) == "$\int^0_1 {x}^{2}\,dx$"
    
def test_latex_limits():
    assert latex(limit(x, x, oo, evaluate=False)) == "$\lim_{x \to \infty}x$"
    