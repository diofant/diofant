import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.printing import pretty

x = Symbol('x')

def test_pretty_basic():
    assert pretty( (x**2) ) == ' 2\nx '
    assert pretty( (x**2 + x + 1)) in ['     2\n1+x+x ']
    assert pretty( oo ) == "oo"

    
def test_pretty_functions():
    assert pretty( (2*x + exp(x)) ) in ['     x\n2*x+e ',' x    \ne +2*x']
    assert pretty( sqrt(2) ) == '  ___\n\\/ 2 '
    assert pretty( sqrt(2+pi) ) == '  ______\n\\/ 2+pi '
    # nesting of square roots
    assert pretty( sqrt((sqrt(x+1))+1) ) == '    ___________\n   /     _____ \n \\/  1+\\/ 1+x  '
    
def test_pretty_integrals():
    f = integrate(log(x), x, evaluate=False)
    assert pretty( f ) in ['/          \n| log(x) dx\n/          ']
    
def test_pretty_limits():
    assert pretty( limit(x, x, oo, evaluate=False) ) == ' lim x\nx->oo '
    assert pretty( limit(x**2, x, 0, evaluate=False) ) == '     2\nlim x \nx->0  '  
