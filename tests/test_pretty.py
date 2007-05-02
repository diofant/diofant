import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.printing import pretty

x = Symbol('x')
y = Symbol('y')

def test_pretty_basic():
    assert pretty( (x**2) ) == ' 2\nx '
    assert pretty( 1-x ) == '1-x'
    assert pretty( 1-2*x ) == '1-2*x'
    assert pretty( (x**2 + x + 1))  == '     2\n1+x+x '
    assert pretty( 1/x ) == '1\n-\nx'
    assert pretty( x/y ) == "x\n-\ny"
    assert pretty( (x-2)/y ) == '-2+x\n----\n y  '
    assert pretty( -x/y ) == '-x\n--\ny '
    assert pretty( oo ) == "oo"

    
def test_pretty_functions():
    assert pretty( (2*x + exp(x)) ) in ['     x\n2*x+e ',' x    \ne +2*x']
    assert pretty( sqrt(2) ) == '  ___\n\\/ 2 '
    assert pretty( sqrt(2+pi) ) == '  ______\n\\/ 2+pi '
    # nesting of square roots
    assert pretty( sqrt((sqrt(x+1))+1) ) == '    ___________\n   /     _____ \n \\/  1+\\/ 1+x  '
    assert pretty( diff(log(x), x, evaluate=False) ) == 'd       \n--log(x)\ndx      '
    assert pretty( diff(log(x), x, evaluate=False) + x ) == '  /d       \\\nx+|--log(x)|\n  \\dx      /'
    
    
def test_pretty_integrals():
    f = integrate(log(x), x, evaluate=False)
    assert pretty( f ) in ['/          \n| log(x) dx\n/          ']
    
def test_pretty_limits():
    assert pretty( limit(x, x, oo, evaluate=False) ) == ' lim x\nx->oo '
    assert pretty( limit(x**2, x, 0, evaluate=False) ) == '     2\nlim x \nx->0  '  
