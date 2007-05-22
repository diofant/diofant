import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.printing.pretty import pretty

x = Symbol('x')
y = Symbol('y')

def test_pretty_basic():
    assert pretty( -Rational(1)/2 ) == '  1\n- -\n  2'
    assert pretty( (x**2) ) == ' 2\nx '
    assert pretty( (x**2 + x + 1))  == '         2\n1 + x + x '
    assert pretty( 1-x ) == '1 - x'
    assert pretty( 1-2*x ) == '1 - 2*x'
    assert pretty( 1-Rational(3,2)*y/x ) == '    3*y\n1 - ---\n    2*x'
    assert pretty( 1/x ) == '1\n-\nx'
    assert pretty( x/y ) == 'x\n-\ny'
    assert pretty( -x/y ) == '-x\n--\ny '
    assert pretty( (x-2)/y ) == '-2 + x\n------\n  y   '
    assert pretty( oo ) == "oo"

    
def test_pretty_functions():
    assert pretty( (2*x + exp(x)) ) in [' x      \ne  + 2*x', '       x\n2*x + e ']
    assert pretty( sqrt(2) ) == '  ___\n\\/ 2 '
    assert pretty( sqrt(2+pi) ) == '  ________\n\\/ 2 + pi '
    # nesting of square roots
    assert pretty( sqrt((sqrt(x+1))+1) ) == '    _______________\n   /       _______ \n \\/  1 + \\/ 1 + x  '
    assert pretty( diff(log(x), x, evaluate=False) ) == 'd       \n--log(x)\ndx      '
    assert pretty( diff(log(x), x, evaluate=False) + x ) == '    /d       \\\nx + |--log(x)|\n    \\dx      /'
    
    
def test_pretty_integrals():
    f_1 = integrate(log(x), x, evaluate=False)
    assert pretty( f_1 ) == '/          \n| log(x) dx\n/          '
    
    f_2 = integrate(x**2, x, evaluate=False)
    assert pretty( f_2 ) == '/      \n|  2   \n| x  dx\n/      '
    
    f_3 = integrate(x**(2**x), x, evaluate=False) # double nesting of pow
    assert pretty( f_3 ) == '/         \n|  / x\\   \n|  \\2 /   \n| x     dx\n/         '
    
def test_pretty_limits():
    assert pretty( limit(x, x, oo, evaluate=False) ) == ' lim x\nx->oo '
    assert pretty( limit(x**2, x, 0, evaluate=False) ) == '     2\nlim x \nx->0  '  
