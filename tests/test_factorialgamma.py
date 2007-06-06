import sys
sys.path.append(".")

import py

from sympy import *
from sympy.modules.factorialgamma import *
from sympy.modules.printing.latex import latex

x = Symbol('x')

def test_factorial():
    assert [factorial(t) for t in [0,1,2,3,4]] == [1,1,2,6,24]
    assert (factorial(-1) == oo) == True
    assert factorial(Rational(1,2)) == Rational(1,2)*sqrt(pi)
    assert factorial(Rational(3,2)) == Rational(3,4)*sqrt(pi)
    assert factorial(Rational(5,2)) == Rational(15,8)*sqrt(pi)
    assert factorial(Rational(-1,2)) == sqrt(pi)
    assert factorial(Rational(-3,2)) == -2*sqrt(pi)
    assert factorial(Rational(-5,2)) == Rational(4,3)*sqrt(pi)
    assert factorial(Rational(-17,2)) == Rational(256,2027025)*sqrt(pi)
    assert latex(factorial(x, evaluate=False)) == "$x!$"
    assert latex(factorial(-4, evaluate=False)) == "$(-4)!$"
    assert latex(factorial(-x, evaluate=False)) == "$(- x)!$"

def test_factorial2():
    assert factorial2(0) == 1
    assert factorial2(2) == 2
    assert factorial2(4) == 2*4
    assert factorial2(6) == 2*4*6
    assert factorial2(1) == 1
    assert factorial2(3) == 3
    assert factorial2(5) == 3*5
    assert factorial2(7) == 3*5*7
    assert (factorial2(-2) == oo) == True
    assert (factorial2(-4) == oo) == True
    assert factorial2(-1) == 1
    assert factorial2(-3) == -1
    assert factorial2(-7) == Rational(-1,15)
    assert factorial2(-9) == Rational(1,105)
    assert latex(factorial2(x, evaluate=False)) == "$x!!$"
    assert latex(factorial2(-4, evaluate=False)) == "$(-4)!!$"
    assert latex(factorial2(-x, evaluate=False)) == "$(- x)!!$"

def test_factorial_quotient():
    assert factorial_quotient(x+5, x+5) == 1
    assert factorial_quotient(x+1, x) == 1+x
    assert factorial_quotient(x+2, x) == (1+x)*(2+x)
    assert factorial_quotient(x+3, x) == (1+x)*(2+x)*(3+x)
    assert factorial_quotient(x-1, x) == (1/x)
    assert factorial_quotient(x-2, x) == 1/(x*(-1+x))
    assert factorial_quotient(x-3, x) == 1/(x*(-1+x)*(-2+x))

def test_gamma():
    assert gamma(0) == oo
    assert gamma(1) == 1
    assert gamma(2) == 1
    assert gamma(3) == 2
    assert gamma(Rational(1,2)) == sqrt(pi)
    assert latex(gamma(3+x)) == "$\Gamma(3+x)$"
