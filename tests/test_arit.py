import sys
sys.path.append(".")

import sympy as g

def testSymbol():
    a=g.Symbol("a")
    b=g.Symbol("b")
    c=g.Symbol("c")
    assert str(a)=="a"
    assert str(b)=="b"
    e=a*b
    assert e==a*b
    assert a*b*b==a*b**2
    assert a*b*b+c==c+a*b**2
    assert a*b*b-c==-c+a*b**2

def testarit():
    a=g.Symbol("a")
    b=g.Symbol("b")
    c=g.Symbol("c")
    p=g.Rational(5)
    e=a*b
    assert e == a*b
    e=a*b+b*a
    assert e == 2*a*b
    e=a*b+b*a+a*b+p*b*a
    assert e == 8*a*b
    e=a*b+b*a+a*b+p*b*a+a
    assert e == a+8*a*b
    e=a+a
    assert e == 2*a
    e=a+b+a
    assert e == b+2*a
    e=a+b*b+a+b*b
    assert e == 2*a+2*b**2
    e=a+g.Rational(2)+b*b+a+b*b+p
    assert e == 7+2*a+2*b**2
    e=(a+b*b+a+b*b)*p
    assert e == 5*(2*a+2*b**2)
    e=(a*b*c+c*b*a+b*a*c)*p
    assert e == 15*a*b*c
    e=(a*b*c+c*b*a+b*a*c)*p-g.Rational(15)*a*b*c
    assert e == g.Rational(0)
    e=g.Rational(50)*(a-a)
    assert e == g.Rational(0)
    e=b*a-b-a*b+b
    assert e == g.Rational(0)
    e=a*b+c**p
    assert e == a*b+c**5
    e=a/b
    assert e == a*b**(-1)
    e=a*2*2
    assert e == 4*a
    e=2+a*2/2
    assert e == 2+a
    e=2-a-2
    assert e == -a
    e=2*a*2
    assert e == 4*a
    e=2/a/2
    assert e == a**(-1)
    e=2**a**2
    assert e == 2**(a**2)

def testdiv():
    a=g.Symbol("a")
    b=g.Symbol("b")
    c=g.Symbol("c")
    e=a/b
    assert e == a*b**(-1)
    e=a/b+c/2
    assert e == a*b**(-1)+g.Rational(1)/2*c
    e=(1-b)/(b-1)
    assert e == (1+-b)*((-1)+b)**(-1)

def testpow():
    a=g.Symbol("a")
    b=g.Symbol("b")
    c=g.Symbol("c")
    n1=g.Rational(1)
    n2=g.Rational(2)
    n5=g.Rational(5)
    e=a*a
    assert e == a**2
    e=a*a*a
    assert e == a**3
    e=a*a*a*a**g.Rational(6)
    assert e == a**9
    e=a*a*a*a**g.Rational(6)-a**g.Rational(9)
    assert e == g.Rational(0)
    e=a**(b+c)*a**(-b)
    assert e == a**c
    e=a**(b+c)*a*a**(-b)*a**(-c)/a
    assert e == g.Rational(1)
    e=a**(b-b)
    assert e == g.Rational(1)
    e=(a-a)**b
    assert e == g.Rational(0)
    e=(a+g.Rational(1)-a)**b
    assert e == g.Rational(1)

    e=(a+b+c)**n2
    assert e == (a+b+c)**2
    assert e.expand() == 2*b*c+2*a*c+2*a*b+a**2+c**2+b**2

    e=(a+b)**n2
    assert e == (a+b)**2
    assert e.expand() == 2*a*b+a**2+b**2

    e=(a+b)**(n1/n2)
    assert e == (a+b)**(g.Rational(1)/2)
    assert e.expand() == (a+b)**(g.Rational(1)/2)

    n=n5**(n1/n2)
    assert n == g.Rational(5)**(g.Rational(1)/2)
    e=n*a*b-n*b*a
    assert e == g.Rational(0)
    e=n*a*b+n*b*a
    assert e == 2*a*b*5**(g.Rational(1)/2)
    assert e.diff(a) == 2*b*5**(g.Rational(1)/2)
    assert e.diff(a) == 2*b*5**(g.Rational(1)/2)
    e=a/b**2
    assert e == a*b**(-2)

def testexpand():
    a=g.Symbol("a")
    b=g.Symbol("b")
    c=g.Symbol("c")
    p=g.Rational(5)
    e=(a+b)*c
    assert e == c*(a+b)
    assert (e.expand()-a*c-b*c) == g.Rational(0)
    e=(a+b)*(a+b)
    assert e == (a+b)**2
    assert e.expand() == 2*a*b+a**2+b**2
    e=(a+b)*(a+b)**g.Rational(2)
    assert e == (a+b)**3
    assert e.expand() == 3*b*a**2+3*a*b**2+a**3+b**3
    assert e.expand() == 3*b*a**2+3*a*b**2+a**3+b**3
    e=(a+b)*(a+c)*(b+c)
    assert e == (a+c)*(a+b)*(b+c)
    assert e.expand() == 2*a*b*c+b*a**2+c*a**2+b*c**2+a*c**2+c*b**2+a*b**2
    e=(a+g.Rational(1))**p
    assert e == (1+a)**5
    assert e.expand() == 1+5*a+10*a**2+10*a**3+5*a**4+a**5
    e=(a+b+c)*(a+c+p)
    assert e == (5+a+c)*(a+b+c)
    assert e.expand() == 5*a+5*b+5*c+2*a*c+b*c+a*b+a**2+c**2
    x=g.Symbol("x")
    s=g.exp(x*x)-1
    e=s.series(x,4)/x**2
    assert e == (x**2+g.Rational(1)/2*x**4)*x**(-2)
    assert e.expand() ==  1+g.Rational(1)/2*x**2

def testncmul():
    A=g.NCSymbol("A")
    B=g.NCSymbol("B")
    C=g.NCSymbol("C")
    b=g.Symbol("b")
    assert A*B != B*A
    assert A*B*C != C*B*A
    assert A*b*B*3*C == 3*b*A*B*C
    assert A*b*B*3*C != 3*b*B*A*C
    assert A*b*B*3*C == 3*A*B*C*b

    assert A+B==B+A
    assert (A+B)*C != C*(A+B)

    assert C*(A+B)*C != C*C*(A+B)

    assert (C*(A+B)).expand() == C*A+C*B
    assert (C*(A+B)).expand() != A*C+B*C

    assert A*A == A**2
    assert (A+B)*(A+B) == (A+B)**2
    assert ((A+B)**2).expand() == A**2 + A*B + B*A +B**2

    assert A**-1  * A == 1
