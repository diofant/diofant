import sys
sys.path.append(".")

import sym as g

def test_hashing_class():
    h1=g.core.hashing.mhash()
    h1.addint(-1)
    h1.addint(-2)

    h2=g.core.hashing.mhash()
    h2.addint(-1)
    h2.addint(-2)
    assert h1.value == h2.value

    h3=g.core.hashing.mhash()
    h3.addint(-2)
    h3.addint(-1)
    assert h1.value != h3.value

def test_basic_class():
    n0=g.Rational(-0)
    n1=g.Rational(-1)
    n2=g.Rational(-2)
    n3=g.Rational(-3)
    assert not n1.hash()==n2.hash()
    assert not n1.hash()==n3.hash()
    assert not n0.hash()==n1.hash()
    x=g.Symbol("x")
    y=g.Symbol("y")
    z=g.Symbol("y1")
    z2=g.Symbol("y1")
    assert not x.hash()==y.hash()
    assert not x.hash()==z.hash()
    assert not y.hash()==z.hash()
    assert z.hash()==z2.hash()

def test_bug():
    from sym.core.hashing import mhash
    m1=mhash()
    m1.addstr("<class 'sym.Symbol.Symbol'>")
    m1.addstr("y")

    m2a=mhash()
    m2a.addstr("<class 'sym.functions.log'>")
    m2a.addint(m1.value)

    m2=mhash()
    m2.addstr("<class 'sym.functions.log'>")
    m2.addint(m2a.value)
    assert m1.value!=m2.value

def test_bug2():
    from sym.core.hashing import mhash
    num=123456

    m2a=mhash()
    m2a.addint(num)

    m2=mhash()
    m2.addint(m2a.value)
    assert num!=m2.value

def test_bug3():
    x=g.Symbol('x')
    y=x*x
    assert y.subs(x,g.Rational(3))==9
    assert y.subs(x,g.Real(3.2))!=9
