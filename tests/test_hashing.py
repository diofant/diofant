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
    m1.addstr("<class 'sym.symbol.Symbol'>")
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

def test_bug4():
    from sym.core.hashing import mhash
    m1=mhash()
    m1.addstr("<class 'sym.core.symbol.Symbol'>")
    m1.addstr("x")

    m2a=mhash()
    m2a.addstr("<class 'sym.co.functions.log'>")
    m2a.addint(m1.value)

    m2=mhash()
    m2.addstr("<class 'sym.co.functions.log'>")
    m2.addint(m2a.value)
    assert m1.value!=m2.value

    h = mhash()
    h.addstr("<class 'sym.core.symbol.Symbol'>")
    h.addstr("x")
    hp=mhash()
    hp.addstr("<class 'sym.core.functions.log'>")
    hp.addint(h.value)
    h2 = mhash()
    h2.addstr("<class 'sym.core.functions.log'>")
    h2.addint(hp.value)
    assert h.value!=h2.value

def test_bug5():
    x=g.Symbol("x")
    assert g.cos(x)/g.sin(x)!=g.sin(x)/g.cos(x)

def test_bug5h():
    from sym.core import hashing
    mhash=hashing.mhash()
    mhash.addstr("<class 'sym.core.numbers.Rational'>")
    mhash.addint(-1)
    mhash.addint(1)
    mhash=mhash.value

    xhash=hashing.mhash()
    xhash.addstr("<class 'sym.core.symbol.Symbol'>")
    xhash.addstr("x")
    xhash=xhash.value

    a1hash=hashing.mhash()
    a1hash.addstr("<class 'sym.modules.trigonometric.cos'>")
    a1hash.addint(xhash)
    a1hash=a1hash.value

    b1hash=hashing.mhash()
    b1hash.addstr("<class 'sym.modules.trigonometric.sin'>")
    b1hash.addint(xhash)
    b1hash=b1hash.value

    b2hash=hashing.mhash()
    b2hash.addstr("<class 'sym.core.power.Pow'>")
    b2hash.add(a1hash)
    b2hash.add(mhash)
    b2hash=b2hash.value

    a2hash=hashing.mhash()
    a2hash.addstr("<class 'sym.core.power.Pow'>")
    a2hash.add(b1hash)
    a2hash.add(mhash)
    a2hash=a2hash.value

    m=hashing.mhash()
    m.addstr("<class 'sym.core.addmul.Mul'>")
    m.add(a1hash)
    m.add(a2hash)

    m2=hashing.mhash()
    m2.addstr("<class 'sym.core.addmul.Mul'>")
    m2.add(b1hash)
    m2.add(b2hash)

    assert m.value!=m2.value

def test_bug5h2():
    from sym.core import hashing
    mhash=hashing.mhash()
    mhash.addstr("<class 'sym.core.numbers.Rational'>")
    mhash.addint(-1)
    mhash.addint(1)
    mhash=mhash.value

    xhash=hashing.mhash()
    xhash.addstr("<class 'sym.core.symbol.Symbol'>")
    xhash.addstr("x")
    xhash=xhash.value

    a1hash=hashing.mhash()
    a1hash.addstr("<class 'sym.modules.trigonometric.cos'>")
    a1hash.addint(xhash)
    a1hash=a1hash.value

    b1hash=hashing.mhash()
    b1hash.addstr("<class 'sym.modules.trigonometric.sin'>")
    b1hash.addint(xhash)
    b1hash=b1hash.value

    b2hash=hashing.mhash()
    b2hash.addstr("<class 'sym.core.power.Pow'>")
    b2hash.addint(a1hash)
    b2hash.addint(mhash)
    b2hash=b2hash.value

    a2hash=hashing.mhash()
    a2hash.addstr("<class 'sym.core.power.Pow'>")
    a2hash.addint(b1hash)
    a2hash.addint(mhash)
    a2hash=a2hash.value

    m=hashing.mhash()
    m.addstr("<class 'sym.core.addmul.Mul'>")
    m.addint(a1hash)
    m.addint(a2hash)

    m2=hashing.mhash()
    m2.addstr("<class 'sym.core.addmul.Mul'>")
    m2.addint(b1hash)
    m2.addint(b2hash)

    assert m.value!=m2.value
