import sys
sys.path.append(".")

import sym as g

def test_hashing_class():
    h1=g.hashing.mhash()
    h1.addint(-1)
    h1.addint(-2)

    h2=g.hashing.mhash()
    h2.addint(-1)
    h2.addint(-2)
    assert h1.value == h2.value

    h3=g.hashing.mhash()
    h3.addint(-2)
    h3.addint(-1)
    assert h1.value != h3.value

def test_basic_class():
    n0=g.rational(-0)
    n1=g.rational(-1)
    n2=g.rational(-2)
    n3=g.rational(-3)
    assert not n1.hash()==n2.hash()
    assert not n1.hash()==n3.hash()
    assert not n0.hash()==n1.hash()
    x=g.symbol("x")
    y=g.symbol("y")
    z=g.symbol("y1")
    z2=g.symbol("y1")
    assert not x.hash()==y.hash()
    assert not x.hash()==z.hash()
    assert not y.hash()==z.hash()
    assert z.hash()==z2.hash()
