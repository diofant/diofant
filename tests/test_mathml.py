import sys
sys.path.append(".")

from sympy import *


x = Symbol('x')
f = integrate(log(x), (x,1,2), evaluate=False)

def test_mathml_1():
    assert f.mathml == "<apply><int/><bvar><ci sympy:assumptions='is_commutative:True;'> x </ci></bvar><lowlimit><cn sympy:assumptions='is_real:True;is_commutative:True;'> 1 </cn></lowlimit><uplimit><cn sympy:assumptions='is_real:True;is_commutative:True;'> 2 </cn></uplimit><apply><log/> <ci sympy:assumptions='is_commutative:True;'> x </ci> </apply></apply>"
    assert (x**2 + x +1/x).mathml in ["<apply><plus/><apply><power/><ci sympy:assumptions='is_commutative:True;'> x </ci><cn sympy:assumptions='is_real:True;is_commutative:True;'> -1 </cn></apply><ci sympy:assumptions='is_commutative:True;'> x </ci><apply><power/><ci sympy:assumptions='is_commutative:True;'> x </ci><cn sympy:assumptions='is_real:True;is_commutative:True;'> 2 </cn></apply></apply>" , 
"<apply><plus/><ci sympy:assumptions='is_commutative:True;'> x </ci><apply><power/><ci sympy:assumptions='is_commutative:True;'> x </ci><cn sympy:assumptions='is_real:True;is_commutative:True;'> 2 </cn></apply><apply><power/><ci sympy:assumptions='is_commutative:True;'> x </ci><cn sympy:assumptions='is_real:True;is_commutative:True;'> -1 </cn></apply></apply>" ]

def test_c2p():
    """This tests some optional routines that depend on libxslt1 (which is optional)"""
    try:
        from sympy.modules.mathml import c2p
             
        #assert c2p(f.mathml) == result
    except ImportError:
        pass
