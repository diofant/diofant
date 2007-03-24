import sys
sys.path.append(".")

from sympy import *
from sympy.modules.mathml import *

x = Symbol('x')
f = integrate(log(x), (x,1,2), evaluate=False)

def test_mathml_1():
    assert f.mathml == '<apply><int/><bvar><ci> x </ci></bvar><lowlimit><cn> 1 </cn></lowlimit><uplimit><cn> 2 </cn></uplimit><apply><log/> <ci> x </ci> </apply></apply>'
    assert (x**2 + x +1/x).mathml == '<apply><plus/><apply><power/><ci> x </ci><cn> 2 </cn></apply><apply><power/><ci> x </ci><cn> -1 </cn></apply><ci> x </ci></apply>'
    
