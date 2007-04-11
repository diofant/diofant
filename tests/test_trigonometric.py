import sys
sys.path.append(".")

from sympy import *

def test_cos():
    cos(pi) == -1
    float(cos(1)) == 0.54030230586813977
    cos(1).evalf() == 0.540302305868171952183697906453

def test_sin():
    sin(pi) == 0
    float(sin(1)) == 0.8414709848078965
