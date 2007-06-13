import sys
sys.path.append(".")

import py

from sympy.modules.specfun import *

def test_import():
    assert factorial(3) == 6