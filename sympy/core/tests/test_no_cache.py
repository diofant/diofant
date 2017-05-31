import imp
import os
import sys

os.environ['USE_CACHE'] = 'no'
imp.reload(sys.modules["sympy.core.cache"])


from sympy.core.cache import print_cache, CACHE
from sympy import Symbol


def test_print_cache():
    assert CACHE == []
    x = Symbol('x')
    assert CACHE == []
