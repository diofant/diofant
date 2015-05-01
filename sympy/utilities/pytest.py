"""py.test hacks"""

from __future__ import print_function, division

import sys
import functools
import os

from sympy.core.compatibility import get_function_name

import py
from py.test import skip, raises

ON_TRAVIS = os.getenv('TRAVIS_BUILD_NUMBER', None)

XFAIL = py.test.mark.xfail
slow = py.test.mark.slow

def SKIP(reason):
    def skipping(func):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            skip(reason)
        return inner

    return skipping
