import warnings

from sympy.core.decorators import deprecated
from sympy.utilities.exceptions import SymPyDeprecationWarning


def test_deprecated():
    @deprecated(useinstead="bar", issue=1234, deprecated_since_version="0.7.2")
    def foo():
        return

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        foo()
        assert len(w) == 1
        assert issubclass(w[-1].category, SymPyDeprecationWarning)
        assert str(w[-1].message) == """\n\n\
foo has been deprecated since SymPy 0.7.2. Use bar instead. See\n\
https://github.com/sympy/sympy/issues/1234 for more info.\n\
"""
