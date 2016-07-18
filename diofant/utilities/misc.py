"""Miscellaneous stuff that doesn't really fit anywhere else."""

import sys
from textwrap import fill, dedent


def filldedent(s, w=70):
    """
    Strips leading and trailing empty lines from a copy of `s`, then dedents,
    fills and returns it.

    Empty line stripping serves to deal with docstrings like this one that
    start with a newline after the initial triple quote, inserting an empty
    line at the beginning of the string."""
    return '\n' + fill(dedent(str(s)).strip('\n'), width=w)


size = sys.maxsize
if size > 2**32:
    ARCH = "64-bit"
else:  # pragma: no cover
    ARCH = "32-bit"  # XXX we don't test this


# XXX: PyPy doesn't support hash randomization
HASH_RANDOMIZATION = getattr(sys.flags, 'hash_randomization', False)


def debug(*args):
    """
    Print ``*args`` if DIOFANT_DEBUG is True, else do nothing.
    """
    from diofant import DIOFANT_DEBUG
    if DIOFANT_DEBUG:
        print(*args, file=sys.stderr)
