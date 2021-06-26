"""Miscellaneous stuff that doesn't really fit anywhere else."""

import textwrap


def filldedent(s, w=70):
    """
    Strips leading and trailing empty lines from a copy of `s`, then dedents,
    fills and returns it.

    Empty line stripping serves to deal with docstrings like this one that
    start with a newline after the initial triple quote, inserting an empty
    line at the beginning of the string.
    """
    return '\n' + textwrap.fill(textwrap.dedent(str(s)).strip('\n'), width=w)
