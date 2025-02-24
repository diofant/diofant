"""Miscellaneous stuff that doesn't really fit anywhere else."""

import textwrap


def as_int(n):
    """
    Convert the argument to a builtin integer.

    The return value is guaranteed to be equal to the input. ValueError is
    raised if the input has a non-integral value.

    Examples
    ========

    >>> 3.0
    3.0
    >>> as_int(3.0)  # convert to int and test for equality
    3
    >>> int(sqrt(10))
    3
    >>> as_int(sqrt(10))
    Traceback (most recent call last):
    ...
    ValueError: ... is not an integer

    """
    try:
        result = int(n)
        if result != n:
            raise TypeError
    except TypeError as exc:
        raise ValueError(f'{n} is not an integer') from exc
    return result


def filldedent(s, w=70):
    """
    Strips leading and trailing empty lines from a copy of `s`, then dedents,
    fills and returns it.

    Empty line stripping serves to deal with docstrings like this one that
    start with a newline after the initial triple quote, inserting an empty
    line at the beginning of the string.
    """
    return '\n' + textwrap.fill(textwrap.dedent(str(s)).strip('\n'), width=w)
