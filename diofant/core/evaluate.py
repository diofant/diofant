from contextlib import contextmanager


global_evaluate = [True]


@contextmanager
def evaluate(x):
    """Control automatic evaluation.

    This context managers controls whether or not all Diofant functions evaluate
    by default.

    Note that much of Diofant expects evaluated expressions.  This functionality
    is experimental and is unlikely to function as intended on large
    expressions.

    Examples
    ========

    >>> x + x
    2*x
    >>> with evaluate(False):
    ...     x + x
    x + x

    """
    old = global_evaluate[0]

    global_evaluate[0] = x
    yield
    global_evaluate[0] = old
