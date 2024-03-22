"""Useful utility decorators."""


def conserve_mpmath_dps(func):
    """After the function finishes, resets the value of mpmath.mp.dps to
    the value it had before the function was run.
    """
    import functools

    import mpmath

    def func_wrapper():
        dps = mpmath.mp.dps
        try:
            func()
        finally:
            mpmath.mp.dps = dps

    func_wrapper = functools.update_wrapper(func_wrapper, func)
    return func_wrapper


class no_attrs_in_subclass:
    """Don't 'inherit' certain attributes from a base class

    >>> class A:
    ...     x = 'test'

    >>> A.x = no_attrs_in_subclass(A, A.x)

    >>> class B(A):
    ...     pass

    >>> hasattr(A, 'x')
    True
    >>> hasattr(B, 'x')
    False
    """

    def __init__(self, cls, f):
        """Initialize self."""
        self.cls = cls
        self.f = f

    def __get__(self, instance, owner=None):
        if owner == self.cls:
            return self.f
        raise AttributeError


def doctest_depends_on(exe=None, modules=None, disable_viewers=None):
    """Adds metadata about the dependencies which need to be met for doctesting
    the docstrings of the decorated objects.
    """

    def depends_on_deco(fn):
        fn._doctest_depends_on = {'exe': exe, 'modules': modules,
                                  'disable_viewers': disable_viewers}
        return fn
    return depends_on_deco
