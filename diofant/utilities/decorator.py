"""Useful utility decorators. """

import inspect
import sys
import types


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
        self.cls = cls
        self.f = f

    def __get__(self, instance, owner=None):
        if owner == self.cls:
            if hasattr(self.f, '__get__'):
                return self.f.__get__(instance, owner)
            return self.f
        raise AttributeError


def doctest_depends_on(exe=None, modules=None, disable_viewers=None):
    """Adds metadata about the dependencies which need to be met for doctesting
    the docstrings of the decorated objects.
    """

    def depends_on_deco(fn):
        fn._doctest_depends_on = {'exe': exe, 'modules': modules,
                                  'disable_viewers': disable_viewers}

        # once we drop py2.5 support and use class decorators this evaluates
        # to True
        if inspect.isclass(fn):
            fn._doctest_depdends_on = no_attrs_in_subclass(fn, fn._doctest_depends_on)
        return fn
    return depends_on_deco


def public(obj):
    """
    Append ``obj``'s name to global ``__all__`` variable (call site).

    By using this decorator on functions or classes you achieve the same goal
    as by filling ``__all__`` variables manually, you just don't have to repeat
    yourself (object's name). You also know if object is public at definition
    site, not at some random location (where ``__all__`` was set).

    Note that in multiple decorator setup (in almost all cases) ``@public``
    decorator must be applied before any other decorators, because it relies
    on the pointer to object's global namespace. If you apply other decorators
    first, ``@public`` may end up modifying the wrong namespace.

    Examples
    ========

    >>> __all__
    Traceback (most recent call last):
    ...
    NameError: name '__all__' is not defined

    >>> @public
    ... def some_function():
    ...     pass

    >>> __all__
    ['some_function']

    """
    if isinstance(obj, types.FunctionType):
        ns = obj.__globals__
        name = obj.__name__
    elif isinstance(obj, type):
        ns = sys.modules[obj.__module__].__dict__
        name = obj.__name__
    else:
        raise TypeError("expected a function or a class, got %s" % obj)

    if "__all__" not in ns:
        ns["__all__"] = [name]
    else:
        ns["__all__"].append(name)

    return obj
