"""Caching facility for Diofant."""

import functools
import os

from .evaluate import global_evaluate


# global cache registry: [] of (item, {})
CACHE: list[object] = []


def print_cache():
    """Print cache content."""
    for item in CACHE:
        print(item.__qualname__, item.cache_info())


def clear_cache():
    """Clear cache content."""
    for item in CACHE:
        item.cache_clear()


USE_CACHE = os.getenv('DIOFANT_USE_CACHE', 'True') == 'True'


def cacheit(f, maxsize=None):
    """Caching decorator.

    The result of cached function must be *immutable*.

    Examples
    ========

    >>> @cacheit
    ... def f(a, b):
    ...     print(a, b)
    ...     return a + b
    >>> f(x, y)
    x y
    x + y
    >>> f(x, y)
    x + y

    """
    if USE_CACHE:
        cfunc = functools.lru_cache(maxsize=maxsize, typed=True)(f)

        def wrapper(*args, **kwargs):
            try:
                if global_evaluate[0] and kwargs.get('evaluate', True):
                    return cfunc(*args, **kwargs)
            except TypeError:
                pass
            return f(*args, **kwargs)

        wrapper.cache_info = cfunc.cache_info
        wrapper.cache_clear = cfunc.cache_clear
        functools.update_wrapper(wrapper, f)

        CACHE.append(wrapper)
        return wrapper
    return f
