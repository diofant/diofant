""" Caching facility for Diofant """

import os

from cachetools import cached

from .evaluate import global_evaluate


# global cache registry: [] of (item, {})
CACHE = []


def print_cache():
    """Print cache content"""

    for item, cache in CACHE:
        item = str(item)

        if cache:
            head = '='*len(item)
            print(head)
            print(item)
            print(head)

        for k, v in list(cache.items()):
            print('  %s : %s' % (k, v))


def clear_cache():
    """Clear cache content"""
    for item, cache in CACHE:
        cache.clear()


def cache_key(*args, **kwargs):
    key = [(x, type(x)) for x in args]
    if kwargs:
        key.extend([(x, kwargs[x], type(kwargs[x])) for x in sorted(kwargs)])
    key.extend([tuple(global_evaluate)])
    return tuple(key)


USE_CACHE = os.getenv('DIOFANT_USE_CACHE', 'True') == 'True'


def cacheit(f):
    """Caching decorator.

    The result of cached function must be *immutable*.

    Examples
    ========

    >>> @cacheit
    ... def f(a, b):
    ...    print(a, b)
    ...    return a + b

    >>> f(x, y)
    x y
    x + y
    >>> f(x, y)
    x + y
    """

    if USE_CACHE:
        f_cache_it_cache = {}
        CACHE.append((f, f_cache_it_cache))
        return cached(f_cache_it_cache, key=cache_key)(f)
    else:
        return f
