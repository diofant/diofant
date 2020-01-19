"""This module contains some general purpose utilities that are used across
Diofant.
"""

from .iterables import (cantor_product, capture, default_sort_key, flatten,
                        group, has_dups, has_variety, numbered_symbols,
                        ordered, postfixes, postorder_traversal, prefixes,
                        sift, subsets, topological_sort, unflatten, variations)
from .lambdify import lambdify
from .misc import filldedent


__all__ = ('cantor_product', 'capture', 'default_sort_key', 'flatten',
           'group', 'has_dups', 'has_variety', 'numbered_symbols',
           'ordered', 'postfixes', 'postorder_traversal', 'prefixes',
           'sift', 'subsets', 'topological_sort', 'unflatten', 'variations',
           'lambdify', 'filldedent')
