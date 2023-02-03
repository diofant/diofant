"""This module contains some general purpose utilities that are used across
Diofant.
"""

from .iterables import (cantor_product, default_sort_key, flatten, group,
                        has_dups, has_variety, numbered_symbols, ordered,
                        postorder_traversal, sift, subsets, topological_sort,
                        unflatten)
from .lambdify import lambdify
from .misc import filldedent


__all__ = ('cantor_product', 'default_sort_key', 'flatten',
           'group', 'has_dups', 'has_variety', 'numbered_symbols',
           'ordered', 'postorder_traversal', 'sift', 'subsets',
           'topological_sort', 'unflatten', 'lambdify', 'filldedent')
