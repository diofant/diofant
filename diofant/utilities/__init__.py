"""This module contains some general purpose utilities that are used across
Diofant.
"""

from .iterables import (cantor_product, default_sort_key, flatten, group,
                        numbered_symbols, ordered, postorder_traversal, sift,
                        subsets, unflatten)
from .lambdify import lambdify
from .misc import filldedent


__all__ = ('cantor_product', 'default_sort_key', 'flatten',
           'group', 'numbered_symbols',
           'ordered', 'postorder_traversal', 'sift', 'subsets',
           'unflatten', 'lambdify', 'filldedent')
