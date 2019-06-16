"""This module contains some general purpose utilities that are used across
Diofant.
"""

from .iterables import (cantor_product, capture, default_sort_key, dict_merge,
                        flatten, group, has_dups, has_variety,
                        numbered_symbols, ordered, postfixes,
                        postorder_traversal, prefixes, reshape, sift, subsets,
                        topological_sort, unflatten, variations)
from .lambdify import lambdify
from .misc import filldedent
