"""This module contains some general purpose utilities that are used across
Diofant.
"""

from .iterables import (flatten, group, take, subsets,  # noqa: F401
                        variations, numbered_symbols, capture, dict_merge,
                        postorder_traversal, prefixes, postfixes, sift,
                        topological_sort, unflatten, has_dups, has_variety,
                        reshape, default_sort_key, ordered, cantor_product)

from .misc import filldedent  # noqa: F401

from .lambdify import lambdify  # noqa: F401
