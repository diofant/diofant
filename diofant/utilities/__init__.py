"""This module contains some general purpose utilities that are used across
Diofant.
"""

from .iterables import (flatten, group, take, subsets,
                        variations, numbered_symbols, capture, dict_merge,
                        postorder_traversal, interactive_traversal,
                        prefixes, postfixes, sift, topological_sort, unflatten,
                        has_dups, has_variety, reshape, default_sort_key, ordered)

from .misc import filldedent

from .lambdify import lambdify

from .decorator import public
