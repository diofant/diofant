"""
Package for handling logical expressions.
"""

from .boolalg import (ITE, And, Equivalent, Implies, Nand, Nor, Not, Or,
                      POSform, SOPform, Xor, bool_map, false, simplify_logic,
                      to_cnf, to_dnf, to_nnf, true)
from .inference import satisfiable


__all__ = ('ITE', 'And', 'Equivalent', 'Implies', 'Nand', 'Nor', 'Not', 'Or',
           'POSform', 'SOPform', 'Xor', 'bool_map', 'false', 'simplify_logic',
           'to_cnf', 'to_dnf', 'to_nnf', 'true', 'satisfiable')
