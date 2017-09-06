"""
Package for handling logical expressions.
"""

from .boolalg import (to_cnf, to_dnf, to_nnf, And, Or, Not, Xor, Nand, Nor, Implies,  # noqa: F401
                      Equivalent, ITE, POSform, SOPform, simplify_logic, bool_map,
                      true, false)
from .inference import satisfiable  # noqa: F401
