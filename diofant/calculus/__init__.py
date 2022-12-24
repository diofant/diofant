"""Some calculus-related methods waiting to find a better place in the
Diofant modules tree.
"""

from .optimization import maximize, minimize
from .order import O, Order
from .residues import residue


__all__ = 'O', 'Order', 'maximize', 'minimize', 'residue'
