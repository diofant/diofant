"""Some calculus-related methods waiting to find a better place in the
Diofant modules tree.
"""

from .optimization import maximize, minimize
from .residues import residue


__all__ = 'maximize', 'minimize', 'residue'
