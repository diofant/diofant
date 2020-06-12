"""A module that handles series: find a limit, order the series etc.
"""
from .limits import Limit, limit
from .order import O, Order
from .residues import residue
from .series import series


__all__ = 'Limit', 'O', 'Order', 'limit', 'residue', 'series'
