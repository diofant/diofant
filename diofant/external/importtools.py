"""Tools to assist importing optional external modules."""

import importlib


def import_module(module):
    """Import and return a module if it is installed."""
    try:
        return importlib.import_module(module)
    except ImportError:
        return
