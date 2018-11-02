#!/usr/bin/env python3
"""
Setup script for Diofant.

This script uses Setuptools (https://setuptools.readthedocs.io/en/latest/).
"""

import re

import setuptools


with open('diofant/__init__.py') as f:
    m = re.search('^__version__ = "([0-9ab.]+(dev[0-9]+)?)".*$', f.read(), re.M)
    if m:
        __version__ = m.group(1)
    else:
        raise RuntimeError("Unable to find version string.")

setuptools.setup(version=__version__)
