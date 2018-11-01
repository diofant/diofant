#!/usr/bin/env python3
"""
Setup script for Diofant.

This script uses Setuptools (https://setuptools.readthedocs.io/en/latest/),
a collection of enhancements to the standard Python distutils.
"""

import re

import setuptools


with open('diofant/__init__.py') as f:
    source = f.read()
    for line in source.splitlines():
        m = re.match('^__version__ = "([0-9ab.]+(dev[0-9]+)?)".*$', line)
        if m:
            __version__ = m.group(1)

setuptools.setup(version=__version__)
