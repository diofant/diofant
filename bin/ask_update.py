#!/usr/bin/env python

""" Update the ask_generated.py file

This must be run each time known_facts is changed

Should be run from sympy root directory

$ python bin/ask_update.py
"""

# hook in-tree SymPy into Python path, if possible
import os
import sys

ask_path = os.path.abspath(__file__)
ask_dir = os.path.dirname(ask_path)
sympy_top = os.path.split(ask_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')

if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

from sympy.assumptions.ask import (compute_known_facts, known_facts,
        known_facts_keys)

with open('sympy/assumptions/ask_generated.py', 'w') as f:
    code = compute_known_facts(known_facts, known_facts_keys)
    f.write(code)
