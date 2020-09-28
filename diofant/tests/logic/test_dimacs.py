"""Various tests on satisfiability using dimacs cnf file syntax
You can find lots of cnf files in
ftp://dimacs.rutgers.edu/pub/challenge/satisfiability/benchmarks/cnf/
"""

import os

from diofant.logic.algorithms.dpll import dpll_satisfiable
from diofant.logic.algorithms.dpll2 import \
    dpll_satisfiable as dpll2_satisfiable
from diofant.logic.utilities.dimacs import load


__all__ = ()


def load_file(location):
    """Loads a boolean expression from a file."""
    location = os.path.dirname(__file__) + '/' + location
    with open(location) as f:
        s = f.read()

    return load(s)


def test_f1():
    assert bool(dpll_satisfiable(load_file('simple_v3_c2.cnf')))


def test_f2():
    assert bool(dpll_satisfiable(load_file('quinn.cnf')))


def test_f3():
    assert bool(dpll_satisfiable(load_file('f3.cnf')))


def test_f4():
    assert not bool(dpll_satisfiable(load_file('hole6.cnf')))


def test_f5():
    assert bool(dpll_satisfiable(load_file('f5.cnf')))


def test_f6():
    assert not bool(dpll2_satisfiable(load_file('aim-50-2_0-no-2.cnf')))
