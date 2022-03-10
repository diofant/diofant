"""
Various tests on satisfiability using the DIMACS CNF file format.

References
==========

* http://archive.dimacs.rutgers.edu/pub/challenge/satisfiability/benchmarks/cnf/

"""

import os

import pytest

from diofant.logic import satisfiable
from diofant.logic.utils import parse_dimacs


__all__ = ()


def load_file(location):
    """Loads a boolean expression from a file."""
    location = os.path.dirname(__file__) + '/' + location + '.cnf'
    with open(location, encoding='utf-8') as f:
        s = f.read()

    return parse_dimacs(s)


def test_parse_dimacs():
    pytest.raises(ValueError, lambda: parse_dimacs('q 1 2 3'))
    pytest.raises(ValueError, lambda: parse_dimacs('p cnf 2 1\n1 ?'))
    pytest.raises(ValueError, lambda: parse_dimacs('p cnf 2 2\n 1'))


@pytest.mark.parametrize('algorithm', ['dpll', 'dpll2'])
@pytest.mark.parametrize('name', ['f3', 'f5', 'quinn', 'simple_v3_c2'])
def test_dimacs_satisfiable(algorithm, name):
    assert bool(satisfiable(load_file(name), algorithm=algorithm))


@pytest.mark.parametrize('algorithm', ['dpll', 'dpll2'])
@pytest.mark.parametrize('name', ['aim-50-2_0-no-2', 'hole6'])
def test_dimacs_not_satisfiable(algorithm, name):
    assert not bool(satisfiable(load_file(name), algorithm=algorithm))
