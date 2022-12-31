"""Some helper functions for the logic module."""

import re

from ..core import Symbol
from . import And, Or


def parse_dimacs(s):
    r"""
    Loads a boolean expression from a string in the DIMACS CNF format.

    The file in the DIMACS CNF format is an ASCII file consisting of
    a two major sections: the preamble and the clauses.

    The preamble contains information about the satisfability problem
    instance, which is contained in lines.  Each line begins with a
    single character (followed by whitespace character): either "c" (a
    comment, rest of line is ignored) or "p" (the problem line).
    The problem line has the following format:

        p dnf VARIABLES CLAUSES

    where VARIABLES and CLAUSES fields contain integers, specifying the
    number of variables (n) and the number of clauses (m), respectively.

    The m clauses appear immediately after the problem line.  The variables
    are assumed to be numbered from 1 up to n.  Each clause should be
    represented by a sequence of numbers, separated by one or more
    whitespace characters.  The non-negated version of a variable i is
    represented by i, the negated version is represented by -i.  Each
    clause is terminated by the value 0.

    Examples
    ========

    >>> parse_dimacs('p cnf 2 1\n1 2 0')
    x1 | x2

    References
    ==========

    * https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps

    """
    pComment = re.compile('^c')
    lines = [_ for _ in s.split('\n') if not pComment.match(_)]

    pStats = re.compile(r'^p\s+cnf\s+(?P<nsymbols>\d+)\s+(?P<nclauses>\d+)$')
    m = pStats.match(lines.pop(0))
    if not m:
        raise ValueError('Bad problem line')

    syms = [Symbol(f'x{i}') for i in range(1, int(m.group('nsymbols')) + 1)]
    clauses = []
    nclauses = int(m.group('nclauses'))
    args = []

    while nclauses:
        try:
            line = lines.pop(0)
        except IndexError as exc:
            raise ValueError('Too few clauses') from exc

        line = re.sub(r'^\s*', '', line)
        line = re.sub(r'\s*$', '', line)
        line = re.sub(r'\s+', ' ', line)

        for lit in line.split(' '):
            try:
                ilit = int(lit)
            except ValueError as exc:
                raise ValueError('Bad clause format') from exc
            if ilit == 0:
                clauses.append(Or(*args))
                nclauses -= 1
                args = []
                continue
            num = abs(ilit)
            sign = ilit >= 0
            sym = syms[num - 1]
            args.append(sym if sign else ~sym)

    return And(*clauses)
