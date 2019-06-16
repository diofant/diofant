"""For reading in DIMACS file format

www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps
"""

import re

from ...core import Symbol
from ..boolalg import And, Or


def load(s):
    r"""Loads a boolean expression from a string.

    Examples
    ========

    >>> load('1')
    cnf_1
    >>> load('1 2')
    cnf_1 | cnf_2
    >>> load('1 \n 2')
    cnf_1 & cnf_2
    >>> load('1 2 \n 3')
    cnf_3 & (cnf_1 | cnf_2)

    """
    clauses = []

    lines = s.split('\n')

    pComment = re.compile('c.*')
    pStats = re.compile(r'p\s*cnf\s*(\d*)\s*(\d*)')

    while len(lines) > 0:
        line = lines.pop(0)

        # Only deal with lines that aren't comments
        if not pComment.match(line):
            m = pStats.match(line)

            if not m:
                nums = line.rstrip('\n').split(' ')
                list = []
                for lit in nums:
                    if lit != '':
                        if int(lit) == 0:
                            continue
                        num = abs(int(lit))
                        sign = True
                        if int(lit) < 0:
                            sign = False

                        if sign:
                            list.append(Symbol("cnf_%s" % num))
                        else:
                            list.append(~Symbol("cnf_%s" % num))

                if len(list) > 0:
                    clauses.append(Or(*list))

    return And(*clauses)
