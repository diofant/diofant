"""For reading in DIMACS file format

www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/satformat.ps

"""

import re

from diofant.core import Symbol
from diofant.logic.boolalg import And, Or


def load(s):
    """Loads a boolean expression from a string.

    Examples
    ========

    >>> from diofant.logic.utilities.dimacs import load
    >>> load('1')
    cnf_1
    >>> load('1 2')
    Or(cnf_1, cnf_2)
    >>> load('1 \\n 2')
    And(cnf_1, cnf_2)
    >>> load('1 2 \\n 3')
    And(Or(cnf_1, cnf_2), cnf_3)
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


def load_file(location):
    """Loads a boolean expression from a file."""
    with open(location) as f:
        s = f.read()

    return load(s)
