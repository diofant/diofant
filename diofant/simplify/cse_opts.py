"""Optimizations of the expression tree representation for better CSE
opportunities.
"""

from ..core import Add, Basic, Integer, Mul, preorder_traversal
from ..utilities.iterables import default_sort_key


def sub_pre(e):
    """Replace y - x with -(x - y) if -1 can be extracted from y - x."""
    reps = [a for a in e.atoms(Add) if a.could_extract_minus_sign()]

    # make it canonical
    reps.sort(key=default_sort_key)

    e = e.xreplace({a: Mul._from_args([Integer(-1), -a]) for a in reps})
    # repeat again for persisting Adds but mark these with a leading 1, -1
    # e.g. y - x -> 1*-1*(x - y)
    assert isinstance(e, Basic)
    negs = {}
    for a in sorted(e.atoms(Add), key=default_sort_key):
        if a in reps or a.could_extract_minus_sign():
            negs[a] = Mul._from_args([Integer(1), Integer(-1), -a])
    return e.xreplace(negs)


def sub_post(e):
    """Replace 1*-1*x with -x."""
    replacements = []
    for node in preorder_traversal(e):
        if isinstance(node, Mul) and node.args[0] == 1 and node.args[1] == -1:
            replacements.append((node, -Mul._from_args(node.args[2:])))
    for node, replacement in replacements:
        e = e.xreplace({node: replacement})

    return e
