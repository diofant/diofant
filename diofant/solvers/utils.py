"""General utility functions for solvers."""

import warnings

from ..core import (expand_mul, expand_multinomial, nan, oo,
                    preorder_traversal, zoo)
from ..core.sympify import sympify
from ..simplify.simplify import posify, simplify


__all__ = 'checksol',


def checksol(f, sol, **flags):
    r"""Checks whether sol is a solution of equations f.

    Examples
    ========

    >>> checksol(x**4 - 1, {x: 1})
    True
    >>> checksol(x**4 - 1, {x: 0})
    False
    >>> checksol(x**2 + y**2 - 5**2, {x: 3, y: 4})
    True

    Returns
    =======

    bool or None
        Return True, if solution satisfy all equations
        in ``f``.  Return False, if a solution doesn't
        satisfy any equation.  Else (i.e. one or more checks
        are inconclusive), return None.

    Parameters
    ==========

    f : Expr or iterable of Expr's
        Equations to substitute solutions in.
    sol : dict of Expr's
        Mapping of symbols to values.
    \*\*flags : dict
        A dictionary of following parameters:

        minimal : bool, optional
            Do a very fast, minimal testing.  Default is False.
        warn : bool, optional
            Show a warning if it could not conclude.  Default is False.
        simplify : bool, optional
            Simplify solution before substituting into function and
            simplify the function before trying specific simplifications.
            Default is True.
        force : bool, optional
           Make positive all symbols without assumptions regarding
           sign.  Default is False.

    """
    minimal = flags.get('minimal', False)

    if not isinstance(sol, dict):
        raise ValueError(f'Expecting dictionary but got {sol}')

    if sol and not f.has(*list(sol)):
        # if f(y) == 0, x=3 does not set f(y) to zero...nor does it not
        if f.is_Number:
            return f.is_zero
        return

    illegal = {nan, zoo, oo, -oo}
    if any(sympify(v).atoms() & illegal for k, v in sol.items()):
        return False

    was = f
    attempt = -1
    while 1:
        attempt += 1
        if attempt == 0:
            val = f.subs(sol)
        elif attempt == 1:
            assert val.free_symbols
            if not val.is_constant(*list(sol), simplify=not minimal):
                return False
            # there are free symbols -- simple expansion might work
            _, val = val.as_content_primitive()
            val = expand_mul(expand_multinomial(val))
        elif attempt == 2:
            if minimal:
                return
            if flags.get('simplify', True):
                for k in sol:
                    sol[k] = simplify(sol[k])
            # start over without the failed expanded form, possibly
            # with a simplified solution
            val = simplify(f.subs(sol))
            if flags.get('force', True):
                val, _ = posify(val)
                # expansion may work now, so try again and check
                exval = expand_mul(expand_multinomial(val))
                if exval.is_number or not exval.free_symbols:
                    # we can decide now
                    val = exval
        else:
            # if there are no radicals and no functions then this can't be
            # zero anymore -- can it?
            pot = preorder_traversal(expand_mul(val))
            seen = set()
            saw_pow_func = False
            for p in pot:
                if p in seen:
                    continue
                seen.add(p)
                if p.is_Pow and not p.exp.is_Integer:
                    saw_pow_func = True
                elif p.is_Function:
                    saw_pow_func = True
                if saw_pow_func:
                    break
            if saw_pow_func is False:
                return False
            break

        if val == was:
            continue
        if val.is_Rational:
            return val == 0
        if val.is_nonzero:
            return False
        if not val.free_symbols:
            return bool(abs(val.evalf(18, strict=False).evalf(12, chop=True)) < 1e-9)
        was = val

    if flags.get('warn', False):
        warnings.warn(f'\n\tWarning: could not verify solution {sol}.')
