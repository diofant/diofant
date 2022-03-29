"""
This module exports all latin and greek letters as Symbols, so you can
conveniently do

    >>> from diofant.abc import x, y

instead of the slightly more clunky-looking

    >>> x, y = symbols('x y')

Caveats
=======

1. As of the time of writing this, the names ``O``, ``S``, ``I``, ``N``,
   ``E``, and ``Q`` are colliding with names defined in Diofant. If you import
   them from both ``diofant.abc`` and ``diofant``, the second import will "win".
   This is an issue only for * imports, which should only be used for
   short-lived code such as interactive sessions and throwaway scripts that do
   not survive until the next Diofant upgrade, where ``diofant`` may contain a
   different set of names.

2. This module does not define symbol names on demand, i.e.
   ``from diofant.abc import foo`` will be reported as an error because
   ``diofant.abc`` does not contain the name ``foo``. To get a symbol named
   `'foo'`, you still need to use ``Symbol('foo')`` or ``symbols('foo')``.
   You can freely mix usage of ``diofant.abc`` and ``Symbol``/``symbols``, though
   sticking with one and only one way to get the symbols does tend to make the
   code more readable.
"""

import string

from .core import Symbol, symbols
from .core.alphabets import greeks


# ##### Symbol definitions #####

a, b, c, d, e, f, g, h, i, j, k, l, m = symbols('a:m')
n, o, p, q, r, s, t, u, v, w, x, y, z = symbols('n:z')

A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A:M')
N, O, P, Q, R, S, T, U, V, W, X, Y, Z = symbols('N:Z')

alpha, beta, gamma, delta = symbols('alpha, beta, gamma, delta')
epsilon, zeta, eta, theta = symbols('epsilon, zeta, eta, theta')
iota, kappa, lamda, mu = symbols('iota, kappa, lamda, mu')
nu, xi, omicron, pi = symbols('nu, xi, omicron, pi')
rho, sigma, tau, upsilon = symbols('rho, sigma, tau, upsilon')
phi, chi, psi, omega = symbols('phi, chi, psi, omega')


# ##### Clashing-symbols diagnostics #####

# We want to know which names in Diofant collide with those in here.
# This is mostly for diagnosing Diofant's namespace during Diofant development.

_latin = list(string.ascii_letters)
# OSINEQ should not be imported as they clash; gamma, pi and zeta clash, too
_greek = list(greeks)  # make a copy, so we can mutate it
# Note: We import lamda since lambda is a reserved keyword in Python
_greek.remove('lambda')
_greek.append('lamda')


def clashing():
    """Return the clashing-symbols dictionaries.

    ``clash1`` defines all the single letter variables that clash with
    Diofant objects; ``clash2`` defines the multi-letter clashing symbols;
    and ``clash`` is the union of both. These can be passed for ``locals``
    during sympification if one desires Symbols rather than the non-Symbol
    objects for those names.

    Examples
    ========

    >>> from diofant.abc import _clash, _clash1, _clash2
    >>> sympify('Q & C', locals=_clash1)
    And(C, Q)
    >>> sympify('pi(x)', locals=_clash2)
    pi(x)
    >>> sympify('pi(C, Q)', locals=_clash)
    pi(C, Q)

    Note: if changes are made to the docstring examples they can only
    be tested after removing "clashing" from the list of deleted items
    at the bottom of this file which removes this function from the
    namespace.

    """
    ns = {}
    exec('from diofant import *', ns)  # pylint: disable=exec-used
    clash1 = {}
    clash2 = {}
    while ns:
        k, _ = ns.popitem()
        if k in _greek:
            clash2[k] = Symbol(k)
            _greek.remove(k)
        elif k in _latin:
            clash1[k] = Symbol(k)
            _latin.remove(k)
    clash = {}
    clash.update(clash1)
    clash.update(clash2)
    return clash1, clash2, clash


_clash1, _clash2, _clash = clashing()

del _latin, _greek, clashing, greeks, symbols, Symbol
