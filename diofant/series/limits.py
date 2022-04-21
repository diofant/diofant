from ..core import (Dummy, Expr, Float, Integer, PoleError, Rational, Symbol,
                    nan, oo, sympify)
from ..core.function import UndefinedFunction
from ..functions import Abs, cos, sign, sin
from .gruntz import limitinf
from .order import Order


def limit(expr, z, z0, dir=-1):
    """
    Compute the directional limit of ``expr`` at the point ``z0``.

    Examples
    ========

    >>> limit(sin(x)/x, x, 0)
    1
    >>> limit(1/x, x, 0)
    oo
    >>> limit(1/x, x, 0, dir=1)
    -oo
    >>> limit(1/x, x, oo)
    0

    See Also
    ========

    Limit

    """
    return Limit(expr, z, z0, dir).doit(deep=False)


def heuristics(e, z, z0, dir):
    e = e.expand()
    if (e.is_Mul or e.is_Add or e.is_Pow or
            (e.is_Function and not isinstance(e.func, UndefinedFunction))):
        r = []
        for a in e.args:
            l = limit(a, z, z0, dir)
            if (l.has(oo) and (l.func not in (sin, cos) and
                               l.is_finite is None)) or isinstance(l, Limit):
                return
            r.append(l)

        if (rv := e.func(*r)) is not nan:
            return rv


class Limit(Expr):
    r"""
    Represents a directional limit of ``expr`` at the point ``z0``.

    Parameters
    ==========

    expr : Expr
        algebraic expression
    z    : Symbol
        variable of the ``expr``
    z0   : Expr
        limit point, `z_0`
    dir  : {-1, 1, "real"}, optional
        For ``dir=-1`` (default) it calculates the limit from the right
        (`z\to z_0 + 0`) and for ``dir=1`` the limit from the left (`z\to
        z_0 - 0`).  If ``dir="real"``, the limit is the bidirectional real
        limit.  For infinite ``z0`` (``oo`` or ``-oo``), the ``dir`` argument
        is determined from the direction of the infinity (i.e.,
        ``dir=1`` for ``oo``).

    Examples
    ========

    >>> Limit(sin(x)/x, x, 0)
    Limit(sin(x)/x, x, 0)
    >>> Limit(1/x, x, 0, dir=1)
    Limit(1/x, x, 0, dir=1)

    """

    def __new__(cls, e, z, z0, dir=-1):
        e, z, z0 = map(sympify, [e, z, z0])

        if z0 is oo:
            dir = Integer(+1)
        elif z0 == -oo:
            dir = Integer(-1)

        if isinstance(dir, str):
            dir = Symbol(dir)
        else:
            dir = sympify(dir)

        if str(dir) not in ('1', '-1', 'real'):
            raise ValueError(f"direction must be '±1' or 'real', not {dir}")

        obj = super().__new__(cls)
        obj._args = (e, z, z0, dir)
        return obj

    @property
    def free_symbols(self):
        e, z, z0 = self.args[:3]
        return (e.free_symbols - z.free_symbols) | z0.free_symbols

    def doit(self, **hints):
        """
        Evaluates limit.

        Notes
        =====

        First we handle some trivial cases (i.e. constant), then try
        Gruntz algorithm (see the :py:mod:`~diofant.series.gruntz` module).

        """
        e, z, z0, dir = self.args

        if hints.get('deep', True):
            e = e.doit(**hints)
            z = z.doit(**hints)
            z0 = z0.doit(**hints)

        if str(dir) == 'real':
            right = limit(e, z, z0)
            left = limit(e, z, z0, 1)
            if not left.equals(right):
                raise PoleError(f'left and right limits for expression {e} at '
                                f'point {z}={z0} seems to be not equal')
            return right

        if (has_Floats := e.has(Float) or z0.has(Float)):
            e = e.subs({k: Rational(k) for k in e.atoms(Float)})
            z0 = z0.subs({k: Rational(k) for k in z0.atoms(Float)})

        if z0.has(z):
            newz = z.as_dummy()
            r = limit(e.subs({z: newz}), newz, z0, dir)
            if isinstance(r, Limit):
                r = r.subs({newz: z})
            return r

        if not e.has(z):
            return e

        if e.has(Order) and (order := e.getO()) and (z, z0) in order.args[1:]:
            order = limit(order.expr, z, z0, dir)
            e = e.removeO() + order

        # Convert to the limit z->oo and use Gruntz algorithm.
        if z0 == -oo:
            e = e.subs({z: -z})
        elif z0 != oo:
            e = e.subs({z: z0 - dir/z})

        # We need a fresh variable with correct assumptions.
        newz = Dummy(z.name, positive=True, finite=True)
        e = e.subs({z: newz})

        if e.is_Boolean or e.is_Relational:
            try:
                has_oo = e.as_set().closure.contains(oo)
            except NotImplementedError:
                return self
            if has_oo.is_Boolean:
                return has_oo
            raise NotImplementedError

        def tr_abs(f):
            s = sign(limit(f.args[0], newz, oo))
            return s*f.args[0] if s in (1, -1) else f

        def tr_Piecewise(f):
            for a, c in f.args:
                if not c.is_Atom:
                    c = c.as_set().closure.contains(oo)
                    if not c.is_Atom:
                        raise NotImplementedError("Parametric limits aren't "
                                                  'supported yet.')
                    if c:
                        break
            return a

        e = e.replace(lambda f: isinstance(f, Abs) and f.has(newz), tr_abs)
        e = e.replace(lambda f: f.is_Piecewise and f.has(newz), tr_Piecewise)

        try:
            r = limitinf(e, newz)
        except (PoleError, ValueError, NotImplementedError):
            r = None
            if hints.get('heuristics', True):
                r = heuristics(*self.args)
            if r is None:
                return self

        if has_Floats:
            r = r.evalf()

        return r
