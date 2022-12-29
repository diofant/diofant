from ..core import Dummy, Expr, Float, PoleError, Rational, nan, oo, sympify
from ..core.function import UndefinedFunction
from ..sets import Reals
from .gruntz import limitinf
from .order import Order


def limit(expr, z, z0, dir=None):
    """
    Compute the directional limit of ``expr`` at the point ``z0``.

    Examples
    ========

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
    from ..functions import cos, sin

    e = e.expand()
    if (e.is_Mul or e.is_Add or e.is_Pow or
            (e.is_Function and not e.is_Piecewise and
             not isinstance(e.func, UndefinedFunction))):
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
    Represents an unevaluated directional limit of ``expr`` at the point ``z0``.

    Parameters
    ==========

    expr : Expr
        algebraic expression
    z    : Symbol
        variable of the ``expr``
    z0   : Expr
        limit point, `z_0`
    dir  : Expr or Reals, optional
        selects the direction (as ``sign(dir)``) to approach the limit point
        if the ``dir`` is an Expr.  For infinite ``z0``, the default value
        is determined from the direction of the infinity (e.g., the limit
        from the left, ``dir=1``, for ``oo``).  Otherwise, the default is
        the limit from the right, ``dir=-1``.   If ``dir=Reals``, the limit
        is the bidirectional real limit.

    Examples
    ========

    >>> Limit(1/x, x, 0, dir=1)
    Limit(1/x, x, 0, dir=1)
    >>> _.doit()
    -oo

    See Also
    ========

    limit

    """

    def __new__(cls, e, z, z0, dir=None):
        from ..functions import sign

        e, z, z0, dir = map(sympify, [e, z, z0, dir])

        if z0.is_infinite:
            dir = sign(z0).simplify()
        elif dir is None:
            dir = Rational(-1)

        if dir == Reals:
            pass
        elif isinstance(dir, Expr) and dir.is_nonzero:
            dir = dir/abs(dir)
        else:
            raise ValueError('direction must be either a nonzero expression '
                             f'or Reals, not {dir}')

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
        Gruntz algorithm (see the :py:mod:`~diofant.calculus.gruntz` module).

        """
        e, z, z0, dir = self.args

        if hints.get('deep', True):
            e = e.doit(**hints)
            z = z.doit(**hints)
            z0 = z0.doit(**hints)

        if dir == Reals:
            right = limit(e, z, z0)
            left = limit(e, z, z0, 1)
            if not left.equals(right):
                raise PoleError(f'left and right limits for the expression {e} '
                                f'at point {z}={z0} seems to be not equal')
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

        # Use a fresh variable to remove assumptions on the dummy variable.
        newz = Dummy('z')
        e, z = e.subs({z: newz}), newz

        if not e.has(z):
            return e

        if e.has(Order) and (order := e.getO()) and (z, z0) in order.args[1:]:
            order = limit(order.expr, z, z0, dir)
            e = e.removeO() + order

        # Convert to the limit z->oo and use Gruntz algorithm.
        e = e.subs({z: dir*z})
        z0 = z0/dir
        if z0 != oo:
            e = e.subs({z: z0 - 1/z})

        # We need a fresh variable with correct assumptions.
        newz = Dummy(z.name, positive=True, finite=True)
        e, z = e.subs({z: newz}), newz

        if e.is_Boolean or e.is_Relational:
            try:
                has_oo = e.as_set().closure.contains(oo)
            except NotImplementedError:
                return self
            if has_oo.is_Boolean:
                return has_oo
            raise NotImplementedError

        try:
            r = limitinf(e, z)
        except (PoleError, ValueError, NotImplementedError):
            r = None
            if hints.get('heuristics', True):
                r = heuristics(*self.args)
            if r is None:
                return self

        if has_Floats:
            r = r.evalf()

        return r
