import math

from ..core import (Add, Dummy, Expr, Integer, Mul, Rational, count_ops,
                    expand_mul, factor_terms)
from ..core.function import _mexpand
from ..core.sympify import sympify
from ..functions import log, root, sign, sqrt
from ..polys import Poly, PolynomialError, cancel, degree
from ..utilities import default_sort_key, ordered
from .powsimp import powdenest


def is_sqrt(expr):
    """Return True if expr is a sqrt, otherwise False."""
    return expr.is_Pow and expr.exp.is_Rational and abs(expr.exp) == Rational(1, 2)


def sqrt_depth(p):
    """Return the maximum depth of any square root argument of p.

    Neither of these square roots contains any other square roots
    so the depth is 1:

    >>> sqrt_depth(1 + sqrt(2)*(1 + sqrt(3)))
    1

    The sqrt(3) is contained within a square root so the depth is
    2:

    >>> sqrt_depth(1 + sqrt(2)*sqrt(1 + sqrt(3)))
    2

    """
    if p.is_Atom:
        return 0
    elif p.is_Add or p.is_Mul:
        return max((sqrt_depth(x) for x in p.args), key=default_sort_key)
    elif is_sqrt(p):
        return sqrt_depth(p.base) + 1
    else:
        return 0


def is_algebraic(p):
    """Return True if p is comprised of only Rationals or square roots
    of Rationals and algebraic operations.

    Examples
    ========

    >>> is_algebraic(sqrt(2)*(3/(sqrt(7) + sqrt(5)*sqrt(2))))
    True
    >>> is_algebraic(sqrt(2)*(3/(sqrt(7) + sqrt(5)*cos(2))))
    False

    """
    if p.is_Rational:
        return True
    elif p.is_Atom:
        return False
    elif is_sqrt(p) or p.is_Pow and p.exp.is_Integer:
        return is_algebraic(p.base)
    elif p.is_Add or p.is_Mul:
        return all(is_algebraic(x) for x in p.args)
    else:
        return False


def _subsets(n):
    """
    Returns all possible subsets of the set (0, 1, ..., n-1) except the
    empty set, listed in reversed lexicographical order according to binary
    representation, so that the case of the fourth root is treated last.

    Examples
    ========

    >>> _subsets(2)
    [[1, 0], [0, 1], [1, 1]]

    """
    if n == 1:
        a = [[1]]
    elif n == 2:
        a = [[1, 0], [0, 1], [1, 1]]
    elif n == 3:
        a = [[1, 0, 0], [0, 1, 0], [1, 1, 0],
             [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    else:
        b = _subsets(n - 1)
        a0 = [x + [0] for x in b]
        a1 = [x + [1] for x in b]
        a = a0 + [[0]*(n - 1) + [1]] + a1
    return a


def sqrtdenest(expr, max_iter=3):
    """Denests sqrts in an expression that contain other square roots
    if possible, otherwise returns the expr unchanged. This is based on the
    algorithms of [1].

    Examples
    ========

    >>> sqrtdenest(sqrt(5 + 2 * sqrt(6)))
    sqrt(2) + sqrt(3)

    See Also
    ========

    unrad

    References
    ==========

    * https://researcher.watson.ibm.com/researcher/files/us-fagin/symb85.pdf

    * D. J. Jeffrey and A. D. Rich, 'Symplifying Square Roots of Square Roots
      by Denesting' (available at http://www.cybertester.com/data/denest.pdf)

    """
    expr = expand_mul(sympify(expr))
    for i in range(max_iter):
        z = _sqrtdenest0(expr)
        if expr == z:
            return expr
        expr = z
    return expr


def _sqrt_match(p):
    """Return [a, b, r] for p.match(a + b*sqrt(r)) where, in addition to
    matching, sqrt(r) also has then maximal sqrt_depth among addends of p.

    Examples
    ========

    >>> _sqrt_match(1 + sqrt(2) + sqrt(2)*sqrt(3) + 2*sqrt(1+sqrt(5)))
    [1 + sqrt(2) + sqrt(6), 2, 1 + sqrt(5)]

    """
    from .radsimp import split_surds

    p = _mexpand(p)
    if p.is_Number:
        res = (p, Integer(0), Integer(0))
    elif p.is_Add:
        pargs = sorted(p.args, key=default_sort_key)
        if all((x**2).is_Rational for x in pargs):
            r, b, a = split_surds(p)
            res = a, b, r
            return list(res)
        # to make the process canonical, the argument is included in the tuple
        # so when the max is selected, it will be the largest arg having a
        # given depth
        v = [(sqrt_depth(x), x, i) for i, x in enumerate(pargs)]
        nmax = max(v, key=default_sort_key)
        if nmax[0] == 0:
            res = []
        else:
            # select r
            depth, _, i = nmax
            r = pargs.pop(i)
            v.pop(i)
            b = Integer(1)
            if r.is_Mul:
                bv = []
                rv = []
                for x in r.args:
                    if sqrt_depth(x) < depth:
                        bv.append(x)
                    else:
                        rv.append(x)
                b = Mul._from_args(bv)
                r = Mul._from_args(rv)
            # collect terms comtaining r
            a1 = []
            b1 = [b]
            for x in v:
                if x[0] < depth:
                    a1.append(x[1])
                else:
                    x1 = x[1]
                    if x1 == r:
                        b1.append(1)
                    else:
                        if x1.is_Mul:
                            x1args = list(x1.args)
                            if r in x1args:
                                x1args.remove(r)
                                b1.append(Mul(*x1args))
                            else:
                                a1.append(x[1])
                        else:
                            a1.append(x[1])
            a = Add(*a1)
            b = Add(*b1)
            res = (a, b, r**2)
    else:
        b, r = p.as_coeff_Mul()
        if is_sqrt(r):
            res = (Integer(0), b, r**2)
        else:
            res = []
    return list(res)


class SqrtdenestStopIteration(StopIteration):
    """Raised when sqrtdenest algorithm can't denest an expression."""


def _sqrtdenest0(expr):
    """Returns expr after denesting its arguments."""
    if is_sqrt(expr):
        n, d = expr.as_numer_denom()
        if d == 1:  # n is a square root
            if n.base.is_Add:
                args = sorted(n.base.args, key=default_sort_key)
                if len(args) > 2 and all((x**2).is_Integer for x in args):
                    try:
                        return _sqrtdenest_rec(n)
                    except SqrtdenestStopIteration:
                        pass
                expr = sqrt(_mexpand(Add(*[_sqrtdenest0(x) for x in args])))
            return _sqrtdenest1(expr)
        else:
            n, d = [_sqrtdenest0(i) for i in (n, d)]
            return n/d
    if isinstance(expr, Expr):
        args = expr.args
        if args:
            return expr.func(*[_sqrtdenest0(a) for a in args])
    return expr


def _sqrtdenest_rec(expr):
    """Helper that denests the square root of three or more surds.

    It returns the denested expression; if it cannot be denested it
    throws SqrtdenestStopIteration

    Algorithm: expr.base is in the extension Q_m = Q(sqrt(r_1),..,sqrt(r_k));
    split expr.base = a + b*sqrt(r_k), where `a` and `b` are on
    Q_(m-1) = Q(sqrt(r_1),..,sqrt(r_(k-1))); then a**2 - b**2*r_k is
    on Q_(m-1); denest sqrt(a**2 - b**2*r_k) and so on.
    See [1], section 6.

    Examples
    ========

    >>> _sqrtdenest_rec(sqrt(-72*sqrt(2) + 158*sqrt(5) + 498))
    -sqrt(10) + sqrt(2) + 9 + 9*sqrt(5)
    >>> w = -6*sqrt(55)-6*sqrt(35)-2*sqrt(22)-2*sqrt(14)+2*sqrt(77)+6*sqrt(10)+65
    >>> _sqrtdenest_rec(sqrt(w))
    -sqrt(11) - sqrt(7) + sqrt(2) + 3*sqrt(5)

    """
    from .radsimp import rad_rationalize, radsimp, split_surds
    if not expr.is_Pow:
        return sqrtdenest(expr)
    if expr.base < 0:
        return sqrt(-1)*_sqrtdenest_rec(sqrt(-expr.base))
    g, a, b = split_surds(expr.base)
    a = a*sqrt(g)
    if a < b:
        a, b = b, a
    c2 = _mexpand(a**2 - b**2)
    if len(c2.args) > 2:
        g, a1, b1 = split_surds(c2)
        a1 = a1*sqrt(g)
        if a1 < b1:
            a1, b1 = b1, a1
        c2_1 = _mexpand(a1**2 - b1**2)
        c_1 = _sqrtdenest_rec(sqrt(c2_1))
        d_1 = _sqrtdenest_rec(sqrt(a1 + c_1))
        num, den = rad_rationalize(b1, d_1)
        c = _mexpand(d_1/sqrt(2) + num/(den*sqrt(2)))
    else:
        c = _sqrtdenest1(sqrt(c2))

    if sqrt_depth(c) > 1:
        raise SqrtdenestStopIteration
    ac = a + c
    if len(ac.args) >= len(expr.args):
        if count_ops(ac) >= count_ops(expr.base):
            raise SqrtdenestStopIteration
    d = sqrtdenest(sqrt(ac))
    if sqrt_depth(d) > 1:
        raise SqrtdenestStopIteration
    num, den = rad_rationalize(b, d)
    r = d/sqrt(2) + num/(den*sqrt(2))
    r = radsimp(r)
    return _mexpand(r)


def _sqrtdenest1(expr, denester=True):
    """Return denested expr after denesting with simpler methods or, that
    failing, using the denester.

    """
    from .simplify import radsimp

    if not is_sqrt(expr):
        return expr

    a = expr.base
    if a.is_Atom:
        return expr
    val = _sqrt_match(a)
    if not val:
        return expr

    a, b, r = val
    # try a quick numeric denesting
    d2 = _mexpand(a**2 - b**2*r)
    if d2.is_Rational:
        if d2.is_positive:
            z = _sqrt_numeric_denest(a, b, r, d2)
            if z is not None:
                return z
        else:
            # fourth root case
            # sqrtdenest(sqrt(3 + 2*sqrt(3))) =
            # sqrt(2)*3**(1/4)/2 + sqrt(2)*3**(3/4)/2
            dr2 = _mexpand(-d2*r)
            dr = sqrt(dr2)
            if dr.is_Rational:
                z = _sqrt_numeric_denest(_mexpand(b*r), a, r, dr2)
                if z is not None:
                    return z/root(r, 4)

    else:
        z = _sqrt_symbolic_denest(a, b, r)
        if z is not None:
            return z

    if not denester or not is_algebraic(expr):
        return expr

    res = sqrt_biquadratic_denest(expr, a, b, r, d2)
    if res:
        return res

    # now call to the denester
    av0 = [a, b, r, d2]
    z = _denester([radsimp(expr**2)], av0, 0, sqrt_depth(expr))[0]
    if av0[1] is None:
        return expr
    if z is not None:
        if sqrt_depth(z) == sqrt_depth(expr) and count_ops(z) > count_ops(expr):
            return expr
        return z
    return expr


def _sqrt_symbolic_denest(a, b, r):
    """Given an expression, sqrt(a + b*sqrt(b)), return the denested
    expression or None.

    Algorithm:
    If r = ra + rb*sqrt(rr), try replacing sqrt(rr) in ``a`` with
    (y**2 - ra)/rb, and if the result is a quadratic, ca*y**2 + cb*y + cc, and
    (cb + b)**2 - 4*ca*cc is 0, then sqrt(a + b*sqrt(r)) can be rewritten as
    sqrt(ca*(sqrt(r) + (cb + b)/(2*ca))**2).

    Examples
    ========

    >>> a, b, r = 16 - 2*sqrt(29), 2, -10*sqrt(29) + 55
    >>> _sqrt_symbolic_denest(a, b, r)
    sqrt(-2*sqrt(29) + 11) + sqrt(5)

    If the expression is numeric, it will be simplified:

    >>> w = sqrt(sqrt(sqrt(3) + 1) + 1) + 1 + sqrt(2)
    >>> sqrtdenest(sqrt((w**2).expand()))
    1 + sqrt(2) + sqrt(1 + sqrt(1 + sqrt(3)))

    Otherwise, it will only be simplified if assumptions allow:

    >>> w = w.subs({sqrt(3): sqrt(x + 3)})
    >>> sqrtdenest(sqrt((w**2).expand()))
    sqrt((sqrt(sqrt(sqrt(x + 3) + 1) + 1) + 1 + sqrt(2))**2)

    Notice that the argument of the sqrt is a square. If x is made positive
    then the sqrt of the square is resolved:

    >>> _.subs({x: Symbol('x', positive=True)})
    sqrt(sqrt(sqrt(x + 3) + 1) + 1) + 1 + sqrt(2)

    """
    a, b, r = map(sympify, (a, b, r))
    rval = _sqrt_match(r)
    if not rval:
        return
    ra, rb, rr = rval
    if rb:
        y = Dummy('y', positive=True)
        try:
            newa = Poly(a.subs({sqrt(rr): (y**2 - ra)/rb}), y)
        except PolynomialError:
            return
        if newa.degree() == 2:
            cc, cb, ca = newa.all_coeffs()
            cb += b
            if _mexpand(cb**2 - 4*ca*cc).equals(0):
                z = sqrt(ca*(sqrt(r) + cb/(2*ca))**2)
                if z.is_number:
                    z = _mexpand(Mul._from_args(z.as_content_primitive()))
                return z


def _sqrt_numeric_denest(a, b, r, d2):
    """Helper that denest expr = a + b*sqrt(r), with d2 = a**2 - b**2*r > 0
    or returns None if not denested.

    """
    from .simplify import radsimp
    depthr = sqrt_depth(r)
    d = sqrt(d2)
    vad = a + d
    # sqrt_depth(res) <= sqrt_depth(vad) + 1
    # sqrt_depth(expr) = depthr + 2
    # there is denesting if sqrt_depth(vad)+1 < depthr + 2
    # if vad**2 is Number there is a fourth root
    if sqrt_depth(vad) < depthr + 1 or (vad**2).is_Rational:
        vad1 = radsimp(1/vad)
        return (sqrt(vad/2) + sign(b)*sqrt((b**2*r*vad1/2).expand())).expand()


def sqrt_biquadratic_denest(expr, a, b, r, d2):
    """Denest expr = sqrt(a + b*sqrt(r))
    where a, b, r are linear combinations of square roots of
    positive rationals on the rationals (SQRR) and r > 0, b != 0,
    d2 = a**2 - b**2*r > 0.

    If it cannot denest it returns None.

    ALGORITHM
    Search for a solution A of type SQRR of the biquadratic equation
    4*A**4 - 4*a*A**2 + b**2*r = 0                               (1)
    sqd = sqrt(a**2 - b**2*r)
    Choosing the sqrt to be positive, the possible solutions are
    A = sqrt(a/2 +/- sqd/2)
    Since a, b, r are SQRR, then a**2 - b**2*r is a SQRR,
    so if sqd can be denested, it is done by
    _sqrtdenest_rec, and the result is a SQRR.
    Similarly for A.
    Examples of solutions (in both cases a and sqd are positive):

      Example of expr with solution sqrt(a/2 + sqd/2) but not
      solution sqrt(a/2 - sqd/2):
      expr = sqrt(-sqrt(15) - sqrt(2)*sqrt(-sqrt(5) + 5) - sqrt(3) + 8)
      a = -sqrt(15) - sqrt(3) + 8; sqd = -2*sqrt(5) - 2 + 4*sqrt(3)

      Example of expr with solution sqrt(a/2 - sqd/2) but not
      solution sqrt(a/2 + sqd/2):
      w = 2 + r2 + r3 + (1 + r3)*sqrt(2 + r2 + 5*r3)
      expr = sqrt((w**2).expand())
      a = 4*sqrt(6) + 8*sqrt(2) + 47 + 28*sqrt(3)
      sqd = 29 + 20*sqrt(3)

    Define B = b/2*A; eq.(1) implies a = A**2 + B**2*r; then
    expr**2 = a + b*sqrt(r) = (A + B*sqrt(r))**2

    Examples
    ========

    >>> z = sqrt((2*sqrt(2) + 4)*sqrt(2 + sqrt(2)) + 5*sqrt(2) + 8)
    >>> a, b, r = _sqrt_match(z**2)
    >>> d2 = a**2 - b**2*r
    >>> sqrt_biquadratic_denest(z, a, b, r, d2)
    sqrt(2) + sqrt(sqrt(2) + 2) + 2

    """
    from .radsimp import rad_rationalize, radsimp
    if r <= 0 or d2 < 0 or not b or sqrt_depth(expr.base) < 2:
        return
    for x in (a, b, r):
        for y in x.args:
            y2 = y**2
            if not y2.is_Integer or not y2.is_positive:
                return
    sqd = _mexpand(sqrtdenest(sqrt(radsimp(d2))))
    if sqrt_depth(sqd) > 1:
        return
    x1, x2 = [a/2 + sqd/2, a/2 - sqd/2]
    # look for a solution A with depth 1
    for x in (x1, x2):
        A = sqrtdenest(sqrt(x))
        if sqrt_depth(A) > 1:
            continue
        Bn, Bd = rad_rationalize(b, _mexpand(2*A))
        B = Bn/Bd
        z = A + B*sqrt(r)
        if z < 0:
            z = -z
        return _mexpand(z)


def _denester(nested, av0, h, max_depth_level):
    """Denests a list of expressions that contain nested square roots.

    Algorithm based on <http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf>.

    It is assumed that all of the elements of 'nested' share the same
    bottom-level radicand. (This is stated in the paper, on page 177, in
    the paragraph immediately preceding the algorithm.)

    When evaluating all of the arguments in parallel, the bottom-level
    radicand only needs to be denested once. This means that calling
    _denester with x arguments results in a recursive invocation with x+1
    arguments; hence _denester has polynomial complexity.

    However, if the arguments were evaluated separately, each call would
    result in two recursive invocations, and the algorithm would have
    exponential complexity.

    This is discussed in the paper in the middle paragraph of page 179.

    """
    from .simplify import radsimp
    if h > max_depth_level:
        return None, None
    if av0[1] is None:
        return None, None
    if (av0[0] is None and
            all(n.is_Number for n in nested)):  # no arguments are nested
        for f in _subsets(len(nested)):  # test subset 'f' of nested
            p = _mexpand(Mul(*[nested[i] for i in range(len(f)) if f[i]]))
            if f.count(1) > 1 and f[-1]:
                p = -p
            sqp = sqrt(p)
            if sqp.is_Rational:
                return sqp, f  # got a perfect square so return its square root.
        # Otherwise, return the radicand from the previous invocation.
        return sqrt(nested[-1]), [0]*len(nested)
    else:
        R = None
        if av0[0] is not None:
            values = [av0[:2]]
            R = av0[2]
            nested2 = [av0[3], R]
            av0[0] = None
        else:
            values = list(filter(None, (_sqrt_match(expr) for expr in nested)))
            for v in values:
                if v[2]:  # Since if b=0, r is not defined
                    if R is not None:
                        if R != v[2]:
                            av0[1] = None
                            return None, None
                    else:
                        R = v[2]
            if R is None:
                # return the radicand from the previous invocation
                return sqrt(nested[-1]), [0]*len(nested)
            nested2 = [_mexpand(v[0]**2) -
                       _mexpand(R*v[1]**2) for v in values] + [R]
        d, f = _denester(nested2, av0, h + 1, max_depth_level)
        if not f:
            return None, None
        if not any(f[i] for i in range(len(nested))):
            v = values[-1]
            return sqrt(v[0] + _mexpand(v[1]*d)), f
        else:
            p = Mul(*[nested[i] for i in range(len(nested)) if f[i]])
            v = _sqrt_match(p)
            if 1 in f and f.index(1) < len(nested) - 1 and f[len(nested) - 1]:
                v[0] = -v[0]
                v[1] = -v[1]
            if not f[len(nested)]:  # Solution denests with square roots
                vad = _mexpand(v[0] + d)
                if vad <= 0:
                    # return the radicand from the previous invocation.
                    return sqrt(nested[-1]), [0]*len(nested)
                if not(sqrt_depth(vad) <= sqrt_depth(R) + 1 or
                       (vad**2).is_Number):
                    av0[1] = None
                    return None, None

                sqvad = _sqrtdenest1(sqrt(vad), denester=False)
                if not (sqrt_depth(sqvad) <= sqrt_depth(R) + 1):
                    av0[1] = None
                    return None, None
                sqvad1 = radsimp(1/sqvad)
                res = _mexpand(sqvad/sqrt(2) + (v[1]*sqrt(R)*sqvad1/sqrt(2)))
                return res, f
            else:  # Solution requires a fourth root
                s2 = _mexpand(v[1]*R) + d
                if s2 <= 0:
                    return sqrt(nested[-1]), [0]*len(nested)
                FR, s = root(_mexpand(R), 4), sqrt(s2)
                return _mexpand(s/(sqrt(2)*FR) + v[0]*FR/(sqrt(2)*s)), f


def unrad(eq, *syms, **flags):
    """Remove radicals with symbolic arguments and return (eq, cov),
    None or raise an error:

    None is returned if there are no radicals to remove.

    NotImplementedError is raised if there are radicals and they cannot be
    removed or if the relationship between the original symbols and the
    change of variable needed to rewrite the system as a polynomial cannot
    be solved.

    Otherwise the tuple, ``(eq, cov)``, is returned where::

        ``eq``, ``cov``
            ``eq`` is an equation without radicals (in the symbol(s) of
            interest) whose solutions are a superset of the solutions to the
            original expression. ``eq`` might be re-written in terms of a new
            variable; the relationship to the original variables is given by
            ``cov`` which is a list containing ``v`` and ``v**p - b`` where
            ``p`` is the power needed to clear the radical and ``b`` is the
            radical now expressed as a polynomial in the symbols of interest.
            For example, for sqrt(2 - x) the tuple would be
            ``(c, c**2 - 2 + x)``. The solutions of ``eq`` will contain
            solutions to the original equation (if there are any).

    ``syms``
        an iterable of symbols which, if provided, will limit the focus of
        radical removal: only radicals with one or more of the symbols of
        interest will be cleared. All free symbols are used if ``syms`` is not
        set.

    ``flags`` are used internally for communication during recursive calls.
    Two options are also recognized::

        ``take``, when defined, is interpreted as a single-argument function
        that returns True if a given Pow should be handled.

    Radicals can be removed from an expression if::

        *   all bases of the radicals are the same; a change of variables is
            done in this case.
        *   if all radicals appear in one term of the expression
        *   there are only 4 terms with sqrt() factors or there are less than
            four terms having sqrt() factors
        *   there are only two terms with radicals

    Examples
    ========

    >>> unrad(sqrt(x)*cbrt(x) + 2)
    (x**5 - 64, [])
    >>> unrad(sqrt(x) + root(x + 1, 3))
    (x**3 - x**2 - 2*x - 1, [])
    >>> eq = sqrt(x) + root(x, 3) - 2
    >>> unrad(eq)
    (_p**3 + _p**2 - 2, [_p, -x + _p**6])

    """
    from ..polys.rootoftools import RootOf
    from ..solvers import solve

    uflags = {'check': False, 'simplify': False}

    def _cov(p, e):
        if not cov:
            cov[:] = [p, e]
        else:
            raise NotImplementedError

    def _canonical(eq, cov):
        if cov:
            # change symbol to vanilla so no solutions are eliminated
            p, e = cov
            rep = {p: Dummy(p.name)}
            eq = eq.xreplace(rep)
            cov = [p.xreplace(rep), e.xreplace(rep)]

        # remove constants and powers of factors since these don't change
        # the location of the root; XXX should factor or factor_terms be used?
        eq = factor_terms(_mexpand(eq.as_numer_denom()[0], recursive=True))
        if eq.is_Mul:
            args = []
            for f in eq.args:
                if f.is_number:
                    continue
                if f.is_Pow and _take(f, True):
                    args.append(f.base)
                else:
                    args.append(f)
            eq = Mul(*args)  # leave as Mul for more efficient solving

        # make the sign canonical
        free = eq.free_symbols
        if len(free) == 1:
            if eq.coeff(free.pop()**degree(eq)).could_extract_minus_sign():
                eq = -eq
        elif eq.could_extract_minus_sign():
            eq = -eq

        return eq, cov

    def _Q(pow):
        # return leading Rational of denominator of Pow's exponent
        c = pow.as_base_exp()[1].as_coeff_Mul()[0]
        if not c.is_Rational:
            return Integer(1)
        return c.denominator

    # define the _take method that will determine whether a term is of interest
    def _take(d, take_int_pow):
        # return True if coefficient of any factor's exponent's den is not 1
        for pow in Mul.make_args(d):
            if not (pow.is_Symbol or pow.is_Pow):
                continue
            b, e = pow.as_base_exp()
            if not b.has(*syms):
                continue
            if not take_int_pow and _Q(pow) == 1:
                continue
            free = pow.free_symbols
            if free.intersection(syms):
                return True
        return False
    _take = flags.setdefault('_take', _take)

    cov, nwas, rpt = [flags.setdefault(k, v) for k, v in
                      sorted({'cov': [], 'n': None, 'rpt': 0}.items())]

    # preconditioning
    eq = powdenest(factor_terms(eq, radical=True))
    eq, d = eq.as_numer_denom()
    eq = _mexpand(eq, recursive=True)
    if eq.is_number:
        return

    syms = set(syms) or eq.free_symbols
    poly = eq.as_poly()
    gens = [g for g in poly.gens if _take(g, True)]
    if not gens:
        return

    # check for trivial case
    # - already a polynomial in integer powers
    if all(_Q(g) == 1 for g in gens):
        return
    # - an exponent has a symbol of interest (don't handle)
    if any(g.as_base_exp()[1].has(*syms) for g in gens):
        return

    def _rads_bases_lcm(poly):
        # if all the bases are the same or all the radicals are in one
        # term, `lcm` will be the lcm of the denominators of the
        # exponents of the radicals
        lcm = 1
        rads = set()
        bases = set()
        for g in poly.gens:
            if not _take(g, False):
                continue
            q = _Q(g)
            if q != 1:
                rads.add(g)
                lcm = math.lcm(lcm, q)
                bases.add(g.base)
        return rads, bases, lcm
    rads, bases, lcm = _rads_bases_lcm(poly)

    if not rads:
        return

    covsym = Dummy('p', nonnegative=True)

    # only keep in syms symbols that actually appear in radicals;
    # and update gens
    newsyms = set()
    for r in rads:
        newsyms.update(syms & r.free_symbols)
    if newsyms != syms:
        syms = newsyms
        gens = [g for g in gens if g.free_symbols & syms]

    # get terms together that have common generators
    drad = dict(zip(rads, range(len(rads))))
    rterms = {(): []}
    args = Add.make_args(poly.as_expr())
    for t in args:
        if _take(t, False):
            common = set(t.as_poly().gens).intersection(rads)
            key = tuple(sorted(drad[i] for i in common))
        else:
            key = ()
        rterms.setdefault(key, []).append(t)
    others = Add(*rterms.pop(()))
    rterms = [Add(*rterms[k]) for k in rterms]

    # the output will depend on the order terms are processed, so
    # make it canonical quickly
    rterms = list(reversed(list(ordered(rterms))))

    ok = False  # we don't have a solution yet
    depth = sqrt_depth(eq)

    if len(rterms) == 1 and not (rterms[0].is_Add and lcm > 2):
        eq = rterms[0]**lcm - ((-others)**lcm)
        ok = True
    else:
        if len(rterms) == 1 and rterms[0].is_Add:
            rterms = list(rterms[0].args)
        if len(bases) == 1:
            b = bases.pop()
            if len(syms) > 1:
                free = b.free_symbols
                x = {g for g in gens if g.is_Symbol} & free
                if not x:
                    x = free
                x = ordered(x)
            else:
                x = syms
            x = list(x)[0]
            try:
                inv = solve(covsym**lcm - b, x, **uflags)
                if not inv or any(isinstance(s[x], RootOf) for s in inv):
                    raise NotImplementedError
                eq = poly.as_expr().subs({b: covsym**lcm}).subs(inv[0])
                _cov(covsym, covsym**lcm - b)
                return _canonical(eq, cov)
            except NotImplementedError:  # pragma: no cover
                pass

        if len(rterms) == 2:
            if not others:
                eq = rterms[0]**lcm - (-rterms[1])**lcm
                ok = True
            elif not log(lcm, 2).is_Integer:
                # the lcm-is-power-of-two case is handled below
                r0, r1 = rterms
                if flags.get('_reverse', False):
                    r1, r0 = r0, r1
                i0 = _rads_bases_lcm(r0.as_poly())
                i1 = _rads1, _bases1, lcm1 = _rads_bases_lcm(r1.as_poly())
                for reverse in range(2):
                    if reverse:
                        i0, i1 = i1, i0
                        r0, r1 = r1, r0
                    _rads1, _, lcm1 = i1
                    _rads1 = Mul(*_rads1)
                    t1 = _rads1**lcm1
                    c = covsym**lcm1 - t1
                    for x in syms:
                        try:
                            sol = solve(c, x, **uflags)
                            if not sol or any(isinstance(s[x], RootOf)
                                              for s in sol):
                                raise NotImplementedError
                            neweq = r0.subs(sol[0]) + covsym*r1/_rads1 + others
                            tmp = unrad(neweq, covsym)
                            if tmp:
                                eq, newcov = tmp
                                if not newcov:
                                    _cov(covsym, c)
                                else:
                                    raise NotImplementedError
                            else:
                                eq = neweq
                                _cov(covsym, c)
                            ok = True
                            break
                        except NotImplementedError:
                            if reverse:
                                raise NotImplementedError(
                                    'no successful change of variable found')
                            else:
                                pass
                    if ok:
                        break
        elif len(rterms) == 3:
            # two cube roots and another with order less than 5
            # (so an analytical solution can be found) or a base
            # that matches one of the cube root bases
            info = [_rads_bases_lcm(i.as_poly()) for i in rterms]
            RAD = 0
            BASES = 1
            LCM = 2
            if info[0][LCM] != 3:
                info.append(info.pop(0))
                rterms.append(rterms.pop(0))
            elif info[1][LCM] != 3:
                info.append(info.pop(1))
                rterms.append(rterms.pop(1))
            if info[0][LCM] == info[1][LCM] == 3:
                if info[1][BASES] != info[2][BASES]:
                    info[0], info[1] = info[1], info[0]
                    rterms[0], rterms[1] = rterms[1], rterms[0]
                if info[1][BASES] == info[2][BASES]:
                    eq = rterms[0]**3 + (rterms[1] + rterms[2] + others)**3
                    ok = True
                elif info[2][LCM] < 5:
                    # a*root(A, 3) + b*root(B, 3) + others = c
                    a, b, c, d, A, B = [Dummy(i) for i in 'abcdAB']
                    # zz represents the unraded expression into which the
                    # specifics for this case are substituted
                    zz = (c - d)*(A**3*a**9 + 3*A**2*B*a**6*b**3 -
                                  3*A**2*a**6*c**3 + 9*A**2*a**6*c**2*d - 9*A**2*a**6*c*d**2 +
                                  3*A**2*a**6*d**3 + 3*A*B**2*a**3*b**6 + 21*A*B*a**3*b**3*c**3 -
                                  63*A*B*a**3*b**3*c**2*d + 63*A*B*a**3*b**3*c*d**2 -
                                  21*A*B*a**3*b**3*d**3 + 3*A*a**3*c**6 - 18*A*a**3*c**5*d +
                                  45*A*a**3*c**4*d**2 - 60*A*a**3*c**3*d**3 + 45*A*a**3*c**2*d**4 -
                                  18*A*a**3*c*d**5 + 3*A*a**3*d**6 + B**3*b**9 - 3*B**2*b**6*c**3 +
                                  9*B**2*b**6*c**2*d - 9*B**2*b**6*c*d**2 + 3*B**2*b**6*d**3 +
                                  3*B*b**3*c**6 - 18*B*b**3*c**5*d + 45*B*b**3*c**4*d**2 -
                                  60*B*b**3*c**3*d**3 + 45*B*b**3*c**2*d**4 - 18*B*b**3*c*d**5 +
                                  3*B*b**3*d**6 - c**9 + 9*c**8*d - 36*c**7*d**2 + 84*c**6*d**3 -
                                  126*c**5*d**4 + 126*c**4*d**5 - 84*c**3*d**6 + 36*c**2*d**7 -
                                  9*c*d**8 + d**9)

                    def _t(i):
                        b = Mul(*info[i][RAD])
                        return cancel(rterms[i]/b), Mul(*info[i][BASES])

                    aa, AA = _t(0)
                    bb, BB = _t(1)
                    cc = -rterms[2]
                    dd = others
                    eq = zz.xreplace(dict(zip(
                        (a, A, b, B, c, d),
                        (aa, AA, bb, BB, cc, dd))))
                    ok = True
        # handle power-of-2 cases
        if not ok:
            if log(lcm, 2).is_Integer and (not others and
                                           len(rterms) == 4 or len(rterms) < 4):
                def _norm2(a, b):
                    return a**2 + b**2 + 2*a*b

                if len(rterms) == 4:
                    # (r0+r1)**2 - (r2+r3)**2
                    r0, r1, r2, r3 = rterms
                    eq = _norm2(r0, r1) - _norm2(r2, r3)
                    ok = True
                elif len(rterms) == 3:
                    # (r1+r2)**2 - (r0+others)**2
                    r0, r1, r2 = rterms
                    eq = _norm2(r1, r2) - _norm2(r0, others)
                    ok = True
                elif len(rterms) == 2:
                    # r0**2 - (r1+others)**2
                    r0, r1 = rterms
                    eq = r0**2 - _norm2(r1, others)
                    ok = True

    new_depth = sqrt_depth(eq) if ok else depth
    rpt += 1  # XXX how many repeats with others unchanging is enough?
    if not ok or (nwas is not None and len(rterms) == nwas and
                  new_depth is not None and new_depth == depth and rpt > 3):
        raise NotImplementedError('Cannot remove all radicals')

    flags.update({'cov': cov, 'n': len(rterms), 'rpt': rpt})
    neq = unrad(eq, *syms, **flags)
    if neq:
        eq, cov = neq
    eq, cov = _canonical(eq, cov)
    return eq, cov
