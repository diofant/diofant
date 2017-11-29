"""Implementation of RootOf class and related tools. """

from math import log as mathlog

from mpmath import findroot, mpc, mpf, workprec
from mpmath.libmp.libmpf import prec_to_dps

from ..core import (Add, Dummy, Expr, Float, I, Integer, Lambda, Rational, S,
                    cacheit, symbols, sympify)
from ..core.evaluate import global_evaluate
from ..core.function import AppliedUndef
from ..domains import QQ
from ..functions import root as _root
from ..functions import sign
from ..utilities import lambdify, public
from .polyerrors import (DomainError, GeneratorsNeeded,
                         MultivariatePolynomialError, PolynomialError)
from .polyfuncs import symmetrize, viete
from .polyroots import (preprocess_roots, roots, roots_binomial, roots_cubic,
                        roots_linear, roots_quadratic, roots_quartic)
from .polytools import Poly, PurePoly, factor
from .rationaltools import together
from .rootisolation import (dup_isolate_complex_roots_sqf,
                            dup_isolate_real_roots_sqf)


def _ispow2(i):
    v = mathlog(i, 2)
    return v == int(v)


_reals_cache = {}
_complexes_cache = {}


@public
class RootOf(Expr):
    """Represents ``k``-th root of a univariate polynomial. """

    is_commutative = True

    def __new__(cls, f, x, index=None, radicals=True, expand=True,
                evaluate=None):
        """Construct a new ``RootOf`` object for ``k``-th root of ``f``. """
        x = sympify(x)

        if index is None and x.is_Integer:
            x, index = None, x
        else:
            index = sympify(index)

        if index is not None and index.is_Integer:
            index = int(index)
        else:
            raise ValueError("expected an integer root index, got %s" % index)

        poly = PurePoly(f, x, expand=expand)

        if not poly.is_univariate:
            raise PolynomialError("only univariate polynomials are allowed")

        degree = poly.degree()
        dom = poly.domain

        if degree <= 0:
            raise PolynomialError("can't construct RootOf object for %s" % f)

        if index < -degree or index >= degree:
            raise IndexError("root index out of [%d, %d] range, got %d" %
                             (-degree, degree - 1, index))
        elif index < 0:
            index += degree

        if not dom.is_ZZ and poly.LC().is_nonzero is False:
            raise NotImplementedError("sorted roots not supported over %s" % dom)

        if evaluate is None:
            evaluate = global_evaluate[0]

        if not evaluate:
            obj = Expr.__new__(cls)

            obj.poly = poly
            obj.index = index

            return obj

        if not dom.is_Exact:
            poly = poly.to_exact()

        roots = cls._roots_trivial(poly, radicals)

        if roots is not None:
            return roots[index]

        coeff, poly = preprocess_roots(poly)
        dom = poly.domain

        if dom.is_ZZ:
            root = cls._indexed_root(poly, index)
        else:
            root = poly, index

        return coeff*cls._postprocess_root(root, radicals)

    @classmethod
    def _new(cls, poly, index):
        """Construct new ``RootOf`` object from raw data. """
        obj = Expr.__new__(cls)

        obj.poly = PurePoly(poly)
        obj.index = index

        try:
            _reals_cache[obj.poly] = _reals_cache[poly]
            _complexes_cache[obj.poly] = _complexes_cache[poly]
        except KeyError:
            pass

        return obj

    def _hashable_content(self):
        return self.poly, self.index

    @property
    def expr(self):
        return self.poly.as_expr()

    @property
    def args(self):
        return self.expr, self.poly.gen, Integer(self.index)

    @property
    def free_symbols(self):
        return self.poly.free_symbols

    def _eval_is_real(self):
        try:
            return self.index < len(_reals_cache[self.poly])
        except KeyError:
            pass
    _eval_is_extended_real = _eval_is_real

    def _eval_is_complex(self):
        if all(_.is_complex for _ in self.poly.coeffs()):
            return True

    def _eval_is_algebraic(self):
        if all(_.is_algebraic for _ in self.poly.coeffs()):
            return True

    def _eval_power(self, expt):
        p = self.poly
        if p.degree() == expt and p.length() == 2 and p.TC():
            return -p.TC()/p.LC()

    def _eval_rewrite_as_Pow(self, e, x, i):
        p = self.poly
        n = p.degree()
        if n == 3:
            return roots_cubic(p)[i]
        elif n == 4:
            return roots_quartic(p)[i]

    @property
    def is_number(self):
        return not self.free_symbols

    @classmethod
    def real_roots(cls, poly, radicals=True):
        """Get real roots of a polynomial. """
        return cls._get_roots("_real_roots", poly, radicals)

    @classmethod
    def all_roots(cls, poly, radicals=True):
        """Get real and complex roots of a polynomial. """
        return cls._get_roots("_all_roots", poly, radicals)

    @classmethod
    def _get_reals_sqf(cls, factor):
        """Compute real root isolating intervals for a square-free polynomial. """
        if factor in _reals_cache:
            real_part = _reals_cache[factor]
        else:
            _reals_cache[factor] = real_part = \
                dup_isolate_real_roots_sqf(
                    factor.rep.rep, factor.rep.domain, blackbox=True)

        return real_part

    @classmethod
    def _get_complexes_sqf(cls, factor):
        """Compute complex root isolating intervals for a square-free polynomial. """
        if factor in _complexes_cache:
            complex_part = _complexes_cache[factor]
        else:
            _complexes_cache[factor] = complex_part = \
                dup_isolate_complex_roots_sqf(
                factor.rep.rep, factor.rep.domain, blackbox=True)
        return complex_part

    @classmethod
    def _get_reals(cls, factors):
        """Compute real root isolating intervals for a list of factors. """
        reals = []

        for factor, k in factors:
            real_part = cls._get_reals_sqf(factor)
            reals.extend([ (root, factor, k) for root in real_part ])

        return reals

    @classmethod
    def _get_complexes(cls, factors):
        """Compute complex root isolating intervals for a list of factors. """
        complexes = []

        for factor, k in factors:
            complex_part = cls._get_complexes_sqf(factor)
            complexes.extend([ (root, factor, k) for root in complex_part ])

        return complexes

    @classmethod
    def _reals_sorted(cls, reals):
        """Make real isolating intervals disjoint and sort roots. """
        cache = {}

        for i, (u, f, k) in enumerate(reals):
            for j, (v, g, m) in enumerate(reals[i + 1:]):
                u, v = u.refine_disjoint(v)
                reals[i + j + 1] = (v, g, m)

            reals[i] = (u, f, k)

        reals = sorted(reals, key=lambda r: r[0].a)

        for root, factor, _ in reals:
            if factor in cache:
                cache[factor].append(root)
            else:
                cache[factor] = [root]

        for factor, roots in cache.items():
            _reals_cache[factor] = roots

        return reals

    @classmethod
    def _separate_imaginary_from_complex(cls, complexes):
        from ..utilities import sift

        def is_imag(c):
            """
            return True if all roots are imaginary (ax**2 + b)
            return False if no roots are imaginary
            return None if 2 roots are imaginary (ax**N
            """
            u, f, k = c
            deg = f.degree()
            if f.length() == 2:
                if deg == 2:
                    return True  # both imag
                elif _ispow2(deg):
                    if f.LC()*f.TC() < 0:
                        return  # 2 are imag
            return False  # none are imag

        # separate according to the function
        sifted = sift(complexes, lambda c: c[1])
        del complexes
        imag = []
        complexes = []
        for f in sifted:
            isift = sift(sifted[f], lambda c: is_imag(c))
            imag.extend(isift.pop(True, []))
            complexes.extend(isift.pop(False, []))
            mixed = isift.pop(None, [])
            assert not isift
            if not mixed:
                continue
            while True:
                # the non-imaginary ones will be on one side or the other
                # of the y-axis
                i = 0
                while i < len(mixed):
                    u, f, k = mixed[i]
                    if u.ax*u.bx > 0:
                        complexes.append(mixed.pop(i))
                    else:
                        i += 1
                if len(mixed) == 2:
                    imag.extend(mixed)
                    break
                # refine
                for i, (u, f, k) in enumerate(mixed):
                    u = u._inner_refine()
                    mixed[i] = u, f, k
        return imag, complexes

    @classmethod
    def _refine_complexes(cls, complexes):
        """return complexes such that no bounding rectangles of non-conjugate
        roots would intersect if slid horizontally or vertically/
        """
        while complexes:  # break when all are distinct
            # get the intervals pairwise-disjoint. If rectangles were drawn around
            # the coordinates of the bounding rectangles, no rectangles would
            # intersect after this procedure
            for i, (u, f, k) in enumerate(complexes):
                for j, (v, g, m) in enumerate(complexes[i + 1:]):
                    u, v = u.refine_disjoint(v)
                    complexes[i + j + 1] = (v, g, m)

                complexes[i] = (u, f, k)
            # Although there are no intersecting rectangles, a given rectangle
            # might intersect another when slid horizontally. We have to refine
            # intervals until this is not true so we can sort the roots
            # unambiguously. Since complex roots come in conjugate pairs, we
            # will always have 2 rectangles above each other but we should not
            # have more than that.
            N = len(complexes)//2 - 1
            # check x (real) parts: there must be N + 1 disjoint x ranges, i.e.
            # the first one must be different from N others
            uu = {(u.ax, u.bx) for u, _, _ in complexes}
            u = uu.pop()
            if sum(u[1] <= v[0] or v[1] <= u[0] for v in uu) < N:
                # refine
                for i, (u, f, k) in enumerate(complexes):
                    u = u._inner_refine()
                    complexes[i] = u, f, k
            else:
                # intervals with identical x-values have disjoint y-values or
                # else they would not be disjoint so there is no need for
                # further checks
                break
        return complexes

    @classmethod
    def _complexes_sorted(cls, complexes):
        """Make complex isolating intervals disjoint and sort roots. """
        if not complexes:
            return []
        cache = {}

        # imaginary roots can cause a problem in terms of sorting since
        # their x-intervals will never refine as distinct from others
        # so we handle them separately
        imag, complexes = cls._separate_imaginary_from_complex(complexes)
        complexes = cls._refine_complexes(complexes)

        # sort imaginary roots
        def key(c):
            """return, for ax**n+b, +/-root(abs(b/a), b) according to the
            apparent sign of the imaginary interval, e.g. if the interval
            were (0, 3) the positive root would be returned.
            """
            u, f, k = c
            r = _root(abs(f.TC()/f.LC()), f.degree())
            if u.ay < 0 or u.by < 0:
                return -r
            return r
        imag = sorted(imag, key=lambda c: key(c))

        # sort complexes and combine with imag
        if complexes:
            # key is (x1, y1) e.g. (1, 2)x(3, 4) -> (1,3)
            complexes = sorted(complexes, key=lambda c: c[0].a)
            # find insertion point for imaginary
            for i, c in enumerate(reversed(complexes)):
                if c[0].bx <= 0:
                    break
            i = len(complexes) - i - 1
            if i:
                i += 1
            complexes = complexes[:i] + imag + complexes[i:]
        else:
            complexes = imag

        # update cache
        for root, factor, _ in complexes:
            if factor in cache:
                cache[factor].append(root)
            else:
                cache[factor] = [root]

        for factor, roots in cache.items():
            _complexes_cache[factor] = roots

        return complexes

    @classmethod
    def _reals_index(cls, reals, index):
        """Map initial real root index to an index in a factor where the root belongs. """
        i = 0

        for j, (_, factor, k) in enumerate(reals):  # pragma: no branch
            if index < i + k:
                poly, index = factor, 0

                for _, factor, _ in reals[:j]:
                    if factor == poly:
                        index += 1

                return poly, index
            else:
                i += k

    @classmethod
    def _complexes_index(cls, complexes, index):
        """Map initial complex root index to an index in a factor where the root belongs. """
        index, i = index, 0

        for j, (_, factor, k) in enumerate(complexes):  # pragma: no branch
            if index < i + k:
                poly, index = factor, 0

                for _, factor, _ in complexes[:j]:
                    if factor == poly:
                        index += 1

                index += len(_reals_cache[poly])

                return poly, index
            else:
                i += k

    @classmethod
    def _count_roots(cls, roots):
        """Count the number of real or complex roots including multiplicities."""
        return sum(k for _, _, k in roots)

    @classmethod
    def _indexed_root(cls, poly, index):
        """Get a root of a composite polynomial by index. """
        (_, factors) = poly.factor_list()

        reals = cls._get_reals(factors)
        reals_count = cls._count_roots(reals)

        if index < reals_count:
            reals = cls._reals_sorted(reals)
            return cls._reals_index(reals, index)
        else:
            complexes = cls._get_complexes(factors)
            complexes = cls._complexes_sorted(complexes)
            return cls._complexes_index(complexes, index - reals_count)

    @classmethod
    def _real_roots(cls, poly):
        """Get real roots of a composite polynomial. """
        (_, factors) = poly.factor_list()

        reals = cls._get_reals(factors)
        reals = cls._reals_sorted(reals)
        reals_count = cls._count_roots(reals)

        roots = []

        for index in range(reals_count):
            roots.append(cls._reals_index(reals, index))

        return roots

    @classmethod
    def _all_roots(cls, poly):
        """Get real and complex roots of a composite polynomial. """

        if not poly.domain.is_ZZ:
            return [(poly, i) for i in range(poly.degree())]

        (_, factors) = poly.factor_list()

        reals = cls._get_reals(factors)
        reals = cls._reals_sorted(reals)
        reals_count = cls._count_roots(reals)

        roots = []

        for index in range(reals_count):
            roots.append(cls._reals_index(reals, index))

        complexes = cls._get_complexes(factors)
        complexes = cls._complexes_sorted(complexes)
        complexes_count = cls._count_roots(complexes)

        for index in range(complexes_count):
            roots.append(cls._complexes_index(complexes, index))

        return roots

    @classmethod
    @cacheit
    def _roots_trivial(cls, poly, radicals):
        """Compute roots in linear, quadratic and binomial cases. """
        n = poly.degree()

        if n == 1:
            return roots_linear(poly)

        if not radicals:
            return

        if n == 2:
            return roots_quadratic(poly)
        elif poly.length() == 2 and poly.coeff_monomial(1):
            if not poly.free_symbols_in_domain:
                return roots_binomial(poly)
            elif all(sign(_) in (-1, 1) for _ in poly.coeffs()):
                lc, tc = poly.LC(), poly.TC()
                x, r = poly.gen, _root(abs(tc/lc), n)
                poly = Poly(x**n + sign(lc*tc), x)
                return [r*_ for _ in cls._roots_trivial(poly, radicals)]

    @classmethod
    def _preprocess_roots(cls, poly):
        """Take heroic measures to make ``poly`` compatible with ``RootOf``. """
        dom = poly.domain

        if not dom.is_Exact:
            poly = poly.to_exact()

        coeff, poly = preprocess_roots(poly)
        dom = poly.domain

        if not dom.is_ZZ and poly.LC().is_nonzero is False:
            raise NotImplementedError("sorted roots not supported over %s" % dom)

        return coeff, poly

    @classmethod
    def _postprocess_root(cls, root, radicals):
        """Return the root if it is trivial or a ``RootOf`` object. """
        poly, index = root
        roots = cls._roots_trivial(poly, radicals)

        if roots is not None:
            return roots[index]
        else:
            return cls._new(poly, index)

    @classmethod
    def _get_roots(cls, method, poly, radicals):
        """Return postprocessed roots of specified kind. """

        coeff, poly = cls._preprocess_roots(poly)
        roots = []

        for root in getattr(cls, method)(poly):
            roots.append(coeff*cls._postprocess_root(root, radicals))

        return roots

    @property
    def interval(self):
        """Internal function for retrieving isolation interval from cache. """
        if self.is_real:
            return _reals_cache[self.poly][self.index]
        else:
            reals_count = len(_reals_cache[self.poly])
            return _complexes_cache[self.poly][self.index - reals_count]

    def _eval_subs(self, old, new):
        if old in self.free_symbols:
            return self.func(self.poly.subs(old, new), *self.args[1:])
        else:
            # don't allow subs to change anything
            return self

    def _eval_evalf(self, prec):
        """Evaluate this complex root to the given precision. """
        with workprec(prec):
            g = self.poly.gen
            if not g.is_Symbol:
                d = Dummy('x')
                func = lambdify(d, self.expr.subs(g, d))
            else:
                func = lambdify(g, self.expr)

            try:
                interval = self.interval
            except KeyError:
                return super(Expr, self)._eval_evalf(prec)

            while True:
                if self.is_extended_real:
                    a = mpf(str(interval.a))
                    b = mpf(str(interval.b))
                    if a == b:
                        root = a
                        break
                    x0 = mpf(str(interval.center))
                else:
                    ax = mpf(str(interval.ax))
                    bx = mpf(str(interval.bx))
                    ay = mpf(str(interval.ay))
                    by = mpf(str(interval.by))
                    x0 = mpc(*map(str, interval.center))
                    if ax == bx and ay == by:
                        root = x0
                        break

                try:
                    root = findroot(func, x0)
                    # If the (real or complex) root is not in the 'interval',
                    # then keep refining the interval. This happens if findroot
                    # accidentally finds a different root outside of this
                    # interval because our initial estimate 'x0' was not close
                    # enough. It is also possible that the secant method will
                    # get trapped by a max/min in the interval; the root
                    # verification by findroot will raise a ValueError in this
                    # case and the interval will then be tightened -- and
                    # eventually the root will be found.
                    if self.is_extended_real:
                        if (a <= root <= b):
                            break
                    elif (ax <= root.real <= bx and ay <= root.imag <= by
                          and (interval.ay > 0 or interval.by < 0)):
                        break
                except ValueError:
                    pass
                interval = interval.refine()

        return (Float._new(root.real._mpf_, prec) +
                I*Float._new(root.imag._mpf_, prec))

    def eval_rational(self, tol):
        """
        Returns a Rational approximation to ``self`` with the tolerance ``tol``.

        The returned instance will be at most 'tol' from the exact root.

        The following example first obtains Rational approximation to 1e-7
        accuracy for all roots of the 4-th order Legendre polynomial, and then
        evaluates it to 5 decimal digits (so all digits will be correct
        including rounding):

        >>> from diofant import Rational, legendre_poly, Symbol
        >>> x = Symbol("x")
        >>> p = legendre_poly(4, x, polys=True)
        >>> roots = [r.eval_rational(Rational(1, 10)**7) for r in p.real_roots()]
        >>> roots = [str(r.n(5)) for r in roots]
        >>> roots
        ['-0.86114', '-0.33998', '0.33998', '0.86114']

        """

        if not self.is_extended_real:
            raise NotImplementedError("eval_rational() only works for real polynomials so far")
        interval = self.interval
        while interval.b - interval.a > tol:
            interval = interval.refine()
        a = Rational(str(interval.a))
        b = Rational(str(interval.b))
        return (a + b)/2

    def _eval_Eq(self, other):
        # RootOf represents a Root, so if other is that root, it should set
        # the expression to zero *and* it should be in the interval of the
        # RootOf instance. It must also be a number that agrees with the
        # is_real value of the RootOf instance.
        if type(self) == type(other):
            return sympify(self.__eq__(other))
        if not (other.is_number and not other.has(AppliedUndef)):
            return S.false
        if not other.is_finite:
            return S.false
        z = self.expr.subs(self.expr.free_symbols.pop(), other).is_zero
        if z is False:  # all roots will make z True but we don't know
                        # whether this is the right root if z is True
            return S.false
        o = other.is_extended_real, other.is_imaginary
        s = self.is_extended_real, self.is_imaginary
        if o != s and None not in o and None not in s:
            return S.false
        i = self.interval
        was = i.a, i.b
        need = [True]*2
        # make sure it would be distinct from others
        while any(need):
            i = i.refine()
            a, b = i.a, i.b
            if need[0] and a != was[0]:
                need[0] = False
            if need[1] and b != was[1]:
                need[1] = False
        re, im = other.as_real_imag()
        if not im:
            if self.is_extended_real:
                a, b = [Rational(str(i)) for i in (a, b)]
                return sympify(a < other and other < b)
            return S.false
        if self.is_extended_real:
            return S.false
        z = r1, r2, i1, i2 = [Rational(str(j)) for j in (
            i.ax, i.bx, i.ay, i.by)]
        return sympify((
            r1 < re and re < r2) and (
            i1 < im and im < i2))


@public
class RootSum(Expr):
    """Represents a sum of all roots of a univariate polynomial. """

    def __new__(cls, expr, func=None, x=None, auto=True, quadratic=False):
        """Construct a new ``RootSum`` instance carrying all roots of a polynomial. """
        coeff, poly = cls._transform(expr, x)

        if not poly.is_univariate:
            raise MultivariatePolynomialError(
                "only univariate polynomials are allowed")

        if func is None:
            func = Lambda(poly.gen, poly.gen)
        else:
            try:
                is_func = func.is_Function
            except AttributeError:
                is_func = False

            if is_func and 1 in func.nargs:
                if not isinstance(func, Lambda):
                    func = Lambda(poly.gen, func(poly.gen))
            else:
                raise ValueError(
                    "expected a univariate function, got %s" % func)

        var, expr = func.variables[0], func.expr

        if coeff is not S.One:
            expr = expr.subs(var, coeff*var)

        deg = poly.degree()

        if not expr.has(var):
            return deg*expr

        if expr.is_Add:
            add_const, expr = expr.as_independent(var)
        else:
            add_const = S.Zero

        if expr.is_Mul:
            mul_const, expr = expr.as_independent(var)
        else:
            mul_const = S.One

        func = Lambda(var, expr)

        rational = cls._is_func_rational(poly, func)
        factors, terms = poly.factor_list()[1], []

        for poly, k in factors:
            if poly.is_linear:
                term = func(roots_linear(poly)[0])
            elif quadratic and poly.is_quadratic:
                term = sum(map(func, roots_quadratic(poly)))
            else:
                if not rational or not auto:
                    term = cls._new(poly, func, auto)
                else:
                    term = cls._rational_case(poly, func)

            terms.append(k*term)

        return mul_const*Add(*terms) + deg*add_const

    @classmethod
    def _new(cls, poly, func, auto=True):
        """Construct new raw ``RootSum`` instance. """
        obj = Expr.__new__(cls)

        obj.poly = poly
        obj.fun = func
        obj.auto = auto

        return obj

    @classmethod
    def new(cls, poly, func, auto=True):
        """Construct new ``RootSum`` instance. """

        rational = cls._is_func_rational(poly, func)

        if not rational or not auto:
            return cls._new(poly, func, auto)
        else:
            return cls._rational_case(poly, func)

    @classmethod
    def _transform(cls, expr, x):
        """Transform an expression to a polynomial. """
        poly = PurePoly(expr, x, greedy=False)
        return preprocess_roots(poly)

    @classmethod
    def _is_func_rational(cls, poly, func):
        """Check if a lambda is areational function. """
        var, expr = func.variables[0], func.expr
        return expr.is_rational_function(var)

    @classmethod
    def _rational_case(cls, poly, func):
        """Handle the rational function case. """
        roots = symbols('r:%d' % poly.degree())
        var, expr = func.variables[0], func.expr

        f = sum(expr.subs(var, r) for r in roots)
        p, q = together(f).as_numer_denom()

        domain = QQ[roots]

        p = p.expand()
        q = q.expand()

        try:
            p = Poly(p, domain=domain, expand=False)
        except GeneratorsNeeded:
            p, p_coeff = None, (p,)
        else:
            p_monom, p_coeff = zip(*p.terms())

        try:
            q = Poly(q, domain=domain, expand=False)
        except GeneratorsNeeded:
            q, q_coeff = None, (q,)
        else:
            q_monom, q_coeff = zip(*q.terms())

        coeffs, mapping = symmetrize(p_coeff + q_coeff, formal=True)
        formulas, values = viete(poly, roots), []

        for (sym, _), (_, val) in zip(mapping, formulas):
            values.append((sym, val))

        for i, (coeff, _) in enumerate(coeffs):
            coeffs[i] = coeff.subs(values)

        n = len(p_coeff)

        p_coeff = coeffs[:n]
        q_coeff = coeffs[n:]

        if p is not None:
            p = Poly(dict(zip(p_monom, p_coeff)), *p.gens).as_expr()
        else:
            (p,) = p_coeff

        if q is not None:
            q = Poly(dict(zip(q_monom, q_coeff)), *q.gens).as_expr()
        else:
            (q,) = q_coeff

        return factor(p/q)

    def _hashable_content(self):
        return self.poly, self.fun

    @property
    def expr(self):
        return self.poly.as_expr()

    @property
    def args(self):
        return self.expr, self.fun, self.poly.gen

    @property
    def free_symbols(self):
        return self.poly.free_symbols | self.fun.free_symbols

    @property
    def is_commutative(self):
        return True

    def doit(self, **hints):
        _roots = roots(self.poly, multiple=True)

        if len(_roots) < self.poly.degree():
            return self
        else:
            return Add(*[ self.fun(r) for r in _roots ])

    def _eval_evalf(self, prec):
        try:
            _roots = self.poly.nroots(n=prec_to_dps(prec))
        except (DomainError, PolynomialError):
            return self
        else:
            return Add(*[ self.fun(r) for r in _roots ])

    def _eval_derivative(self, x):
        var, expr = self.fun.args
        func = Lambda(var, expr.diff(x))
        return self.new(self.poly, func, self.auto)
