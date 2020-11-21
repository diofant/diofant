"""Real and complex elements."""

from mpmath.ctx_mp_python import PythonMPContext, _constant, _mpc, _mpf
from mpmath.libmp import (MPZ_ONE, finf, fnan, fninf, fone, from_float,
                          from_int, from_str, fzero, int_types, mpf_mul,
                          round_nearest, to_rational)
from mpmath.rational import mpq

from .domainelement import DomainElement


class RealElement(_mpf, DomainElement):
    """An element of a real domain."""

    def _set_mpf(self, val):
        self.__mpf__ = val

    _mpf_ = property(lambda self: self.__mpf__, _set_mpf)

    @property
    def parent(self):
        return self.context._parent

    @property
    def numerator(self):
        return self

    @property
    def denominator(self):
        return self.parent.one

    def __reduce__(self):
        return self.parent.__call__, (self._mpf_,)


class ComplexElement(_mpc, DomainElement):
    """An element of a complex domain."""

    def _set_mpc(self, val):
        self.__mpc__ = val

    _mpc_ = property(lambda self: self.__mpc__, _set_mpc)

    @property
    def parent(self):
        return self.context._parent

    @property
    def numerator(self):
        return self

    @property
    def denominator(self):
        return self.parent.one

    def __reduce__(self):
        return self.parent.__call__, (*self._mpc_,)


class MPContext(PythonMPContext):
    """Base class to keep mpmath evaluation context."""

    def __init__(self, prec=53, dps=None, tol=None):
        new = object.__new__

        self._prec_rounding = [prec, round_nearest]

        if dps is None:
            self._set_prec(prec)
        else:
            self._set_dps(dps)

        self.mpf = type('RealElement', (RealElement,), {})
        self.mpc = type('ComplexElement', (ComplexElement,), {})
        self.mpf._ctxdata = [self.mpf, new, self._prec_rounding]
        self.mpc._ctxdata = [self.mpc, new, self._prec_rounding]
        self.mpf.context = self
        self.mpc.context = self
        self.constant = type('constant', (_constant,), {})
        self.constant._ctxdata = [self.mpf, new, self._prec_rounding]
        self.constant.context = self

        self.types = [self.mpf, self.mpc, self.constant]
        self.trap_complex = True
        self.pretty = True

        if tol is None:
            self.tol = self._make_tol()
        elif tol is False:
            self.tol = fzero
        else:
            self.tol = self._convert_tol(tol)

        self.tolerance = self.make_mpf(self.tol)

        if not self.tolerance:
            self.max_denom = 1000000
        else:
            self.max_denom = int(1/self.tolerance)

        self.zero = self.make_mpf(fzero)
        self.one = self.make_mpf(fone)
        self.j = self.make_mpc((fzero, fone))
        self.inf = self.make_mpf(finf)
        self.ninf = self.make_mpf(fninf)
        self.nan = self.make_mpf(fnan)

    def __eq__(self, other):
        return (isinstance(other, MPContext) and self.prec == other.prec and
                self.dps == other.dps and self.tolerance == other.tolerance)

    def _make_tol(self):
        hundred = (0, 25, 2, 5)
        eps = (0, MPZ_ONE, 1 - self.prec, 1)
        return mpf_mul(hundred, eps)

    def _convert_tol(self, tol):
        if isinstance(tol, int_types):
            return from_int(tol)
        if isinstance(tol, float):
            return from_float(tol)
        if hasattr(tol, '_mpf_'):
            return tol._mpf_
        prec, rounding = self._prec_rounding
        if isinstance(tol, str):
            return from_str(tol, prec, rounding)
        raise ValueError(f'expected a real number, got {tol}')

    def _convert_fallback(self, x, strings):
        raise TypeError('cannot create mpf from ' + str(x))

    @property
    def _str_digits(self):
        return self._dps

    def to_rational(self, s, limit=True):
        p, q = to_rational(s._mpf_)

        if not limit or q <= self.max_denom:
            return p, q

        p0, q0, p1, q1 = 0, 1, 1, 0
        n, d = p, q

        while True:
            a = n//d
            q2 = q0 + a*q1
            if q2 > self.max_denom:
                break
            p0, q0, p1, q1 = p1, q1, p0 + a*p1, q2
            n, d = d, n - a*d

        k = (self.max_denom - q0)//q1

        number = mpq(p, q)
        bound1 = mpq(p0 + k*p1, q0 + k*q1)
        bound2 = mpq(p1, q1)

        if not bound2 or not bound1:
            return p, q
        elif abs(bound2 - number) <= abs(bound1 - number):
            return bound2._mpq_
        else:
            return bound1._mpq_

    def almosteq(self, s, t, rel_eps=None, abs_eps=None):
        t = self.convert(t)
        if abs_eps is None and rel_eps is None:
            rel_eps = abs_eps = self.tolerance
        if abs_eps is None:
            abs_eps = self.convert(rel_eps)
        elif rel_eps is None:
            rel_eps = self.convert(abs_eps)
        diff = abs(s - t)
        if diff <= abs_eps:
            return True
        abss = abs(s)
        abst = abs(t)
        if abss < abst:
            err = diff/abst
        else:
            err = diff/abss
        return err <= rel_eps
