"""Elliptic integrals."""

from ...core import Function, I, Integer, Rational, oo, pi, zoo
from ...core.function import ArgumentIndexError
from ..elementary.complexes import sign
from ..elementary.hyperbolic import atanh
from ..elementary.miscellaneous import sqrt
from ..elementary.trigonometric import sin, tan
from .gamma_functions import gamma
from .hyper import hyper, meijerg


class elliptic_k(Function):
    r"""
    The complete elliptic integral of the first kind, defined by

    .. math:: K(m) = F\left(\tfrac{\pi}{2}\middle| m\right)

    where `F\left(z\middle| m\right)` is the Legendre incomplete
    elliptic integral of the first kind.

    The function `K(m)` is a single-valued function on the complex
    plane with branch cut along the interval `(1, \infty)`.

    Examples
    ========

    >>> elliptic_k(0)
    pi/2
    >>> elliptic_k(1.0 + I)
    1.50923695405127 + 0.625146415202697*I
    >>> elliptic_k(m).series(m, n=3)
    pi/2 + pi*m/8 + 9*pi*m**2/128 + O(m**3)

    References
    ==========

    * https://en.wikipedia.org/wiki/Elliptic_integrals
    * http://functions.wolfram.com/EllipticIntegrals/EllipticK

    See Also
    ========

    elliptic_f

    """

    @classmethod
    def eval(cls, m):
        if m == 0:
            return pi/2
        elif m == Rational(1, 2):
            return 8*pi**Rational(3, 2)/gamma(-Rational(1, 4))**2
        elif m == 1:
            return zoo
        elif m == -1:
            return gamma(Rational(1, 4))**2/(4*sqrt(2*pi))
        elif m in (oo, -oo, I*oo, I*-oo, zoo):
            return Integer(0)

    def fdiff(self, argindex=1):
        m = self.args[0]
        if argindex == 1:
            return (elliptic_e(m) - (1 - m)*elliptic_k(m))/(2*m*(1 - m))
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        m = self.args[0]
        if (m.is_extended_real and (m - 1).is_positive) is False:
            return self.func(m.conjugate())

    def _eval_nseries(self, x, n, logx):
        from ...simplify import hyperexpand
        return hyperexpand(self.rewrite(hyper)._eval_nseries(x, n=n, logx=logx))

    def _eval_rewrite_as_hyper(self, m):
        return (pi/2)*hyper((Rational(1, 2), Rational(1, 2)), (1,), m)

    def _eval_rewrite_as_meijerg(self, m):
        return meijerg(((Rational(1, 2), Rational(1, 2)), []), ((0,), (0,)), -m)/2


class elliptic_f(Function):
    r"""
    The Legendre incomplete elliptic integral of the first
    kind, defined by

    .. math:: F\left(z\middle| m\right) =
              \int_0^z \frac{dt}{\sqrt{1 - m \sin^2 t}}

    This function reduces to a complete elliptic integral of
    the first kind, `K(m)`, when `z = \pi/2`.

    Examples
    ========

    >>> elliptic_f(z, m).series(z)
    z + z**5*(3*m**2/40 - m/30) + m*z**3/6 + O(z**6)
    >>> elliptic_f(3.0 + I/2, 1.0 + I)
    2.909449841483 + 1.74720545502474*I

    References
    ==========

    * https://en.wikipedia.org/wiki/Elliptic_integrals
    * http://functions.wolfram.com/EllipticIntegrals/EllipticF

    See Also
    ========

    elliptic_k

    """

    @classmethod
    def eval(cls, z, m):
        k = 2*z/pi
        if m.is_zero:
            return z
        elif z.is_zero:
            return Integer(0)
        elif k.is_integer:
            return k*elliptic_k(m)
        elif m in (oo, -oo):
            return Integer(0)
        elif z.could_extract_minus_sign():
            return -elliptic_f(-z, m)

    def fdiff(self, argindex=1):
        z, m = self.args
        fm = sqrt(1 - m*sin(z)**2)
        if argindex == 1:
            return 1/fm
        elif argindex == 2:
            return (elliptic_e(z, m)/(2*m*(1 - m)) - elliptic_f(z, m)/(2*m) -
                    sin(2*z)/(4*(1 - m)*fm))
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        z, m = self.args
        if (m.is_extended_real and (m - 1).is_positive) is False:
            return self.func(z.conjugate(), m.conjugate())


class elliptic_e(Function):
    r"""
    Called with two arguments `z` and `m`, evaluates the
    incomplete elliptic integral of the second kind, defined by

    .. math:: E\left(z\middle| m\right) = \int_0^z \sqrt{1 - m \sin^2 t} dt

    Called with a single argument `m`, evaluates the Legendre complete
    elliptic integral of the second kind

    .. math:: E(m) = E\left(\tfrac{\pi}{2}\middle| m\right)

    The function `E(m)` is a single-valued function on the complex
    plane with branch cut along the interval `(1, \infty)`.

    Examples
    ========

    >>> elliptic_e(z, m).series(z)
    z + z**5*(-m**2/40 + m/30) - m*z**3/6 + O(z**6)
    >>> elliptic_e(m).series(m, n=4)
    pi/2 - pi*m/8 - 3*pi*m**2/128 - 5*pi*m**3/512 + O(m**4)
    >>> elliptic_e(1 + I, 2 - I/2).evalf()
    1.55203744279187 + 0.290764986058437*I
    >>> elliptic_e(0)
    pi/2
    >>> elliptic_e(2.0 - I)
    0.991052601328069 + 0.81879421395609*I

    References
    ==========

    * https://en.wikipedia.org/wiki/Elliptic_integrals
    * http://functions.wolfram.com/EllipticIntegrals/EllipticE2
    * http://functions.wolfram.com/EllipticIntegrals/EllipticE

    """

    @classmethod
    def eval(cls, z, m=None):
        if m is not None:
            k = 2*z/pi
            if m.is_zero:
                return z
            if z.is_zero:
                return Integer(0)
            elif k.is_integer:
                return k*elliptic_e(m)
            elif m in (oo, -oo):
                return zoo
            elif z.could_extract_minus_sign():
                return -elliptic_e(-z, m)
        else:
            m = z
            if m.is_zero:
                return pi/2
            elif m == 1:
                return Integer(1)
            elif m is oo:
                return I*oo
            elif m == -oo:
                return oo
            elif m is zoo:
                return zoo

    def fdiff(self, argindex=1):
        if len(self.args) == 2:
            z, m = self.args
            if argindex == 1:
                return sqrt(1 - m*sin(z)**2)
            elif argindex == 2:
                return (elliptic_e(z, m) - elliptic_f(z, m))/(2*m)
            else:
                raise ArgumentIndexError(self, argindex)
        else:
            m = self.args[0]
            if argindex == 1:
                return (elliptic_e(m) - elliptic_k(m))/(2*m)
            else:
                raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        if len(self.args) == 2:
            z, m = self.args
            if (m.is_extended_real and (m - 1).is_positive) is False:
                return self.func(z.conjugate(), m.conjugate())
        else:
            m = self.args[0]
            if (m.is_extended_real and (m - 1).is_positive) is False:
                return self.func(m.conjugate())

    def _eval_nseries(self, x, n, logx):
        from ...simplify import hyperexpand
        if len(self.args) == 1:
            return hyperexpand(self.rewrite(hyper)._eval_nseries(x, n=n, logx=logx))
        return super()._eval_nseries(x, n=n, logx=logx)

    def _eval_rewrite_as_hyper(self, *args):
        if len(args) == 1:
            m = args[0]
            return (pi/2)*hyper((-Rational(1, 2), Rational(1, 2)), (1,), m)

    def _eval_rewrite_as_meijerg(self, *args):
        if len(args) == 1:
            m = args[0]
            return -meijerg(((Rational(1, 2), Rational(3, 2)), []),
                            ((0,), (0,)), -m)/4


class elliptic_pi(Function):
    r"""
    Called with three arguments `n`, `z` and `m`, evaluates the
    Legendre incomplete elliptic integral of the third kind, defined by

    .. math:: \Pi\left(n; z\middle| m\right) = \int_0^z \frac{dt}
              {\left(1 - n \sin^2 t\right) \sqrt{1 - m \sin^2 t}}

    Called with two arguments `n` and `m`, evaluates the complete
    elliptic integral of the third kind:

    .. math:: \Pi\left(n\middle| m\right) =
              \Pi\left(n; \tfrac{\pi}{2}\middle| m\right)

    Examples
    ========

    >>> elliptic_pi(n, z, m).series(z, n=4)
    z + z**3*(m/6 + n/3) + O(z**4)
    >>> elliptic_pi(0.5 + I, 1.0 - I, 1.2)
    2.50232379629182 - 0.760939574180767*I
    >>> elliptic_pi(0, 0)
    pi/2
    >>> elliptic_pi(1.0 - I/3, 2.0 + I)
    3.29136443417283 + 0.32555634906645*I

    References
    ==========

    * https://en.wikipedia.org/wiki/Elliptic_integrals
    * http://functions.wolfram.com/EllipticIntegrals/EllipticPi3
    * http://functions.wolfram.com/EllipticIntegrals/EllipticPi

    """

    @classmethod
    def eval(cls, n, m, z=None):
        if z is not None:
            n, z, m = n, m, z
            k = 2*z/pi
            if n == 0:
                return elliptic_f(z, m)
            elif n == 1:
                return (elliptic_f(z, m) +
                        (sqrt(1 - m*sin(z)**2)*tan(z) -
                         elliptic_e(z, m))/(1 - m))
            elif k.is_integer:
                return k*elliptic_pi(n, m)
            elif m == 0:
                return atanh(sqrt(n - 1)*tan(z))/sqrt(n - 1)
            elif n == m:
                return (elliptic_f(z, n) - elliptic_pi(1, z, n) +
                        tan(z)/sqrt(1 - n*sin(z)**2))
            elif n in (oo, -oo):
                return Integer(0)
            elif m in (oo, -oo):
                return Integer(0)
            elif z.could_extract_minus_sign():
                return -elliptic_pi(n, -z, m)
        else:
            if n == 0:
                return elliptic_k(m)
            elif n == 1:
                return zoo
            elif m == 0:
                return pi/(2*sqrt(1 - n))
            elif m == 1:
                return -oo/sign(n - 1)
            elif n == m:
                return elliptic_e(n)/(1 - n)
            elif n in (oo, -oo):
                return Integer(0)
            elif m in (oo, -oo):
                return Integer(0)

    def _eval_conjugate(self):
        if len(self.args) == 3:
            n, z, m = self.args
            if (n.is_extended_real and (n - 1).is_positive) is False and \
               (m.is_extended_real and (m - 1).is_positive) is False:
                return self.func(n.conjugate(), z.conjugate(), m.conjugate())
        else:
            n, m = self.args
            return self.func(n.conjugate(), m.conjugate())

    def fdiff(self, argindex=1):
        if len(self.args) == 3:
            n, z, m = self.args
            fm, fn = sqrt(1 - m*sin(z)**2), 1 - n*sin(z)**2
            if argindex == 1:
                return (elliptic_e(z, m) + (m - n)*elliptic_f(z, m)/n +
                        (n**2 - m)*elliptic_pi(n, z, m)/n -
                        n*fm*sin(2*z)/(2*fn))/(2*(m - n)*(n - 1))
            elif argindex == 2:
                return 1/(fm*fn)
            elif argindex == 3:
                return (elliptic_e(z, m)/(m - 1) +
                        elliptic_pi(n, z, m) -
                        m*sin(2*z)/(2*(m - 1)*fm))/(2*(n - m))
            else:
                raise ArgumentIndexError(self, argindex)
        else:
            n, m = self.args
            if argindex == 1:
                return (elliptic_e(m) + (m - n)*elliptic_k(m)/n +
                        (n**2 - m)*elliptic_pi(n, m)/n)/(2*(m - n)*(n - 1))
            elif argindex == 2:
                return (elliptic_e(m)/(m - 1) + elliptic_pi(n, m))/(2*(n - m))
            else:
                raise ArgumentIndexError(self, argindex)
