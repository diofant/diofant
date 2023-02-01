from ...core import Dummy, Function, I, pi
from ...core.function import ArgumentIndexError
from ..combinatorial.factorials import factorial
from ..elementary.complexes import Abs
from ..elementary.exponential import exp
from ..elementary.miscellaneous import sqrt
from ..elementary.trigonometric import cos, cot, sin
from .polynomials import assoc_legendre


_x = Dummy('dummy_for_spherical_harmonics')


class Ynm(Function):
    r"""
    Spherical harmonics defined as

    .. math::
        Y_n^m(\theta, \varphi) := \sqrt{\frac{(2n+1)(n-m)!}{4\pi(n+m)!}}
                                  \exp(i m \varphi)
                                  \mathrm{P}_n^m\left(\cos(\theta)\right)

    Ynm() gives the spherical harmonic function of order `n` and `m`
    in `\theta` and `\varphi`, `Y_n^m(\theta, \varphi)`. The four
    parameters are as follows: `n \geq 0` an integer and `m` an integer
    such that `-n \leq m \leq n` holds. The two angles are real-valued
    with `\theta \in [0, \pi]` and `\varphi \in [0, 2\pi]`.

    Examples
    ========

    >>> theta = Symbol('theta')
    >>> phi = Symbol('phi')

    >>> Ynm(n, m, theta, phi)
    Ynm(n, m, theta, phi)

    Several symmetries are known, for the order

    >>> theta = Symbol('theta')
    >>> phi = Symbol('phi')

    >>> Ynm(n, -m, theta, phi)
    (-1)**m*E**(-2*I*m*phi)*Ynm(n, m, theta, phi)

    as well as for the angles

    >>> theta = Symbol('theta')
    >>> phi = Symbol('phi')

    >>> Ynm(n, m, -theta, phi)
    Ynm(n, m, theta, phi)

    >>> Ynm(n, m, theta, -phi)
    E**(-2*I*m*phi)*Ynm(n, m, theta, phi)

    For specific integers n and m we can evaluate the harmonics
    to more useful expressions

    >>> simplify(Ynm(0, 0, theta, phi).expand(func=True))
    1/(2*sqrt(pi))

    >>> simplify(Ynm(1, -1, theta, phi).expand(func=True))
    sqrt(6)*E**(-I*phi)*sin(theta)/(4*sqrt(pi))

    >>> simplify(Ynm(1, 0, theta, phi).expand(func=True))
    sqrt(3)*cos(theta)/(2*sqrt(pi))

    >>> simplify(Ynm(1, 1, theta, phi).expand(func=True))
    -sqrt(6)*E**(I*phi)*sin(theta)/(4*sqrt(pi))

    >>> simplify(Ynm(2, -2, theta, phi).expand(func=True))
    sqrt(30)*E**(-2*I*phi)*sin(theta)**2/(8*sqrt(pi))

    >>> simplify(Ynm(2, -1, theta, phi).expand(func=True))
    sqrt(30)*E**(-I*phi)*sin(2*theta)/(8*sqrt(pi))

    >>> simplify(Ynm(2, 0, theta, phi).expand(func=True))
    sqrt(5)*(3*cos(theta)**2 - 1)/(4*sqrt(pi))

    >>> simplify(Ynm(2, 1, theta, phi).expand(func=True))
    -sqrt(30)*E**(I*phi)*sin(2*theta)/(8*sqrt(pi))

    >>> simplify(Ynm(2, 2, theta, phi).expand(func=True))
    sqrt(30)*E**(2*I*phi)*sin(theta)**2/(8*sqrt(pi))

    We can differentiate the functions with respect
    to both angles

    >>> theta = Symbol('theta')
    >>> phi = Symbol('phi')

    >>> diff(Ynm(n, m, theta, phi), theta)
    m*cot(theta)*Ynm(n, m, theta, phi) + E**(-I*phi)*sqrt((-m + n)*(m + n + 1))*Ynm(n, m + 1, theta, phi)

    >>> diff(Ynm(n, m, theta, phi), phi)
    I*m*Ynm(n, m, theta, phi)

    Further we can compute the complex conjugation

    >>> theta = Symbol('theta')
    >>> phi = Symbol('phi')
    >>> m = Symbol('m')

    >>> conjugate(Ynm(n, m, theta, phi))
    (-1)**(2*m)*E**(-2*I*m*phi)*Ynm(n, m, theta, phi)

    To get back the well known expressions in spherical
    coordinates we use full expansion

    >>> theta = Symbol('theta')
    >>> phi = Symbol('phi')

    >>> expand_func(Ynm(n, m, theta, phi))
    E**(I*m*phi)*sqrt((2*n + 1)*factorial(-m + n)/factorial(m + n))*assoc_legendre(n, m, cos(theta))/(2*sqrt(pi))

    See Also
    ========

    diofant.functions.special.spherical_harmonics.Ynm_c
    diofant.functions.special.spherical_harmonics.Znm

    References
    ==========

    * https://en.wikipedia.org/wiki/Spherical_harmonics
    * https://mathworld.wolfram.com/SphericalHarmonic.html
    * http://functions.wolfram.com/Polynomials/SphericalHarmonicY/
    * https://dlmf.nist.gov/14.30

    """

    @classmethod
    def eval(cls, n, m, theta, phi):
        # Handle negative index m and arguments theta, phi
        if m.could_extract_minus_sign():
            m = -m
            return (-1)**m * exp(-2*I*m*phi) * Ynm(n, m, theta, phi)
        if theta.could_extract_minus_sign():
            theta = -theta
            return Ynm(n, m, theta, phi)
        if phi.could_extract_minus_sign():
            phi = -phi
            return exp(-2*I*m*phi) * Ynm(n, m, theta, phi)

        # TODO Add more simplififcation here

    def _eval_expand_func(self, **hints):
        n, m, theta, phi = self.args
        rv = (sqrt((2*n + 1)/(4*pi) * factorial(n - m)/factorial(n + m)) *
              exp(I*m*phi) * assoc_legendre(n, m, cos(theta)))
        # We can do this because of the range of theta
        return rv.subs({sqrt(-cos(theta)**2 + 1): sin(theta)})

    def fdiff(self, argindex=4):
        if argindex == 3:
            # Diff wrt theta
            n, m, theta, phi = self.args
            return (m * cot(theta) * Ynm(n, m, theta, phi) +
                    sqrt((n - m)*(n + m + 1)) * exp(-I*phi) * Ynm(n, m + 1, theta, phi))
        if argindex == 4:
            # Diff wrt phi
            n, m, theta, phi = self.args
            return I * m * Ynm(n, m, theta, phi)
        # diff wrt n, m, etc
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_sin(self, n, m, theta, phi):
        return self.rewrite(cos)

    def _eval_rewrite_as_cos(self, n, m, theta, phi):
        # This method can be expensive due to extensive use of simplification!
        from ...simplify import simplify, trigsimp

        # TODO: Make sure n \in N
        # TODO: Assert |m| <= n ortherwise we should return 0
        term = simplify(self.expand(func=True))
        # We can do this because of the range of theta
        term = term.xreplace({Abs(sin(theta)): sin(theta)})
        return simplify(trigsimp(term))

    def _eval_conjugate(self):
        # TODO: Make sure theta \in R and phi \in R
        n, m, theta, phi = self.args
        return (-1)**m * self.func(n, -m, theta, phi)

    def as_real_imag(self, deep=True, **hints):
        # TODO: Handle deep and hints
        n, m, theta, phi = self.args
        re = (sqrt((2*n + 1)/(4*pi) * factorial(n - m)/factorial(n + m)) *
              cos(m*phi) * assoc_legendre(n, m, cos(theta)))
        im = (sqrt((2*n + 1)/(4*pi) * factorial(n - m)/factorial(n + m)) *
              sin(m*phi) * assoc_legendre(n, m, cos(theta)))
        return re, im


def Ynm_c(n, m, theta, phi):
    r"""Conjugate spherical harmonics defined as

    .. math::
        \overline{Y_n^m(\theta, \varphi)} := (-1)^m Y_n^{-m}(\theta, \varphi)

    See Also
    ========

    diofant.functions.special.spherical_harmonics.Ynm
    diofant.functions.special.spherical_harmonics.Znm

    References
    ==========

    * https://en.wikipedia.org/wiki/Spherical_harmonics
    * https://mathworld.wolfram.com/SphericalHarmonic.html
    * http://functions.wolfram.com/Polynomials/SphericalHarmonicY/

    """
    from .. import conjugate
    return conjugate(Ynm(n, m, theta, phi))


class Znm(Function):
    r"""
    Real spherical harmonics defined as

    .. math::

        Z_n^m(\theta, \varphi) :=
        \begin{cases}
          \frac{Y_n^m(\theta, \varphi) + \overline{Y_n^m(\theta, \varphi)}}{\sqrt{2}} &\quad m > 0 \\
          Y_n^m(\theta, \varphi) &\quad m = 0 \\
          \frac{Y_n^m(\theta, \varphi) - \overline{Y_n^m(\theta, \varphi)}}{i \sqrt{2}} &\quad m < 0 \\
        \end{cases}

    which gives in simplified form

    .. math::

        Z_n^m(\theta, \varphi) =
        \begin{cases}
          \frac{Y_n^m(\theta, \varphi) + (-1)^m Y_n^{-m}(\theta, \varphi)}{\sqrt{2}} &\quad m > 0 \\
          Y_n^m(\theta, \varphi) &\quad m = 0 \\
          \frac{Y_n^m(\theta, \varphi) - (-1)^m Y_n^{-m}(\theta, \varphi)}{i \sqrt{2}} &\quad m < 0 \\
        \end{cases}

    See Also
    ========

    diofant.functions.special.spherical_harmonics.Ynm
    diofant.functions.special.spherical_harmonics.Ynm_c

    References
    ==========

    * https://en.wikipedia.org/wiki/Spherical_harmonics
    * https://mathworld.wolfram.com/SphericalHarmonic.html
    * http://functions.wolfram.com/Polynomials/SphericalHarmonicY/

    """

    @classmethod
    def eval(cls, n, m, theta, phi):
        if m.is_positive:
            zz = (Ynm(n, m, theta, phi) + Ynm_c(n, m, theta, phi))/sqrt(2)
            return zz
        if m.is_zero:
            return Ynm(n, m, theta, phi)
        if m.is_negative:
            zz = (Ynm(n, m, theta, phi) - Ynm_c(n, m, theta, phi))/(sqrt(2)*I)
            return zz
