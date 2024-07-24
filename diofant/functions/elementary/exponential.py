from ...core import Add, E, Function, I, Integer, Mul, exp, log, oo, pi
from ...core.function import ArgumentIndexError, _coeff_isneg
from .miscellaneous import sqrt


class exp_polar(Function):
    r"""
    Represent a 'polar number' (see g-function Sphinx documentation).

    ``exp_polar`` represents the function
    `Exp: \mathbb{C} \rightarrow \mathcal{S}`, sending the complex number
    `z = a + bi` to the polar number `r = exp(a), \theta = b`. It is one of
    the main functions to construct polar numbers.

    The main difference is that polar numbers don't "wrap around" at `2 \pi`:

    >>> exp(2*pi*I)
    1
    >>> exp_polar(2*pi*I)
    exp_polar(2*I*pi)

    apart from that they behave mostly like classical complex numbers:

    >>> exp_polar(2)*exp_polar(3)
    exp_polar(5)

    See also
    ========

    diofant.simplify.powsimp.powsimp
    diofant.functions.elementary.complexes.polar_lift
    diofant.functions.elementary.complexes.periodic_argument
    diofant.functions.elementary.complexes.principal_branch

    """

    is_polar = True
    is_comparable = False  # cannot be evalf'd

    unbranched = True

    def _eval_as_numer_denom(self):
        """
        Returns this with a positive exponent as a 2-tuple (a fraction).

        Examples
        ========

        >>> exp(-x).as_numer_denom()
        (1, E**x)
        >>> exp(x).as_numer_denom()
        (E**x, 1)

        """
        # this should be the same as Pow.as_numer_denom wrt
        # exponent handling
        exp = self.exp
        neg_exp = exp.is_negative
        if not neg_exp and not (-exp).is_negative:
            neg_exp = _coeff_isneg(exp)
        if neg_exp:
            return Integer(1), self.func(-exp)
        return self, Integer(1)

    @property
    def exp(self):
        """Returns the exponent of the function."""
        return self.args[0]

    def _eval_conjugate(self):
        return self.func(self.exp.conjugate())

    def _eval_is_finite(self):
        arg = self.exp
        if arg.is_infinite:
            if arg.is_positive:
                return False
            if arg.is_negative:
                return True
        if arg.is_finite:
            return True

    def _eval_is_rational(self):
        if self.exp == 0:
            return True
        if self.exp.is_rational and self.exp.is_nonzero:
            return False

    def _eval_is_zero(self):
        if self.exp.is_infinite and self.exp.is_negative:
            return True

    def _eval_expand_power_exp(self, **hints):
        arg = self.exp
        if arg.is_Add and arg.is_commutative:
            expr = 1
            for x in arg.args:
                expr *= self.func(x)
            return expr
        return self.func(arg)

    def _eval_Abs(self):
        from ...core import expand_mul
        return sqrt(expand_mul(self * self.conjugate()))

    def _eval_evalf(self, prec):
        """Careful! any evalf of polar numbers is flaky."""
        from .complexes import im, re
        i = im(self.exp)
        try:
            bad = (i <= -pi or i > pi)
        except TypeError:
            bad = True
        if bad:
            return self  # cannot evalf for this argument
        res = exp(self.exp).evalf(prec, strict=False)
        if i > 0 > im(res):
            # i ~ pi, but exp(I*i) evaluated to argument slightly bigger than pi
            return re(res)
        return res

    def _eval_power(self, other):
        return self.func(self.exp*other)

    def _eval_is_extended_real(self):
        if self.exp.is_extended_real:
            return True

    def as_base_exp(self):
        if self.exp == 0:
            return super().as_base_exp()
        return self.func(1), Mul(*self.args)


class LambertW(Function):
    r"""
    The Lambert W function `W(z)` is defined as the inverse
    function of `w \exp(w)`.

    In other words, the value of `W(z)` is such that `z = W(z) \exp(W(z))`
    for any complex number `z`.  The Lambert W function is a multivalued
    function with infinitely many branches `W_k(z)`, indexed by
    `k \in \mathbb{Z}`.  Each branch gives a different solution `w`
    of the equation `z = w \exp(w)`.

    The Lambert W function has two partially real branches: the
    principal branch (`k = 0`) is real for real `z > -1/e`, and the
    `k = -1` branch is real for `-1/e < z < 0`. All branches except
    `k = 0` have a logarithmic singularity at `z = 0`.

    Examples
    ========

    >>> LambertW(1.2)
    0.635564016364870
    >>> LambertW(1.2, -1).evalf()
    -1.34747534407696 - 4.41624341514535*I
    >>> LambertW(-1).is_real
    False

    References
    ==========

    * https://en.wikipedia.org/wiki/Lambert_W_function

    """

    @classmethod
    def eval(cls, x, k=None):
        if k == 0:
            return cls(x)
        if k is None:
            k = Integer(0)

        if k == 0:
            if x == 0:
                return Integer(0)
            if x is E:
                return Integer(1)
            if x == -1/E:
                return Integer(-1)
            if x == -log(2)/2:
                return -log(2)
            if x is oo:
                return oo

        if k.is_nonzero:
            if x == 0:
                return -oo
        if k == -1:
            if x == -pi/2:
                return -I*pi/2
            if x == -1/E:
                return Integer(-1)
            if x == -2*exp(-2):
                return -Integer(2)

    def fdiff(self, argindex=1):
        """Return the first derivative of this function."""
        x = self.args[0]

        if len(self.args) == 1:
            if argindex == 1:
                return LambertW(x)/(x*(1 + LambertW(x)))
            raise ArgumentIndexError(self, argindex)
        k = self.args[1]
        if argindex == 1:
            return LambertW(x, k)/(x*(1 + LambertW(x, k)))
        raise ArgumentIndexError(self, argindex)

    def _eval_is_extended_real(self):
        x = self.args[0]
        if len(self.args) == 1:
            k = Integer(0)
        else:
            k = self.args[1]
        if k.is_zero:
            if (x + 1/E).is_positive:
                return True
            if (x + 1/E).is_nonpositive:
                return False
        elif (k + 1).is_zero:
            if x.is_negative and (x + 1/E).is_nonnegative:
                return True
            if x.is_nonpositive or (x + 1/E).is_positive:
                return False

    def _eval_is_algebraic(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if self.args[0].is_nonzero and self.args[0].is_algebraic:
                return False
        else:
            return s.is_algebraic

    def _eval_nseries(self, x, n, logx):
        if len(self.args) == 1:
            from ...calculus import Order
            from .. import factorial
            x = self.args[0]
            o = Order(x**n, x)
            l = Integer(0)
            if n > 0:
                l += Add(*[Integer(-k)**(k - 1)*x**k/factorial(k)
                           for k in range(1, n)])
            return l + o
        return super()._eval_nseries(x, n, logx)
