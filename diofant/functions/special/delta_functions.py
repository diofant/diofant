from ...core import Eq, Function, Integer, Rational, diff
from ...core.function import ArgumentIndexError
from ...core.sympify import sympify
from ...polys.polyerrors import PolynomialError
from ..elementary.complexes import im, sign
from ..elementary.piecewise import Piecewise


###############################################################################
# ############################## DELTA FUNCTION ############################# #
###############################################################################


class DiracDelta(Function):
    """
    The DiracDelta function and its derivatives.

    DiracDelta function has the following properties:

    1) ``diff(Heaviside(x),x) = DiracDelta(x)``
    2) ``integrate(DiracDelta(x-a)*f(x),(x,-oo,oo)) = f(a)`` and
       ``integrate(DiracDelta(x-a)*f(x),(x,a-e,a+e)) = f(a)``
    3) ``DiracDelta(x) = 0`` for all ``x != 0``
    4) ``DiracDelta(g(x)) = Sum_i(DiracDelta(x-x_i)/abs(g'(x_i)))``
       Where ``x_i``-s are the roots of ``g``

    Derivatives of ``k``-th order of DiracDelta have the following property:

    5) ``DiracDelta(x,k) = 0``, for all ``x != 0``

    See Also
    ========

    diofant.functions.special.delta_functions.Heaviside
    diofant.simplify.simplify.simplify
    diofant.functions.special.delta_functions.DiracDelta.is_simple
    diofant.functions.special.tensor_functions.KroneckerDelta

    References
    ==========

    * https://mathworld.wolfram.com/DeltaFunction.html

    """

    is_commutative = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            # I didn't know if there is a better way to handle default arguments
            k = 0
            if len(self.args) > 1:
                k = self.args[1]
            return self.func(self.args[0], k + 1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg, k=0):
        k = sympify(k)
        if not k.is_Integer or k.is_negative:
            raise ValueError('Error: the second argument of DiracDelta must be '
                             f'a non-negative integer, {k} given instead.')
        arg = sympify(arg)
        if arg.is_nonzero:
            return Integer(0)

    def simplify(self, x):
        """simplify(self, x)

        Compute a simplified representation of the function using
        property number 4.

        x can be:

        - a symbol

        Examples
        ========

        >>> DiracDelta(x*y).simplify(x)
        DiracDelta(x)/Abs(y)
        >>> DiracDelta(x*y).simplify(y)
        DiracDelta(y)/Abs(x)

        >>> DiracDelta(x**2 + x - 2).simplify(x)
        DiracDelta(x - 1)/3 + DiracDelta(x + 2)/3

        See Also
        ========

        diofant.functions.special.delta_functions.DiracDelta.is_simple
        diofant.functions.special.delta_functions.DiracDelta

        """
        from ...polys import roots

        if not self.args[0].has(x) or (len(self.args) > 1 and self.args[1] != 0 ):
            return self
        try:
            argroots = roots(self.args[0], x)
            result = 0
            valid = True
            darg = abs(diff(self.args[0], x))
            for r, m in argroots.items():
                if r.is_extended_real is not False and m == 1:
                    result += self.func(x - r)/darg.subs({x: r})
                else:
                    # don't handle non-real and if m != 1 then
                    # a polynomial will have a zero in the derivative (darg)
                    # at r
                    valid = False
                    break
            if valid:
                return result
        except PolynomialError:
            pass
        return self

    def is_simple(self, x):
        """is_simple(self, x)

        Tells whether the argument(args[0]) of DiracDelta is a linear
        expression in x.

        x can be:

        - a symbol

        Examples
        ========

        >>> DiracDelta(x*y).is_simple(x)
        True
        >>> DiracDelta(x*y).is_simple(y)
        True

        >>> DiracDelta(x**2+x-2).is_simple(x)
        False

        >>> DiracDelta(cos(x)).is_simple(x)
        False

        See Also
        ========

        diofant.simplify.simplify.simplify
        diofant.functions.special.delta_functions.DiracDelta

        """
        p = self.args[0].as_poly(x)
        if p:
            return p.degree() == 1
        return False

    @staticmethod
    def _latex_no_arg(printer):
        return r'\delta'

    def _eval_adjoint(self):
        return self

    def _eval_conjugate(self):
        return self

    def _eval_transpose(self):
        return self


###############################################################################
# ############################ HEAVISIDE FUNCTION ########################### #
###############################################################################


class Heaviside(Function):
    r"""Heaviside step function

    .. math ::
        H(x) = \left\{\begin{matrix}0, x < 0\\
                      1/2, x = 0\\ 1, x > 0 \end{matrix}\right.

    See Also
    ========

    diofant.functions.special.delta_functions.DiracDelta

    References
    ==========

    * https://en.wikipedia.org/wiki/Heaviside_step_function

    """

    is_real = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return DiracDelta(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)
        if im(arg).is_nonzero:
            raise ValueError(f'Function defined only for Real Values. Complex part: {im(arg)!r}  found in {arg!r} .' )
        elif arg.is_negative:
            return Integer(0)
        elif arg.is_zero:
            return Rational(1, 2)
        elif arg.is_positive:
            return Integer(1)

    def _eval_rewrite_as_Piecewise(self, arg):
        if arg.is_extended_real:
            return Piecewise((1, arg > 0), (Rational(1, 2), Eq(arg, 0)), (0, True))

    def _eval_rewrite_as_sign(self, arg):
        if arg.is_extended_real:
            return (sign(arg)+1)/2
