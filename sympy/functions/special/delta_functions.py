from sympy.core import S, sympify, diff
from sympy.core.function import Function, ArgumentIndexError
from sympy.core.relational import Eq
from sympy.polys.polyerrors import PolynomialError
from sympy.functions.elementary.complexes import im, sign
from sympy.functions.elementary.piecewise import Piecewise

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

    sympy.functions.special.delta_functions.Heaviside
    sympy.simplify.simplify.simplify
    sympy.functions.special.delta_functions.DiracDelta.is_simple
    sympy.functions.special.tensor_functions.KroneckerDelta

    References
    ==========

    .. [1] http://mathworld.wolfram.com/DeltaFunction.html
    """

    is_extended_real = True

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
            raise ValueError("Error: the second argument of DiracDelta must be \
            a non-negative integer, %s given instead." % (k,))
        arg = sympify(arg)
        if arg is S.NaN:
            return S.NaN
        if arg.is_positive or arg.is_negative:
            return S.Zero

    def simplify(self, x):
        """simplify(self, x)

           Compute a simplified representation of the function using
           property number 4.

           x can be:

           - a symbol

           Examples
           ========

           >>> from sympy import DiracDelta
           >>> from sympy.abc import x, y

           >>> DiracDelta(x*y).simplify(x)
           DiracDelta(x)/Abs(y)
           >>> DiracDelta(x*y).simplify(y)
           DiracDelta(y)/Abs(x)

           >>> DiracDelta(x**2 + x - 2).simplify(x)
           DiracDelta(x - 1)/3 + DiracDelta(x + 2)/3

           See Also
           ========

           sympy.functions.special.delta_functions.DiracDelta.is_simple
           sympy.functions.special.delta_functions.DiracDelta

        """
        from sympy.polys.polyroots import roots

        if not self.args[0].has(x) or (len(self.args) > 1 and self.args[1] != 0 ):
            return self
        try:
            argroots = roots(self.args[0], x)
            result = 0
            valid = True
            darg = abs(diff(self.args[0], x))
            for r, m in argroots.items():
                if r.is_extended_real is not False and m == 1:
                    result += self.func(x - r)/darg.subs(x, r)
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

           >>> from sympy import DiracDelta, cos
           >>> from sympy.abc import x, y

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

           sympy.simplify.simplify.simplify
           sympy.functions.special.delta_functions.DiracDelta
        """
        p = self.args[0].as_poly(x)
        if p:
            return p.degree() == 1
        return False

    @staticmethod
    def _latex_no_arg(printer):
        return r'\delta'


###############################################################################
# ############################ HEAVISIDE FUNCTION ########################### #
###############################################################################


class Heaviside(Function):
    """Heaviside Piecewise function

    Heaviside function has the following properties [*]_:

    1) ``diff(Heaviside(x),x) = DiracDelta(x)``
                        ``( 0, if x < 0``
    2) ``Heaviside(x) = < ( 1/2 if x==0 [*]``
                        ``( 1, if x > 0``

    .. [*] Regarding to the value at 0, Mathematica defines ``H(0) = 1``,
           but Maple uses ``H(0) = undefined``

    I think is better to have H(0) = 1/2, due to the following::

        integrate(DiracDelta(x), x) = Heaviside(x)
        integrate(DiracDelta(x), (x, -oo, oo)) = 1

    and since DiracDelta is a symmetric function,
    ``integrate(DiracDelta(x), (x, 0, oo))`` should be 1/2 (which is what
    Maple returns).

    If we take ``Heaviside(0) = 1/2``, we would have
    ``integrate(DiracDelta(x), (x, 0, oo)) = ``
    ``Heaviside(oo) - Heaviside(0) = 1 - 1/2 = 1/2``
    and
    ``integrate(DiracDelta(x), (x, -oo, 0)) = ``
    ``Heaviside(0) - Heaviside(-oo) = 1/2 - 0 = 1/2``

    If we consider, instead ``Heaviside(0) = 1``, we would have
    ``integrate(DiracDelta(x), (x, 0, oo)) = Heaviside(oo) - Heaviside(0) = 0``
    and
    ``integrate(DiracDelta(x), (x, -oo, 0)) = Heaviside(0) - Heaviside(-oo) = 1``

    See Also
    ========

    sympy.functions.special.delta_functions.DiracDelta

    References
    ==========

    .. [1] http://mathworld.wolfram.com/HeavisideStepFunction.html

    """

    is_real = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            # property number 1
            return DiracDelta(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)
        if arg is S.NaN:
            return S.NaN
        elif im(arg).is_nonzero:
            raise ValueError("Function defined only for Real Values. Complex part: %s  found in %s ." % (repr(im(arg)), repr(arg)) )
        elif arg.is_negative:
            return S.Zero
        elif arg.is_zero:
            return S.Half
        elif arg.is_positive:
            return S.One

    def _eval_rewrite_as_Piecewise(self, arg):
        if arg.is_extended_real:
            return Piecewise((1, arg > 0), (S(1)/2, Eq(arg, 0)), (0, True))

    def _eval_rewrite_as_sign(self, arg):
        if arg.is_extended_real:
            return (sign(arg)+1)/2
