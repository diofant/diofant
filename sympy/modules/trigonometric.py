from sympy.core.functions import Function, exp, sqrt
from sympy.core.numbers import Number, Real, Rational, pi, I, oo
from sympy import Symbol, Add, Mul
from simplify import ratsimp

import decimal
import math

def pi_ratmatch(x):
    """match x = Rational(p,q) * pi"""

    # try most common cases (?)
    if not x.is_number:
        return None

    a = Symbol('a')
    coeff = x.match(a*pi, [a])

    if coeff is None:
        return None

    arg = coeff[a]

    if isinstance(arg, Rational):
        return arg

    return None



class sin(Function):
    """
    Usage
    =====
      sin(x) -> Returns the sine of x (measured in radians)

    Notes
    =====
        sin(x) will evaluate automatically in the case x
        is a multiple of pi, pi/2, pi/3, pi/4 and pi/6.

    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> sin(x**2).diff(x)
        2*x*cos(x**2)
        >>> sin(1).diff(x)
        0
        >>> sin(pi)
        0
        >>> sin(pi/2)
        1
        >>> sin(pi/6)
        1/2

    See also
    ========
       L{cos}, L{tan}

       External links
       --------------

         U{Definitions in trigonometry<http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html>}
    """

    def derivative(self):
        return cos(self._args)

    @property
    def is_bounded(self):
        return True

    @property
    def is_real(self):
        return self._args.is_real

    def eval(self):
        # FIXME move to class level.
        # doesn't work there because of cyclic import dependency
        # NB: pi is ommited
        cst_table = {
            Rational(1):    Rational(0),
            Rational(1,2):  Rational(1),
            Rational(1,3):  Rational(1,2)*sqrt(3),
            Rational(1,4):  Rational(1,2)*sqrt(2),
            Rational(1,6):  Rational(1,2)
        }

        # sin is odd
        if self._args.is_number and self._args < 0:
            return -sin(-self._args)

        arg = pi_ratmatch(self._args)

        if arg is not None:
            p, q = arg.p, arg.q
            p = p % (2*q)   # sin period is 2*pi

            m = (-1)**(p//q)    # p < q --> +1

            try:
                return m * cst_table[ Rational(1,q) ]
            except KeyError:
                return self


        elif isinstance(self._args, Mul):
            if isinstance(self._args[0], Number) and self._args[0] < 0:
                return -sin(-self._args)

        return self

    def evalf(self, precision=18):
        if not self._args.is_number:
            raise ValueError("Argument can't be a symbolic value")
        if precision <= 18:
            return math.sin(self._args)
        decimal.getcontext().prec = precision + 2
        x = Real(self._args)
        i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
        while s != lasts:
            lasts = s
            i += 2
            fact *= i * (i-1)
            num *= x * x
            sign *= -1
            s += num / fact * sign
        decimal.getcontext().prec = precision - 2
        return s

    def evalc(self):
        x, y = self._args.get_re_im()
        sinh = (exp(y)-exp(-y))/2
        cosh = (exp(y)+exp(-y))/2
        return sin(x)*cosh + I*cos(x)*sinh

    def expand(self):
        if isinstance(self._args, Add):
            left = self._args[0]
            right = self._args[1:]
            if len(right) == 1:
                right = right[0]
            else:
                right = Add(*right)
            t1 = sin(left)*cos(right).expand()
            t2 = cos(left)*sin(right).expand()
            return (t1 + t2).expand()
        elif isinstance(self._args, Mul):
            n = self._args[0]
            if isinstance(n, Rational) and n.is_integer:
                # sin(nx) = 2 sin[(n-1)x] cos x - sin[(n-2)x]
                x = Mul(*self._args[1:])
                sign = 1
                if n < 0:
                    n = -n
                    sign = -1
                t1 = 2*sin((n-1)*x).expand()*cos(x).expand()
                t2 = sin((n-2)*x).expand()
                return (sign*(t1 - t2)).expand()
            else:
                return self
        else:
            return self

class cos(Function):
    """
    Usage
    =====
      cos(x) -> Returns the cosine of x (measured in radians)

    Notes
    =====
        cos(x) will evaluate automatically in the case x
        is a multiple of pi, pi/2, pi/3, pi/4 and pi/6.

    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> cos(x**2).diff(x)
        -2*sin(x**2)*x
        >>> cos(1).diff(x)
        0
        >>> cos(pi)
        -1
        >>> cos(pi/2)
        0
        >>> cos(2*pi/3)
        -1/2

    See also
    ========
       L{sin}, L{tan}

       External links
       --------------

         U{Definitions in trigonometry<http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html>}
    """


    def derivative(self):
        return -sin(self._args)

    @property
    def is_bounded(self):
        return True

    @property
    def is_real(self):
        return self._args.is_real

    def eval(self):
        # FIXME move to class level.
        # doesn't work there because of cyclic import dependency
        # NB: pi is ommited
        cst_table = {
            Rational(1):    Rational(1),    # XXX should be -1, but ok
            Rational(1,2):  Rational(0),
            Rational(1,3):  Rational(1,2),
            Rational(1,4):  Rational(1,2)*sqrt(2),
            Rational(1,6):  Rational(1,2)*sqrt(3)
        }
  
        # cos is even
        if self._args.is_number and self._args < 0:
            return cos(-self._args)
  
        arg = pi_ratmatch(self._args)
  
        if arg is not None:
            p, q = arg.p, arg.q
            p = p % (2*q)   # cos period is 2*pi

            octant = (2*p//q) % 4
            if octant in (1,2):
                m = -1
            else:
                m = +1
  
            try:
                return m * cst_table[ Rational(1,q) ]
            except KeyError:
                return self
        elif isinstance(self._args, Mul):
            if isinstance(self._args[0], Number) and self._args[0] < 0:
                return cos(-self._args)

        return self

    def evalf(self, precision=18):
        if not self._args.is_number:
            raise ValueError("Argument can't be a symbolic value")
        if precision <= 18:
            return math.cos(self._args)
        decimal.getcontext().prec = precision + 2
        x = Real(self._args)
        i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
        while s != lasts:
            lasts = s
            i += 2
            fact *= i * (i-1)
            num *= x * x
            sign *= -1
            s += num / fact * sign
        decimal.getcontext().prec = precision - 2
        return s

    def evalc(self):
        x, y = self._args.get_re_im()
        sinh = (exp(y)-exp(-y))/2
        cosh = (exp(y)+exp(-y))/2
        return cos(x)*cosh - I*sin(x)*sinh

    def expand(self):
        if isinstance(self._args, Add):
            left = self._args[0]
            right = self._args[1:]
            if len(right) == 1:
                right = right[0]
            else:
                right = Add(*right)
            t1 = cos(left)*cos(right).expand()
            t2 = sin(left)*sin(right).expand()
            return (t1 - t2).expand()
        elif isinstance(self._args, Mul):
            n = self._args[0]
            if isinstance(n, Rational) and n.is_integer:
                # cos nx = 2 cos[(n-1)x] cos x - cos[(n-2)x]
                n = self._args[0]
                x = Mul(*self._args[1:])
                if n < 0:
                    n = -n
                t1 = 2*cos((n-1)*x).expand()*cos(x).expand()
                t2 = cos((n-2)*x).expand()
                return (t1 - t2).expand()
            else:
                return self
        else:
            return self

class tan(Function):
    """
    Usage
    =====
      tan(x) -> Returns the tangent of x (measured in radians)

    Notes
    =====
        tan(x) will evaluate automatically in the case x is a
        multiple of pi.

    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> tan(x**2).diff(x)
        2*x*cos(x**2)**(-2)
        >>> tan(1).diff(x)
        0

    See also
    ========
       L{sin}, L{tan}

       External links
       --------------

         U{Definitions in trigonometry<http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html>}
    """


    def derivative(self):
        return Rational(1) / (cos(self._args)**2)

    @property
    def is_bounded(self):
        return False

    @property
    def is_real(self):
        return self._args.is_real

    def eval(self):
        s, c = sin(self._args).eval(), cos(self._args).eval()

        # XXX this is hackish
        # XXX what we really want is to check whether it is primitive constant
        if (not isinstance(s, sin)) and (not isinstance(c, cos)) and c != 0:
            return s/c

        if isinstance(self._args, Number) and self._args < 0:
            return -tan(-self._args)
        elif isinstance(self._args, Mul):
            if isinstance(self._args[0], Number) and self._args[0] < 0:
                return -tan(-self._args)
        return self

    def evalf(self):
        return sin(self._args).evalf() / cos(self._args).evalf()

    def expand(self):
        def expand_fraction(num, den):
            """
            A function to check to see if a fraction is of the form (a/b)/(c/d)
            and then expands (a*d) and (b*c) separately so that if the
            numerator is an Add instance, the fraction isn't broken up
            into multiple instances.
            """
            from sympy import Add,Mul,Pow
            if isinstance(den, Mul):
                a,b = den.getab()
                if isinstance(a, Pow) and a.exp == -1:
                    a,b = b,a
                if isinstance(b, Pow) and b.exp == -1:
                    ret = Rational(0)
                    den = a
                    if isinstance(num, Add):
                        ret = Rational(0)
                        for x in num:
                            ret += (x * b.base)
                        num = ret
                    else:
                        num *= b
            return num.expand() / den.expand()

        if isinstance(self._args, Add):
            left = self._args[0]
            right = self._args[1:]
            if len(right) == 1:
                right = right[0]
            else:
                right = Add(*right)
            a = tan(left).expand()
            b = tan(right).expand()
            t1 = (a + b).expand()
            t2 = ratsimp( (1 - a*b).expand() )
            return expand_fraction(t1, t2)
        elif isinstance(self._args, Mul) and isinstance(self._args[0], Rational):
            n = self._args[0]
            if isinstance(n, Rational) and n.is_integer:
                # tan nx = (tan[(n-1)x] + tan[x]) / (1 - tan[(n-1)x] tan[x])
                x = Mul(*self._args[1:])
                sign = 1
                if n < 0:
                    n = -n
                    sign = -1

                a = tan((n-1)*x).expand()
                b = tan(x).expand()
                t1 = (sign*(a + b)).expand()
                t2 = ratsimp( (1 - a*b).expand() )
                return expand_fraction(t1, t2)
            else:
                return self
        else:
            return self

def sec(x):
    return 1/cos(x)

def csc(x):
    return 1/sin(x)

def cot(x):
    return 1/tan(x)

class asin(Function):
    """Return the arc sine of x (measured in radians)"""

    def derivative(self):
        return sqrt(1-self._args**2)**(-1)


    def eval(self):
        # FIXME move to class level.
        # doesn't work there because of cyclic import dependency
        cst_table = {
            Rational(1):    pi/2,
            sqrt(3)/2:      pi/3,
            sqrt(2)/2:      pi/4,
            Rational(1,2):  pi/6,
            Rational(0):    Rational(0),
        }

        # asin is odd
        if self._args.is_number and self._args < 0:
            return -asin(-self._args)

        try:
            return cst_table[ self._args ]
        except KeyError:
            return self

class acos(Function):
    """Return the arc sine of x (measured in radians)"""

    def derivative(self):
        return - sqrt(1-self._args**2)**(-1)

    def eval(self):
        # FIXME move to class level.
        # doesn't work there because of cyclic import dependency
        cst_table = {
            Rational(1):            Rational(0),
            Rational(1,2)*sqrt(3):  pi/6,
            Rational(1,2)*sqrt(2):  pi/4,
            Rational(1,2):          pi/3,
            Rational(0):            pi/2
        }

        print 'acos: ', self._args

        # acos(-x) = pi - acos(x)
        if self._args.is_number and self._args < 0:
            return pi - acos(-self._args)

        try:
            return cst_table[ self._args ]
        except KeyError:
            return self

class atan(Function):
    """Return the arc tangent of x (measured in radians)"""

    def derivative(self):
        return Rational(1) / (1+(self._args)**2)

    def series(self, x, n=6):
        #expanding around infinity can actually be achieved for atan:
        if self._args == 1/x:
            return pi/2
        return Function.series(self, x, n)

    def eval(self):
        # FIXME move to class level.
        # doesn't work there because of cyclic import dependency
        cst_table = {
            oo:             pi/2,
            sqrt(3):        pi/3,
            Rational(1):    pi/4,
            1/sqrt(3):      pi/6,
            Rational(0):    Rational(0)
        }

        # atan is odd
        if self._args.is_number and self._args < 0:
            return -atan(-self._args)
        elif isinstance(self._args, Mul):
            if isinstance(self._args[0], Number) and self._args[0] < 0:
                return -atan(-self._args)

        try:
            return cst_table[ self._args ]
        except KeyError:
            return self
