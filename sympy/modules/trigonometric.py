from sympy.core.functions import Function, exp, sqrt
from sympy.core.numbers import Real, Rational, pi, I
from sympy.core import Symbol

import decimal
import math

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

    def eval(self):
        if self._args.is_number:
            if self._args == 0:
                return Rational(0)
            else:
                a = Symbol('a')

                coeff = self._args.match(a*pi, [a])

                if coeff != None:
                    arg = coeff [a]

                    if isinstance(arg, int):
                        return Rational(0)
                    elif isinstance(arg, Rational):
                        if arg.is_integer:
                            return Rational(0)
                        else:
                            if arg.q == 2:
                                result = Rational(1)
                            elif arg.q == 3:
                                result = Rational(1, 2)*sqrt(3)
                            elif arg.q == 4:
                                result = Rational(1, 2)*sqrt(2)
                            elif arg.q == 6:
                                result = Rational(1, 2)
                            else:
                                return self

                            if (arg.p / arg.q) % 2 == 1:
                                result *= -1

                            return result

        return self

    def evalf(self, precision=28):
        if not self._args.is_number:
            raise ValueError("Argument can't be a symbolic value")
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
    
class cos(Function):
    """
    Usage
    =====
      cos(x) -> Returns the cosine of x (measured in radians)
        
    Notes
    =====
        cos(x) will evaluate automatically in the case x is a 
        multiple of pi.
    
    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> cos(x**2).diff(x)
        -2*x*sin(x**2)
        >>> cos(1).diff(x)
        0
        >>> cos(pi)
        -1
        
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
    
    def eval(self):
        if not self._args.is_number:
             return self
        # case self._args is a number 
        a = 2*self._args / pi
        if a - int(float(a)) == 0:
            # arg is a multiple of pi
            a_mod4 = int(a) % 4
            if a_mod4 == 1 or a_mod4 == 3:
                return Rational(0)
            else:
                if a_mod4 % 4 == 0: 
                    return Rational(1)
                elif a_mod4 == 2:
                    return Rational(-1)
        return self
    
    def evalf(self, precision=28):
        if not self._args.is_number:
            raise ValueError("Argument can't be a symbolic value")
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
        
    def eval(self):
        return self

    def evalf(self):
        return sin(self._args).evalf() / cos(self._args).evalf()

class asin(Function):
    """Return the arc sine of x (measured in radians)"""
    
    def derivative(self):
        return sqrt(1-self._args**2)**(-1)
    
    def eval(self):
        return self
    
def acos(Function):
    """Return the arc sine of x (measured in radians)"""
    
    def derivative(self):
        return - sqrt(1-self._args**2)**(-1)
    
    def eval(self):
        return self
    
class atan(Function):
    """Return the arc tangent of x (measured in radians)
    """

    def derivative(self):
        return Rational(1) / (1+(self._args)**2)
        
    def eval(self):
        if self._args == 0:
            return Rational(0)
        return self
