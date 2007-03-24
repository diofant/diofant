from sympy.core.functions import Function, exp
from sympy.core.numbers import Real, Rational, pi, I
import decimal

class sin(Function):
    """
    Usage
    =====
      sin(x) -> Returns the sine of x (measured in radians)
        
    Notes
    =====
        sin(x) will evaluate automatically in the case x is a 
        multiple of pi.
    
    Further examples
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> sin(x**2).diff(x)
        2*cos(x**2)*x
        >>> sin(1).diff(x)
        0
        >>> sin(pi)
        0
        
    See also
    ========
       L{cos}, L{tan}
       
       External links
       --------------
         
         U{Definitions in trigonometry<http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html>}
    """
    
    def derivative(self):
        return cos(self._args)

    def bounded(self):
        return True
        
    def eval(self):
        if not self._args.isnumber():
             return self
        a = 2*self._args / pi
        if a - int(float(a)) == 0:
            # arg is a multiple of pi
            a_mod4 = int(a) % 4
            if a_mod4 == 0 or a_mod4 == 2:
                return Rational(0)
            else:
                if a_mod4 % 4 == 1: 
	                   return Rational(1)
                elif a_mod4 == 3:
                    return Rational(-1)
        return self
    
    def evalf(self, precision=28):
        if not self._args.isnumber():
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
        return +s

    def evalc(self):
        x, y = self._args.get_re_im()
        sinh = (exp(y)-exp(-y))/2
        cosh = (exp(y)+exp(-y))/2
        return sin(x)*cosh + I*cos(x)*sinh
    
class cos(Function):
    """Return the cosine of x (measured in radians)
    """
    
    def derivative(self):
        return -sin(self._args)

    def bounded(self):
        return True
    
    def eval(self):
        if not self._args.isnumber():
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
        if not self._args.isnumber():
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
        return +s

    def evalc(self):
        x, y = self._args.get_re_im()
        sinh = (exp(y)-exp(-y))/2
        cosh = (exp(y)+exp(-y))/2
        return cos(x)*cosh - I*sin(x)*sinh

class tan(Function):
    """Return the tangent of x (measured in radians)
    """
    
    def derivative(self):
        return Rational(1) / (cos(self._args)**2)
        
    def eval(self):
        return sin(self._args) / cos(self._args)

    
class arctan(Function):
    """Return the tangent of x (measured in radians)
    """
    
    def derivative(self):
        return Rational(1) / (1+(self._args)**2)
        
    def eval(self):
        if self._args == 0:
            return Rational(0)
        return self
