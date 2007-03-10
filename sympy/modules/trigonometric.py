from sympy.core.functions import Function
from sympy.core.numbers import Real, Rational, pi
import decimal

class sin(Function):
    """
    Usage
        sin(x) -> Returns the sine of x (measured in radians)
        
    Notes
        sin(x) will evaluate automatically in the case x is a 
        multiple of pi.
    
    Further examples
        >>> sin(x**2).diff(x)
        2*cos(x^2)*x
        
        >>> sin(1).diff(x)
        0
        >>> sin(pi)
        0
        
    See also: 
       cos, tan
       Definitions in trigonometry: http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html
    """
    
    def getname(self):
        return "sin"
        
    def derivative(self):
        return cos(self.arg)

    def bounded(self):
        return True
        
    def eval(self):
        if not self.arg.isnumber():
             return self
        a = 2*self.arg / pi
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
        if not self.arg.isnumber():
            raise ValueError("Argument can't be a symbolic value")
        decimal.getcontext().prec = precision + 2
        x = Real(self.arg)
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
    
class cos(Function):
    """Return the cosine of x (measured in radians)
    """
    
    def getname(self):
        return "cos"
        
    def derivative(self):
        return -sin(self.arg)

    def bounded(self):
        return True
    
    def eval(self):
        if not self.arg.isnumber():
             return self
        # case self.arg is a number 
        a = 2*self.arg / pi
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
        if not self.arg.isnumber():
            raise ValueError("Argument can't be a symbolic value")
        decimal.getcontext().prec = precision + 2
        x = Real(self.arg)
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

class tan(Function):
    """Return the tangent of x (measured in radians)
    """
    
    def getname(self):
        return "tan"
        
    def derivative(self):
        return Rational(1) / (cos(self.arg)**2)
        
    def eval(self):
        return sin(self.arg) / cos(self.arg)

    
