from sympy.core import hashing
from sympy.core.basic import Basic
import decimal

class Number(Basic):
    """Represents any kind of number in sympy.


    Floating point numbers are represented by the Real class.
    Integer numbers (of any size), together with rational numbers (again, there
    is no limit on their size) are represented by the Rational class. 

    If you want to represent for example 1+sqrt(2), then you need to do:

    Rational(1) + Rational(2)**( Rational(1)/2 )
    """
    
    mathml_tag = "cn"
    
    def __init__(self):
        Basic.__init__(self, is_commutative = True)
        
    def __int__(self):
        raise NotImplementedError
    
    def __len__(self):
        return 1
    
    def __float__(self):
        return float(self.evalf())
    
    def __abs__(self):
        from functions import abs_
        return abs_(self)
    
    def diff(self,sym):
        return Rational(0)
    
    def evalf(self):
        raise NotImplementedError

    def evalc(self):
        return self
    
class Infinity(Number):
    """
    Usage
    =====
        Represents mathematical infinity. 
        
    Notes
    =====
        Cannot be used in expressions like infty/infty, but can be used in some
        very simple expressions like 1*infty
          
        Should be used only as a Symbol, for example results of limits, integration limits etc.
        Can however be used in comparisons, like infty!=1, or infty!=x**3
          
    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> limit(x, x, infty)
        Infinity
    """
    
    def __init__(self, sign=1):
        Basic.__init__(self, 
                       is_real = False, 
                       is_commutative = False, 
                       )
        
    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        return self._mhash.value
    
    def sign(self):
        return self._sign
    
    def __lt__(self, num):
        if self.sympify(num).is_number:
            if self._sign == -1:
                return True
            else:
                return False
    
    def __gt__(self, num):
        return not self.__lt__(num)

infty = Infinity()

class Real(Number):
    """Represents a floating point number. It is capable of representing
    arbitrary-precision floating-point numbers

    Usage:

    Real(3.5)   .... 3.5 (the 3.5 was converted from a python float)
    Real("3.0000000000000005")
    
    """
    
    
    def __init__(self,num):
        Basic.__init__(self, 
                        is_real = True, 
                        is_commutative = True, 
                        )
        if isinstance(num,str):
            num = decimal.Decimal(num)
        if isinstance(num, decimal.Decimal):
            self.num = num
        elif isinstance(num, Real):
            self.num = num.evalf()
        else:
            self.num = decimal.Decimal(str(float(num)))
        
    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash=hashing.mhash()
        self._mhash.addstr(str(type(self)))
        self._mhash.addfloat(self.num)
        return self._mhash.value
        
    def __str__(self):
        if self.num < 0:
            f = "(%s)"
        else:
            f = "%s"
        return f % (str(self.num))

    def __float__(self):
        return float(self.num)
        
    def __int__(self):
        return int(self.evalf())
    
    def __add__(self,a):
        if isnumber(a):
            if isinstance(a, Real):
                return Real(self.num + a.num)
            else:
                return Real(self.num + decimal.Decimal(str(float(a))))
        else:
            assert isinstance(a, Basic)
            from addmul import Add
            return Add(self, a)
        
    def __mul__(self,a):
        if Basic.sympify(a).is_number:
            return Real(self.num * decimal.Decimal(str(float(a))))
            #FIXME: too many boxing-unboxing
        else:
            assert isinstance(a, Basic)
            from addmul import Mul
            return Mul(self, a)
        
    def __pow__(self,a):
        from power import Pow
        return Pow(self, a)
        
    def __rpow__(self, a):
        from power import Pow
        return Pow(a, self)
        
    def isone(self):
        if self.num == 1:
            return True
        else:
            return False
        
    @property
    def is_integer(self):
        return int(self) - self.evalf() == 0
        
    def evalf(self):
        #evalf() should return either a float or an exception
        return self.num


class Rational(Number):
    """Represents integers and rational numbers (p/q) of any size.

    Thanks to support of long ints in Python. 

    Usage:

    Rational(3)      ... 3
    Rational(1,2)    ... 1/2
    """
    
    def __init__(self,*args):
        Basic.__init__(self, 
                       is_real = True, 
                       is_commutative = True, 
                       )
        if len(args)==1:
            p = args[0]
            q = 1 
        elif len(args)==2:
            p = args[0]
            q = args[1]
        else:
            raise "invalid number of arguments"
        assert (isinstance(p, int) or isinstance(p, long)) and \
               (isinstance(q, int) or isinstance(q, long))
        assert q != 0
        s = sign(p)*sign(q)
        p = abs(p)
        q = abs(q)
        c = self.gcd(p,q)
        self.p = p/c*s
        self.q = q/c
        
    def sign(self):
        return sign(self.p)*sign(self.q)
        
    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        self._mhash.addint(self.p)
        self._mhash.addint(self.q)
        return self._mhash.value
        
    def gcd(self,a,b):
        """Primitive algorithm for a greatest common divisor of "a" and "b"."""
        while b:
            a, b = b, a % b
        return a
        
    def __str__(self):
        if self.q == 1:
            f = "%d"
            return f % (self.p)
        else:
            f = "%d/%d"
            return f % (self.p,self.q)

    def __mul__(self,a):
        a = self.sympify(a)
        if isinstance(a, Rational):
            return Rational(self.p * a.p, self.q * a.q)
        elif isinstance(a, int) or isinstance(a, long):
            return Rational(self.p * a, self.q)
        elif isinstance(a, Real):
            return a.__mul__(self)
        else:
            from addmul import Mul
            return Mul(self, a)
    
    def __rmul__(self, a):
        return self.__mul__(a)
    
    def __div__(self, a):
        #TODO: move to Mul.eval
        if isinstance(a, int):
            return Rational(self.p, self.q *a)
        return self * (a**Rational(-1))
        
    def __rdiv__(self, a):
        #TODO: move to Mul.eval
        if isinstance(a, int):
            return Rational(self.q * a, self.p )
        return self * (a**Rational(-1))
    
    def __add__(self,a):
        a=self.sympify(a)
        if isinstance(a, Rational):
            return Rational(self.p*a.q+self.q*a.p,self.q*a.q)
        elif isinstance(a, int) or isinstance(a, long):
            return Rational(self.p + a*self.q, self.q)
        elif isinstance(a, Real):
            return a.__add__(self)
        else:
            from addmul import Add
            return Add(self, a)
        
    def __pow__(self,a):
        """Returns the self to the power of "a"
        """
        from power import Pow
        return Pow(self, a)
    
    def __rpow__(self, a):  
        """Returns "a" to the power of self
        """
        from power import Pow
        return Pow(a, self)
    

    def __int__(self):
        assert self.is_integer
        return self.p
    
    def iszero(self):
        return self.p == 0 
        
    def isone(self):
        return self.p == 1 and self.q == 1
        
    def isminusone(self):
        return self.p == -1 and self.q == 1
        
    @property
    def is_integer(self):
        """Returns True if the current number is an integer
        and False otherwise. 
        
        Examples
        ========
            >>> Rational(1).is_integer
            True
            >>> Rational(1,2).is_integer
            False
            
        """
        return self.q == 1
        
    def evalf(self):
        return decimal.Decimal(self.p)/self.q
        
    def diff(self,sym):
        return Rational(0)

    def match(self, pattern, syms):
        from symbol import Symbol
        if isinstance(pattern, Symbol):
            return {syms[syms.index(pattern)]: self}
        if isinstance(pattern, Rational):
            if self==pattern:
                return {}
        return None
   

class Constant(Basic):
    """Mathematical constant abstract class.
    
    Is the base class for constatns such as pi or e
    """
    
    def __init__(self):
        Basic.__init__(self, is_commutative = True)
    
    def __call__(self, precision=28):
        return self.evalf(precision)
       
    def eval(self):
        return self
 
    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        return self._mhash.value

    def diff(self,sym):
        return Rational(0)

    def __mod__(self, a):
        raise NotImplementedError

    def __rmod__(self, a):
            raise NotImplementedError

class ImaginaryUnit(Constant):
    """Imaginary unit "i"."""

    def __init__(self):
        Basic.__init__(self, 
                       is_real = False, 
                       is_commutative = True, 
                       )

    def __str__(self):
        return "I"
    
    def evalf(self):
        """Evaluate to a float. By convention, will return 0, 
        which means that evalf() of a complex number will mean 
        the projection of the complex plane to the real line. 
        For example:
        >>> (1-2*I).evalf()
        1.00
        >>> (-2+1*I).evalf()
        (-2.0)
        """
        return Rational(0)

    def evalc(self):
        return self

I = ImaginaryUnit()

class ConstPi(Constant):
    """
    
    Usage
    ===== 
           pi -> Returns the mathematical constant pi 
           pi() -> Returns a numerical aproximation for pi
           
    Notes
    =====
        Can have an option precision (integer) for the number of digits 
        that will be returned. Default is set to 28
       
        pi() is a shortcut for pi.evalf()
    
    Further examples
    ================
        >>> pi
        pi

        >>> pi()
        3.14159265358979323846264338

        >>> pi(precision=200)
        3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303820

    """
    
    def __init__(self):
        Basic.__init__(self,
                       is_commutative = True, 
                       is_real = True, 
                       )
    
    def evalf(self, precision=28):
        """Compute Pi to the current precision.

        >>> print pi.evalf()
        3.14159265358979323846264338
        
        """
        _pi_str = '3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068'
        if precision <= 100:
            # cache for small precision 
            return Real(_pi_str[:precision])
        #for arbitrary precision, we use series
        # FIXME: better algorithms are known
        decimal.getcontext().prec = precision + 2  # extra digits for intermediate steps
        three = decimal.Decimal(3)      # substitute "three=3.0" for regular floats
        lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
        while s != lasts:
            lasts = s
            n, na = n+na, na+8
            d, da = d+da, da+32
            t = (t * n) / d
            s += t
        decimal.getcontext().prec -= 2
        return Real(+s)               # unary plus applies the new precision
        # this was a recipe taken from http://docs.python.org/lib/decimal-recipes.html
        # don't know how fiable it is


    def __str__(self):
        return "pi"

pi=ConstPi()

def isnumber(x):
    """DEPRECATED"""
    #don't use this function. Use x.is_number instead
    #everything in sympy should be subclasses of Basic anyway.

    #if you need the testig for int, float, etc., just do it locally in your
    #class, or even better, call Basic.sympify(x).is_number.
    #so that all the code which converts from python to sympy is localised in 
    #sympify
    from numbers import Number
    from basic import Basic
    from decimal import Decimal
    if isinstance(x, (Number, int, float, long, Decimal)):
        return True
    assert isinstance(x, Basic)
    return x.is_number

def sign(x):
    """Return the sign of x, that is, 
    1 if x is positive, 0 if x == 0 and -1 if x is negative
    """
    if x < 0: return -1
    elif x==0: return 0
    else: return 1
