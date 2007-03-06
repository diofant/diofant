import hashing
from basic import Basic,AutomaticEvaluationType
import utils 

class Number(Basic):
    """Represents any kind of number in sympy.


    Floating point numbers are represented by the Real class.
    Integer numbers (of any size), together with rational numbers (again, there
    is no limit on their size) are represented by the Rational class. 

    If you want to represent for example 1+sqrt(2), then you need to do:

    Rational(1) + Rational(2)**( Rational(1)/2 )
    """

#    __metaclass__ = AutomaticEvaluationType
    
    def __init__(self):
        Basic.__init__(self)
        self.evaluated = True
        
    def __int__(self):
        raise NotImplementedError
        
    def __float__(self):
        return self.evalf()
    
    def __abs__(self):
        from functions import abs_
        return abs_(self)
    
    def diff(self,sym):
        return Rational(0)
    
    def evalf(self):
        return self.eval()
    
class Infinity(Number):
    """Infinity. Cannot be used in expressions like 1+infty.  
    Only as a Symbol, for example results of limits, integration limits etc.
    Can however be used in comparisons, like infty!=1, or infty!=x**3
    
    this class represents all kinds of infinity, i.e. both +-infty.
    """
    
    def __init__(self):
        Number.__init__(self)
        self._sign=1
        
    def __str__(self):
        return "Inf"
    
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        return self.mhash.value
    
    def sign(self):
        return self._sign

infty=Infinity()

class Real(Number):
    """Represents a floating point number.

    Currently, it supports python floats only.

    Usage:

    Real(3.5)   .... 3.5 (the 3.5 was converted from a python float)
    Real("3.5") .... 3.5 (currently, the 3.5 is also a python float,
            but in the future, we could use some other library)
    """
    
    def __init__(self,num):
        Number.__init__(self)
        if isinstance(num,str):
            num = float(num)
        assert isinstance(num,float) or isinstance(num,int)
        self.num = float(num)
        
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addfloat(self.num)
        return self.mhash.value
        
    def __str__(self):
        if self.num < 0:
            f = "(%r)"
        else:
            f = "%r"
        return f % (self.num)
    
    def __float__(self):
        return float(self.num)
        
    def __int__(self):
        return int(self.evalf())
    
    def __add__(self,a):
        if utils.isnumber(a):
            return Real(self.num + float(a))
        else:
            assert isinstance(a, Basic)
            from addmul import Add
            return Add(self, a)
        
    def __mul__(self,a):
        if utils.isnumber(a):
            return Real(self.num * float(a))
        else:
            assert isinstance(a, Basic)
            from addmul import Mul
            return Mul(self, a)
        
    def __pow__(self,a):
        if utils.isnumber(a):
            return Real(self.num ** float(a))
        else:
            assert isinstance(a, Basic)
            from power import Pow
            return Pow(self, a)
        
    def __rpow__(self, a):
        if utils.isnumber(a):
            return float(a) ** self.num
        else:
            assert isinstance(a, Basic)
            from power import Pow
            return Pow(a, self)
        
    def __lt__(self, a):
        return self.num < a
    
    def __gt__(self, a):
        return self.num > a
        
    def iszero(self):
        if self.num == 0:
            return True
        else: 
            return False
        
    def isone(self):
        if self.num == 1:
            return True
        else:
            return False
        
    def isinteger(self):
        return False
        
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
        Number.__init__(self)
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
        s = utils.sign(p)*utils.sign(q)
        p = abs(p)
        q = abs(q)
        c = self.gcd(p,q)
        self.p = p/c*s
        self.q = q/c
        
    def __lt__(self,a):
        """Compares two Rational numbers."""
        return self.evalf() < float(a)
        
    def sign(self):
        return utils.sign(self.p)*utils.sign(self.q)
        
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.p)
        self.mhash.addint(self.q)
        return self.mhash.value
        
    def gcd(self,a,b):
        """Primitive algorithm for a greatest common divisor of "a" and "b"."""
        while b!=0:
            c = a % b
            a,b=b,c
        return a
        
    def __str__(self):
        if self.q == 1:
            if self.p < 0:
                f = "(%d)"
            else:
                f = "%d"
            return f % (self.p)
        else:
            if self.p < 0:
                f = "(%d/%d)"
            else:
                f = "%d/%d"
            return f % (self.p,self.q)
            
    def __mul__(self,a):
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
        if isinstance(a, int):
            return Rational(self.p, self.q *a)
        return self * (a**Rational(-1))
        
    def __rdiv__(self, a):
        if isinstance(a, int):
            return Rational(self.q * a, self.p )
        return self * (a**Rational(-1))
    
    def __add__(self,a):
        if isinstance(a, Rational):
            return Rational(self.p*a.q+self.q*a.p,self.q*a.q)
        elif isinstance(a, int) or isinstance(a, long):
            return Rational(self.p + a*self.q, self.q)
        else:
            from addmul import Add
            return Add(self, a)
        
    def __pow__(self,a):
        """Returns the self to the power of "a"
        """
        from power import Pow, pole_error
    
        if utils.isnumber(a):
            if self.p == 0:
                if a < 0:
                    # 0 ** a = undefined, where a <= 0 
                    raise pole_error("pow::eval(): Division by 0.")
                elif a == 0:
                    return Rational(1)
                    #FIXME : mathematically wrong but needed for limits.py
                else:
                    # 0 ** a = 0, where a > 0
                    return Rational(0)
            elif isinstance(a, Rational):
                if a.q == 1:
                    if a.p > 0:
                        return Rational(self.p ** a.p, self.q ** a.p)
                    else:
                        return Rational(self.q**(-a.p),self.p**(-a.p))
        return Pow(self, self.sympify(a))
            
    def __rpow__(self, a):  
        """Returns "a" to the power of self
        """
        from power import Pow
        if self.p == 0:
            return Rational(1)
        elif a == 0:
            return Rational(0)
        if self.q != 1:
            #if self is an integer
            if hasattr(a, 'evalf'):
                return Pow(a.evalf(), self)
            else:
                return Pow(self.sympify(a), self)
        elif isinstance(a, Rational):
            if self.p > 0:
                return Rational(a.p ** self.p, a.q ** self.p)
            else:
                return Rational(a.q ** (-self.p), a.p ** (-self.p))
        elif isinstance(a, int):
            if self.p > 0:
                return Rational(a ** self.p)
            else:
                return Rational(1, a ** -(self.p))
        return Pow(a, self )
    
    def iszero(self):
        return self.p == 0 
        
    def isone(self):
        return self.p == 1 and self.q == 1
        
    def isminusone(self):
        return self.p == -1 and self.q == 1
        
    def isinteger(self):
        return self.q == 1
        
    def getinteger(self):
        assert self.isinteger()
        return self.p
        
    def evalf(self):
        return float(self.p)/self.q
        
    def diff(self,sym):
        return Rational(0)
    

class ImaginaryUnit(Basic):
    """Imaginary unit "i"."""

    def __str__(self):
        return "i"

    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        return self.mhash.value

I=ImaginaryUnit()

class Constant(Basic):
    """Mathematical constant abstract class."""

    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        return self.mhash.value

class ConstPi(Constant):

    def __str__(self):
        return "pi"

pi=ConstPi()
