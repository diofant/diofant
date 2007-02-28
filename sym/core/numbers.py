import hashing
from basic import Basic
import utils 

class Number(Basic):
    
    def __init__(self):
        Basic.__init__(self)
        self.evaluated = True
        
    def addnumber(self,a):
        if isinstance(a,Real):
            return Real(self.evalf()+a.evalf())
        return self.add(a)
    
    def mulnumber(self,a):
        if isinstance(a,Real):
            return Real(self.evalf()*a.evalf())
        return self.mul(a)
    
    def pownumber(self,a):
        if isinstance(a,Real):
            return Real(self.evalf()**a.evalf())
        return self.pow(a)
    
    def diff(self,sym):
        return Rational(0)

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
        return "inf"
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        return self.mhash.value
    def sign(self):
        return self._sign

infty=Infinity()

class Real(Number):
    def __init__(self,num):
        Number.__init__(self)
        if isinstance(num,str):
            num=float(num)
        assert isinstance(num,float) or isinstance(num,int)
        self.num=num
        
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addfloat(self.num)
        return self.mhash.value
        
    def __str__(self):
        if self.num < 0:
            f="(%r)"
        else:
            f="%r"
        return f%(self.num)
        
    def add(self,a):
        return Real(self.num+a.evalf())
        
    def mul(self,a):
        return Real(self.num*a.evalf())
        
    def pow(self,a):
        return Real(self.num**a.evalf())
        
    def iszero(self):
        return False
        
    def isone(self):
        return False
        
    def isinteger(self):
        return False
        
    def evalf(self):
        #evalf() should return either a float or an exception
        return self.num

class Rational(Number):
    
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
        assert (isinstance(p,int) or isinstance(p,long)) and \
                (isinstance(q,int) or isinstance(q,long))
        assert q!=0
        s = utils.sign(p)*utils.sign(q)
        p = abs(p)
        q = abs(q)
        c = self.gcd(p,q)
        self.p = p/c*s
        self.q = q/c
        
    def __lt__(self,a):
        import basic
        a=basic.c(a)
        assert isinstance(a,Rational)
        return self.p * a.q < self.q * a.p
        
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
        while b!=0:
            c = a % b
            a = b
            b = c
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
            
    def mul(self,a):
        return Rational(self.p*a.p,self.q*a.q)
        
    def add(self,a):
        return Rational(self.p*a.q+self.q*a.p,self.q*a.q)
        
    def pow(self,a):
        assert a.q==1
        if a.p > 0:
            return Rational(self.p**a.p,self.q**a.p)
        else:
            return Rational(self.q**(-a.p),self.p**(-a.p))
            
    def iszero(self):
        return self.p==0 
        
    def isone(self):
        return self.p==1 and self.q==1
        
    def isminusone(self):
        return self.p==-1 and self.q==1
        
    def isinteger(self):
        return self.q==1
        
    def getinteger(self):
        assert self.isinteger()
        return self.p
        
    def evalf(self):
        return float(self.p)/self.q
        
    def diff(self,sym):
        return Rational(0)
