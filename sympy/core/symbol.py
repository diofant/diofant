import hashing
from basic import Basic
from numbers import Rational
#from prettyprint import StringPict

dummycount=0

class Symbol(Basic):
    """
    Assumptions::
       is_real = True
       is_commutative = True

    You can override the default assumptions in the constructor::
       >>> A = Symbol('A', is_commutative = False)
    """
    
    mathml_tag = "ci"

    def __init__(self, name, dummy=False, *args, **kwargs):
        """if dummy==True, then this Symbol is totally unique, i.e.::
        
        >>> Symbol("x") == Symbol("x")
        True
        
        but with the dummy variable ::
        
        >>> Symbol("x", dummy = True) == Symbol("x", dummy = True)
        False

        """
        
        self._assumptions = {
                         'is_commutative' : True, 
                         }
        
        for k in kwargs.keys():
            self._assumptions[k] = kwargs[k]
        
        Basic.__init__(self, **self._assumptions)
        self.name = name
        self.dummy = dummy
        if dummy:
            global dummycount
            dummycount+=1
            self.dummycount=dummycount

    def __str__(self):
        return str(self.name)
    
    def __getitem__(self, iter):
        return (self,)[iter]

    def __len__(self):
        return 1

    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        self._mhash.addstr(self.name)
        if self.dummy:
            self._mhash.addint(self.dummycount)
        return self._mhash.value

    def diff(self,sym):
        if not self.is_commutative:
            raise NotImplementedError("Differentiation of non-commutative objects. " + \
                                      + "Doesn't have a meaning.")
        if self.isequal(sym):
            return Rational(1)
        else:
            return Rational(0)

    def evalc(self):
        if self.is_real:
            return self
        raise NotImplementedError

class NCSymbol(Symbol):

    @property
    def is_commutative(self):
        return False

    def diff(self,sym):
        raise NotImplementedError("Doesn't have a meaning.")

class Order(Basic):

    def __init__(self, f):
        """O(f) at the point x = 0"""
        Basic.__init__(self)
        self._args = [self.sympify(f)]

    def eval(self):
        from addmul import Mul, Add
        from numbers import Number
        f = self[0]
        if isinstance(f, Mul):
            if isinstance(f[0],Number):
                assert len(f) == 2
                return Order(f[1])
        if isinstance(f, Add):
            if isinstance(f[0],Number):
                assert len(f) == 2
                return Order(f[1])
        return self

    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        self._mhash.addint(self[0].hash())
        return self._mhash.value

    def __str__(self):
        return "O(%s)"%str(self[0])

    @staticmethod
    def muleval(x, y):
        if isinstance(x, Order) and isinstance(y, Order):
            return Order(x[0]*y[0])
        if isinstance(y, Order):
            return Order(x*y[0])
        return None

    @staticmethod
    def addeval(x, y):
        if isinstance(x, Order) and isinstance(y, Order):
            raise NotImplementedError()
        if isinstance(x, Order):
            return Order(x[0]+y)
        if isinstance(y, Order):
            return Order(x+y[0])
        return None
