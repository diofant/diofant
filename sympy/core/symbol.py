
from basic import Basic
from numbers import Rational
#from prettyprint import StringPict

dummycount = 0

class Symbol(Basic):
    """
    Assumptions::
       is_real = True
       is_commutative = True

    You can override the default assumptions in the constructor::
       >>> A = Symbol('A', is_commutative = False)
       >>> B = Symbol('B', is_commutative = False)
       >>> A*B != B*A
       True
       >>> A*B*2 == 2*A*B # multiplication by scalars is commutative
       True
    """
    
    mathml_tag = "ci"
    
    dummy_num = 0

    def __init__(self, name, *args, **kwargs):
        """if is_dummy==True, then this Symbol is totally unique, i.e.::
        
        >>> Symbol("x") == Symbol("x")
        True
        
        but with the dummy variable ::
        
        >>> Symbol("x", is_dummy = True) == Symbol("x", is_dummy = True)
        False

        """
        
        self._assumptions = {
                         'is_commutative' : True,
                         }
        
        for k in kwargs.keys():
            self._assumptions[k] = kwargs[k]
        
        Basic.__init__(self, **self._assumptions)
        self.name = name
        if self.is_dummy:
            global dummycount
            self.dummy_num = dummycount
            dummycount += 1

    def __str__(self):
        if not self.is_dummy:
            return str(self.name)
        else:
            # if x is dummy
            return str(self.name + '__' + str(self.dummy_num))
    
    def __getitem__(self, iter):
        return (self,)[iter]

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

    def match(self, pattern, syms):
        if self == pattern:
            return {}
        if len(syms) == 1:
            if pattern == syms[0]:
                return {syms[0]: self}
            if self == pattern:
                return {}
        if isinstance(pattern, Symbol):
            try:
                return {syms[syms.index(pattern)]: self}
            except ValueError:
                pass
        from addmul import Mul
        if isinstance(pattern, Mul):
            return Mul(Rational(1),self,evaluate = False).match(pattern,syms)
        return None

class Order(Basic):

    def __init__(self, f):
        """O(f) at the point x = 0"""
        Basic.__init__(self)
        self._args = [self.sympify(f)]

    def eval(self):
        from addmul import Mul, Add
        from numbers import Real, Rational
        f = self[0]
        if isinstance(f, Mul):
            if isinstance(f[0], (Real, Rational)):
                assert len(f[:]) == 2
                return Order(f[1])
        if isinstance(f, Add):
            if isinstance(f[0], (Real, Rational)):
                assert len(f[:]) == 2
                return Order(f[1])
        return self

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
