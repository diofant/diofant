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
