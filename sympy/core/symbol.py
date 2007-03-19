import hashing
from basic import Basic
from numbers import Rational
#from prettyprint import StringPict

dummycount=0

class Symbol(Basic):
    """
    Assumptions: 
       is_real = True
       is_commutative = True
       
   You can override the default assumptions in the constructor:
   >>> A = Symbol('A', is_commutative = False)
   A
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

    def print_tex(self):
        return str(self.name_tex)

    def print_pretty(self):
		return StringPict(self.print_sympy())
        
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addstr(self.name)
        if self.dummy:
            self.mhash.addint(self.dummycount)
        return self.mhash.value

    def diff(self,sym):
        if self.isequal(sym):
            return Rational(1)
        else:
            return Rational(0)

    def evalc(self):
        if self.is_real:
            return self
        raise NotImplementedError

class NCSymbol(Symbol):

    def commutative(self):
        return False

    def diff(self,sym):
        raise NotImplementedError("Doesn't have a meaning.")
