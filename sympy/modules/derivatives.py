
from sympy.core import Basic, hashing

#currently, the Derivative class is in core.functions
class xDerivative(Basic):

    def __init__(self,f,x):
        Basic.__init__(self)
        self.f=self.sympify(f)
        self.x=self.sympify(x)

    def doit(self):
        return self.f.diff(self.x)

    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.f.hash())
        self.mhash.addint(self.x.hash())
        return self.mhash.value

def derivate(f, x, times):
    pass
