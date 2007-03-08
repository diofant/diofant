import hashing
from basic import Basic
from numbers import Rational

dummycount=0

class Symbol(Basic):

    def __init__(self,name,dummy=False):
        """if dummy==True, then this Symbol is totally unique, i.e.:
            Symbol("x")==Symbol("x")
        but:
            Symbol("x",True)!=Symbol("x",True)

        """
        Basic.__init__(self)
        self.name=name
        self.dummy=dummy
        if dummy:
            global dummycount
            dummycount+=1
            self.dummycount=dummycount

    def __str__(self):
        return str(self.name)

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
