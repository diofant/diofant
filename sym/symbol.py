import hashing
from basic import basic
from numbers import rational

dummycount=0

class symbol(basic):

    def __init__(self,name,dummy=False):
        """if dummy==True, then this symbol is totally unique, i.e.:
            symbol("x")==symbol("x")
        but:
            symbol("x",True)!=symbol("x",True)

        """
        basic.__init__(self)
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
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addstr(self.name)
        if self.dummy:
            self.mhash.addint(self.dummycount)
        return self.mhash.value

    def diff(self,sym):
        if self.isequal(sym):
            return rational(1)
        else:
            return rational(0)
