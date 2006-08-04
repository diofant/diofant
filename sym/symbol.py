import hashing
from basic import basic
from numbers import rational

class symbol(basic):

    def __init__(self,name):
        basic.__init__(self)
        self.name=name

    def __str__(self):
        return str(self.name)

    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addstr(self.name)
        return self.mhash.value

    def diff(self,sym):
        if self.isequal(sym):
            return rational(1)
        else:
            return rational(0)
