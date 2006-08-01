class mhash(object):
    """Provides hashing.
    
    usage:
        h=mhash()
        h.addint(15)
        h.addstr("hey man!")
        h.addfloat(45.3)

        print h.value

    mhash.value contains 4 byte hash value. The number is platform
    dependent due to the usage of the python build in (platform dependent)
    hash() function.

    hash value depends on the order of objects added.
    """

    i2p31=2**31
    i2p32=2**32

    def __init__(self):
        self.value=0x3456

    def trimlong(self,x):
        x=x & 0xFFFFFFFFL
        if x>=mhash.i2p31: x-=mhash.i2p32
        return int(x)

    def add(self,item):
        self.value=self.trimlong(1000003*self.value)^item

    def addstr(self,x):
        self.add(hash(x))

    def addint(self,x):
        if x==-1: self.add(-1)
        else: self.add(hash(x))

    def addfloat(self,x):
        if x==-1: self.add(-1)
        else: self.add(hash(x))
