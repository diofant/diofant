class mhash(object):
    """Provides hashing.
    
    usage:
        h=mhash()
        h.addint(15)
        h.addstr("hey man!")
        h.addfloat(45.3)

        print h.value

    mhash.value contains 4 byte hash value, which depends on the order of
    objects added. The value can be platform independent or dependent as set by
    the "platform_independent" variable.
    """

    i2p31=2**31
    i2p32=2**32
    platform_independent=False
    #platform_independent=True

    def __init__(self):
        self.value=0x3456

    def trimlong(self,x):
        x=x & 0xFFFFFFFFL
        if x>=mhash.i2p31: x-=mhash.i2p32
        return int(x)

    def add(self,item):
        self.value=self.trimlong(1000003*self.value)^item

    def addstr(self,x):
        if self.platform_independent:
            for l in x:
                self.addint(ord(l))
        else: #faster
            #hash only last 5 letters
            self.add(hash(x[-5:]))
            #self.add(hash(x))

    def addint(self,x):
        self.add(x+3)

    def addfloat(self,x):
        if self.platform_independent:
            raise "I don't know how to hash floats in platform independent way"
        else:
            self.add(hash(x)+3)
