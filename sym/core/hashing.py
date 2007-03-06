class mhash(object):
    """Provides hashing.
    
    usage:
        h=mhash()
        h.addint(15)
        h.addstr("hey man!")
        h.addfloat(45.3)

        print h.value

    mhash.value contains the hash value, which depends on the order of
    objects added. The value can be platform independent or dependent as set by
    the "platform_independent" variable.

    algorithm = 0 .... an algorithm used in python, the result is a 4 byte
                        hash value
              = 1 .... an algorithm produced by Professor Daniel J. Bernstein
                see: http://www.partow.net/programming/hashfunctions/index.html
                        the result is generally a long int hash value
    """

    i2p31=2**31
    i2p32=2**32
    platform_independent=False
    #platform_independent=True
    algorithm=0

    def __init__(self):
        if self.algorithm==0:
            self.value=0x3456
        elif self.algorithm==1:
            self.value=5381

    def trimlong(self,x):
        """Trims an unsigned integer"""
        x=x & 0xFFFFFFFFL
        if x>=mhash.i2p31: x-=mhash.i2p32
        return int(x)

    def add(self,item):
        """Adds any integer to the hash value."""
        if self.algorithm==0:
            self.value = self.trimlong(1000003*self.value)^item
        elif self.algorithm==1:
            self.value = ((self.value << 5) + self.value) + item

    def addstr(self,x):
        if self.platform_independent:
            for l in x:
                self.addint(ord(l))
        else: #faster
            if self.algorithm==0:
                #hash only last 5 letters, because the python algorithm is not
                #very strong
                self.add(hash(x[-5:]))
            elif self.algorithm==1:
                self.add(hash(x))

    def addint(self,x):
        """We are not using self.add directly for some reason I forget
        about."""
        self.add(x+3)

    def addfloat(self,x):
        if self.platform_independent:
            raise "I don't know how to hash floats in platform independent way"
        else:
            self.add(hash(x)+3)
