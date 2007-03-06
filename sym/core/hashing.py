"""This module provides hashing.

usage:
    h=mhash()
    h.addint(15)
    h.addstr("hey man!")
    h.addfloat(45.3)

    print h.value

mhash.value contains the hash value which is generally a long int, that
depends on the order of objects added. 
"""

def mhash():
    #return mmd5()
    #return mbernstein()
    return mpython()

class mpython(object):
    """
    An algorithm used in python, the result is a 4 byte hash value.

    It is not very strong.
    """
    i2p31=2**31
    i2p32=2**32

    def __init__(self):
        self.value=0x3456

    def trimlong(self,x):
        """Trims an unsigned integer"""
        x=x & 0xFFFFFFFFL
        if x>=self.i2p31: x-=self.i2p32
        return int(x)

    def add(self,item):
        """Adds any integer to the hash value."""
        self.value = self.trimlong(1000003*self.value)^item

    def addstr(self,x):
        #hash only last 5 letters, because the python algorithm is not
        #very strong
        self.add(hash(x[-5:]))

    def addint(self,x):
        """We are not using self.add directly for some reason I forget
        about."""
        self.add(x+3)

    def addfloat(self,x):
        self.add(hash(x)+3)

class mbernstein(object):
    """
    An algorithm produced by Professor Daniel J. Bernstein

    http://www.partow.net/programming/hashfunctions/index.html

    The result is generally a long int hash value. Unfortunately, this 
    algorithm is not strong enough, it doesn't pass 
    the test test_bug5h2() in test_hashing.py.
    """
    def __init__(self):
        self.value=5381

    def add(self,item):
        """Adds any integer to the hash value."""
        self.value = ((self.value << 5) + self.value) + item

    def addstr(self,x):
        self.add(hash(x))

    def addint(self,x):
        """We are not using self.add directly for some reason I forget
        about."""
        self.add(x+3)

    def addfloat(self,x):
        self.add(hash(x)+3)

class mmd5(object):
    """
    MD5 algorithm, this works fine (so far :).
    """

    def __init__(self):
        import md5
        self.md5=md5.new("heja")
        self.value=int(self.md5.hexdigest(),16)

    def add(self,item):
        """Adds any integer to the hash value."""
        self.md5.update(str(item))
        self.value = int(self.md5.hexdigest(),16)

    def addstr(self,x):
        self.add(x)

    def addint(self,x):
        """We are not using self.add directly for some reason I forget
        about."""
        self.add(x+3)

    def addfloat(self,x):
        self.add(hash(x)+3)
