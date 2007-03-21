"""
This module provides hashing.

usage:
    >>> h = mhash()
    >>> h.addint(15)
    >>> h.addstr("hey man!")
    >>> h.addfloat(45.3)
    >>> print h.value
    1550928569
    
mhash.value contains the hash value which is generally a long int, that
depends on the order of objects added. 
"""

def mhash():
    #return Mmd5()
    #return Mbernstein()
    return Mpython()

class HashAlgorithm(object):
    """Abstract class for varios implementatios of a hash algorithm"""

    def addint(self, num):
        """add an integer to the hash"""
        self.add(num)

    def addfloat(self, num):
        """Add a foat to the hash"""
        self.add(hash(num))
        
    def add(self, other):
        """Add an object to the hash"""
        raise NotImplementedError

class Mpython(HashAlgorithm):
    """
    An algorithm used in python, the result is a 4 byte hash value.

    It is not very strong.
    """
    i2p31 = 2**31
    i2p32 = 2**32

    def __init__(self):
        self.value = 0x3456
        HashAlgorithm.__init__(self)

    def trimlong(self, x):
        """Trims an unsigned integer"""
        x = x & 0xFFFFFFFFL
        if x >= self.i2p31: x -= self.i2p32
        return int(x)

    def add(self, item):
        """Adds any integer to the hash value."""
        self.value = self.trimlong(1000003*self.value)^item

    def addint(self, num):
        """The hash function needs this:"""
        self.add(num+3)

    def addstr(self, x):
        """Adds a string to the has value"""
        #hash only last 5 letters, because the python algorithm is not
        #very strong
        self.add(hash(x[-5:]))

class Mbernstein(HashAlgorithm):
    """
    An algorithm produced by Professor Daniel J. Bernstein

    http://www.partow.net/programming/hashfunctions/index.html

    The result is generally a long int hash value. Unfortunately, this 
    algorithm is not strong enough, it doesn't pass 
    the test test_bug5h2() in test_hashing.py.
    """
    def __init__(self):
        self.value = 5381
        HashAlgorithm.__init__(self)

    def add(self, item):
        """Adds any integer to the hash value."""
        self.value = ((self.value << 5) + self.value) + item

    def addstr(self, x):
        self.add(hash(x))


class Mmd5(HashAlgorithm):
    """
    MD5 algorithm, this works fine (so far :).
    """

    def __init__(self):
        import md5
        self.md5 = md5.new("heja")
        self.value=int(self.md5.hexdigest(),16)
        HashAlgorithm.__init__(self)

    def add(self, item):
        """Adds any integer to the hash value."""
        self.md5.update(str(item))
        self.value = int(self.md5.hexdigest(),16)

    def addstr(self, st):
        """Adds an string to the hash value"""
        self.add(st)

