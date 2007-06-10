
from sympy.core.basic import Basic

from sympy.modules.concrete.structures import BitField

class PrimeGenerator(object):
    """Container for prime numbers. Uses sieve of Eratosthenes
       algorithm over bitmaps to compute primes, which makes it
       to consume little memory. This container is appropriate
       to cache primes. If primality of few, especially big,
       integers is needed, better use test algorithms.

       TODO: Probably it would be better to use sieve of Atkin.
    """

    def __init__(self, n):
        self.storage = self.sieve(n)

    @staticmethod
    def sieve(n):
        """Generates bitmap in which bits, in 0..n range, representing
           prime numbers are set and those representing composite numbers
           are made unset. This functions is a static method so it can
           be used alone or a it can be ised in pair with primes(p, q)
           or is_prime(p) functions if an instance of PrimeGenerator
           was created.

           >>> primes = PrimeGenerator.sieve(10)

           >>> primes[4]
           False
           >>> primes[7]
           True

        """

        primes = BitField(n+1, True)

        # sorry, zero and one aren't primes
        primes[0] = primes[1] = False

        for p in range(2, n+1):
            if primes[p]:
                for i in range(2*p, n+1, p):
                    primes[i] = False

        return primes

    def primes(self, p, q):
        """Returns a list containg prime numbers from p..q range so
           effectively it will return indexes of active bits stored
           in internal primes bitmap. If the bitmap is too small,
           this function will expand it automatically.

           >>> gen = PrimeGenerator(20)

           >>> gen.primes(2, 10)
           [2, 3, 5, 7]
           >>> gen.primes(7, 30)
           [7, 11, 13, 17, 19, 23, 29]

        """

        if p < 0 or q < p:
            raise ValueError("Invalid bounds")
        else:
            # not enough primes, generate more
            if q >= len(self.storage):
                self.storage = self.sieve(q)

            primes = []

            for i in range(p, q+1):
                if self.storage[i]:
                    primes.append(i)

            return primes

    def is_prime(self, p):
        if p <= 1:
            return False
        else:
            # not enough primes, generate more
            if p >= len(self.storage):
                self.storage = self.sieve(p)

            return self.storage[p]

def is_prime(p):
    """Test primality of a number using simple trial division
       algorithm. This is just playground for real probabilistic
       and deterministic tests.

       TODO: Implement AKS algorithm.
    """

    if not isinstance(p, int):
        p = Basic.sympify(p)

        if p.is_integer:
            p = int(p)
        else:
            raise TypeError("Needs integer argument")

    if p == 2:
        return True
    elif p < 2 or p & 1 == 0:
        return False
    else:
        t, s = 2, 4

        while s <= p:
            if p % t == 0:
                return False
            else:
                t, s = t+1, s+2*t-1

        return True

