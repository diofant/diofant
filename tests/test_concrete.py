import sys
sys.path.append(".")

from sympy.modules.concrete.primes import *

def test_structures():
    bits = BitField(11027)

    assert len(bits) == 11027

    for i in [0, 17, 1024, 1027, 11026]:
        bits[i] = True

    assert bits[0] == True
    assert bits[1] == False

    assert bits[17] == True

    assert bits[1023] == False
    assert bits[1024] == True
    assert bits[1025] == False
    assert bits[1026] == False
    assert bits[1027] == True
    assert bits[1028] == False

    assert bits[11025] == False
    assert bits[11026] == True

    assert len([ bit for bit in bits if bit ]) == 5

def test_primes():
    g = PrimeGenerator(123)

    assert g.primes(2, 20) == [2, 3, 5, 7, 11, 13, 17, 19]
    assert g.primes(7, 30) == [7, 11, 13, 17, 19, 23, 29]

    assert g.is_prime(-29) == is_prime(-29) == False

    assert g.is_prime(0) == is_prime(0) == False
    assert g.is_prime(1) == is_prime(1) == False
    assert g.is_prime(2) == is_prime(2) == True

    assert g.is_prime(17) == is_prime(17) == True
    assert g.is_prime(123) == is_prime(123) == False

    assert is_prime(32582657) == True
    assert is_prime(2147483647) == True

