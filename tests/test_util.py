import sys
sys.path.append(".")

import sym as s

def test_sign():
    assert s.utils.sign(-2)  == s.utils.sign(s.Rational(-2)) 
    assert s.utils.sign(s.Rational(-2)) == -1
#    assert s.utils.sign(s.Real(0)*1) == 0 # fails because of issue 21
    assert s.utils.sign(s.Rational(-1)*s.Rational(-1)) == 1


