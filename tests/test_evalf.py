import sys
sys.path.append(".")

import sym as s

def eq(a,b):
    return abs(a-b)<0.0001

def testeval():
    e=s.ln(3)/s.ln(2)-1
    assert eq(e.evalf(),0.58496)
