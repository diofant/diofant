import sys
sys.path.append("..")

import sym
e=sym.Rational(2)**50/sym.Rational(10)**50
print e
print e.eval()
