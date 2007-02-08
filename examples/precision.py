import sys
sys.path.append("..")

import sym
e=sym.rational(2)**50/sym.rational(10)**50
print e
print e.eval()
