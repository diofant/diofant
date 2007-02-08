import sys
sys.path.append("..")

import sym
a=sym.symbol('a')
b=sym.symbol('b')
e=sym.ln((a+b)**5)
print e
print e.eval()
e=sym.exp(e)
print e
print e.eval()
e=sym.ln(sym.exp((a+b)**5))
print e
print e.eval()
