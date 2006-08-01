#copy of the maple code, http://www.cybertester.com/data/gruntz.pdf,
#page 131 (appendix A)
import sym as s

def limit(e,z,z0):
    print "e:",e
    print "z:",z
    print "z0:",z0
    x=s.symbol("x")
    e0=e.subs(z,z0+1/x)
    e1=MrvLeadTerm(Simplify(e0),x)
    if e1[2].isequal(s.rational(0)): r=e1[0]#r=e1[0].expand()
    elif signum(0,e1[2],0)==0: r=e1[0].expand()
    elif signum(0,e1[2],0)==1: r=0
    elif signum(0,e1[2],0)==-1: r=Sign(e1[0])*infinity
    else: raise "cannot determine the sign of %s"%(e1[2])
    return r

def signum(a,b,c):
    return mapping(b,[(s.rational(3),1)])

def mapping(b,m):
    for x in m:
        if b.eval().isequal(x[0].eval()): return x[1]
    raise "%s not found in the mapping"%b

def has(e,x):
    return not e.diff(x).isequal(s.rational(0))
def Simplify(e):
    return e.expand().eval()
def MrvLeadTerm(e,x):
    if not has(e,x): return (e,s.rational(1),s.rational(0))
    return Mrv(e,x)
def Mrv(e,x):
    if not has(e,x): return []
    elif isinstance(e,s.add):
    return mapping(e,[(x),(1,2,s.rational(3)))])

y=s.symbol("y")
e=(s.exp(y)-1)/y
print limit(e,y,s.rational(0))
