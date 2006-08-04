#copy of the maple code, http://www.cybertester.com/data/gruntz.pdf,
#page 131 (appendix A)

#need to fix the failing assert at the bottom, in the mrv
import sym as s

def limitinf(e,x):
    "computes a limit x->inf of e(x), returns 1..inf,-1..-inf,0...in between"
    return mapping(e,[
        (s.exp(x),1),
        (-s.exp(x),-1),
        (x,1),
        (x+1/x,1),
        (x+s.exp(-x),1),
        (x+s.exp(-s.exp(x)),1),
        (-x,-1),
        ])

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
    assert isinstance(b,s.number)
    return b.sign()

def mapping(b,m):
    for x in m:
        if b.eval().isequal(x[0].eval()): return x[1]
    raise "%s not found in the mapping"%b

def has(e,x):
    return not e.diff(x).isequal(s.rational(0))
def Simplify(e):
    return e.expand().eval()
def MrvLeadTerm(e,x):
    e=e.eval()
    #print "MrvLeadTerm:",e
    return mapping(e,(
        ((x+s.exp(-x))/(-x),(e,s.rational(1),s.rational(-1))),
        ((x+s.exp(-s.exp(x)))/x,(e,s.rational(1),s.rational(1))),
        ((-s.exp(x))/x,(e,s.rational(1),s.rational(1))),
        ))
    if not has(e,x): return (e,s.rational(1),s.rational(0))
    return Mrv(e,x)
def Mrv(e,x):
    "Returns the most rapidly varying (mrv) subexpressions of 'e'"
    if not has(e,x): return []
    elif e.isequal(x): return [x]
    elif isinstance(e,s.mul): 
        a,b=e.getab()
        return Max(Mrv(a,x),Mrv(b,x),x)
    elif isinstance(e,s.add): 
        a,b=e.getab()
        return Max(Mrv(a,x),Mrv(b,x),x)
    elif isinstance(e,s.pow) and isinstance(e.b,s.number):
        return Mrv(e.a,x)
    elif isinstance(e,s.exp): 
        if limitinf(e.arg,x) in [-1,1]:
            return Max([e],Mrv(e.arg,x),x)
        else:
            return Mrv([e.arg])
    raise "unimplemented in Mrv: %s"%e
def Max(f,g,x):
    """Computes the maximum of two sets of expressions f and g, which 
    are in the same comparability class, i.e. max() compares (two elements of)
    f and g and returns the set, which is in the higher comparability class
    of the union of both, if they have the same order of variation.

    page 40 (47).
    """
    #print "max:",f,g
    if f==[]: return g
    elif g==[]: return f
    elif intersect(f,g): return union(f,g)
    elif member(x,f): return g
    elif member(x,g): return f
    else:
        c=Compare(f[0],g[0],x)
        if c==">": return f
        elif c=="<": return g
        else: return union(f,g)
    raise "max error",f,g

def Compare(a,b,x):
    #print "compare:",a,b
    c=MrvLeadTerm(s.ln(a)/s.ln(b),x)
    d=signum(0,c[2],0)
    if d==-1: return ">"
    elif d==1: return "<"
    elif d==0: return "="
    raise "compare error"

def intersect(a,b):
    for x in a:
        if member(x,b): return True
    return False

def member(x,a):
    for y in a:
        if x.isequal(y): return True
    return False

def union(a,b):
    z=a[:]
    for x in b:
        if not member(x,a):
            z.append(x)
    return z

def eq(a,b):
    print a,b
    assert len(a)==len(b)
    for x,y in zip(a,b):
        if not x==y:
            print "not equal:",x,y
            assert False

x=s.symbol("y")
#e=(s.exp(y)-1)/y
#print limit(e,y,s.rational(0))
eq(Mrv(s.exp(x+1/x),x),[s.exp(x+1/x)])
eq(Mrv(s.exp(x+s.exp(-s.exp(x))),x),[s.exp(-s.exp(x))] )
eq(Mrv(s.exp(x+s.exp(-x)),x),[s.exp(x+s.exp(-x)),s.exp(x)])
