#copy of the maple code, http://www.cybertester.com/data/gruntz.pdf,
#page 131 (appendix A)

"""
f > g : f is greater than any power of g

f < g : f is lower than any power of g

f ~ g : f and g are both bounded above and below by suitable integral powers of
the other

1 ... the lowest comparability class

f > g .... f is more rapidly varying than g = f goes to infinity faster than g
"""

#need to fix the failing assert at the bottom, in the mrv
import sym as s

def intersect(a,b):
    for x in a:
        if member(x,b): return True
    return False

def member(x,a):
    for y in a:
        if x == y: return True
    return False

def union(a,b):
    z=a[:]
    for x in b:
        if not member(x,a):
            z.append(x)
    return z

def limitinf_manual(e,x):
    "computes a limit x->inf of e(x), returns 1..inf,-1..-inf,0...in between"
    return mapping(e,[
        (s.exp(x),1),
        (-s.exp(x),-1),
        (x,1),
        (-x**2,-1),
        (x+1/x,1),
        (x+s.exp(-x),1),
        (x+s.exp(-s.exp(x)),1),
        (-x,-1),
        (1/x,0),
        (s.exp(-x**2)+x,1),
        (-s.exp(-x)+x**(-1),0),
        ])

def leadterm(series,x):
    def domul(x):
        if len(x)>1:
            return s.mul(x)
        return x[0]
    def extract(t,x):
        if not has(t,x):
            return t,s.rational(0)
        if isinstance(t,s.pow):
            return  s.rational(1),  t.b
        elif isinstance(t,s.symbol):
            return  s.rational(1),  s.rational(1)
        assert isinstance(t,s.mul)
        for i,a in enumerate(t.args):
            if has(a,x):
                if isinstance(a,s.pow):
                    return  domul(t.args[:i]+t.args[i+1:]),  a.b
                elif isinstance(a,s.symbol):
                    return  domul(t.args[:i]+t.args[i+1:]),  s.rational(1)
                assert False
        return t,s.rational(0)
    assert isinstance(series,s.add)
    lowest=(0,(s.rational(10)**10).eval())
    for t in series.args:
        t2=extract(t,x)
        if t2[1]<lowest[1]:
            lowest=t2
    return lowest

def limit(e,z,z0):
    """Currently only limit z->z0+"""
    x=s.symbol("x")
    e0=e.subs(z,z0+1/x)
    return limitinf(e0,x)

def limitinf(e,x):
    """Limit e(x) for x-> infty"""
    e1=mrvleadterm(e.eval(),x)
    if e1[2].isequal(s.rational(0)): r=e1[0]
    elif signum(e1[2])==0: r=e1[0]
    elif signum(e1[2])==1: r=0
    elif signum(e1[2])==-1: r=Sign(e1[0])*infinity
    else: raise "cannot determine the sign of %s"%(e1[2])

    if has(r,x):
        return limitinf(r,x)
    else:
        return r

def signum(a):
    """Returns a sign of an expression at x->infinity"""
    assert isinstance(a,s.number)
    return a.sign()

def mapping(b,m):
    for x in m:
        if b==x[0]: return x[1]
    raise "%s not found in the mapping"%b.eval()

def has(e,x):
    return not e.diff(x).isequal(s.rational(0))

def mrvleadterm(e,x):
    """Returns (coeff,mrv_var,exponent) for e
    
    
    coeff*mrv_var^exponent ~ "e" for x->infinity
    """
    e=e.eval()
    if not has(e,x): return (e,s.rational(1),s.rational(0))
    Omega=mrv(e,x)
    assert len(Omega)==1
    wexpr=Omega[0]
    w=s.symbol("w")
    if wexpr==x:
        f2=e.subs(wexpr,1/w)
    else:
        f2=e.subs(wexpr,w)
    ser=f2.series(w,3)
    lterm=leadterm(ser,w)
    return lterm[0],wexpr,lterm[1]

def mrv(e,x):
    "Returns the most rapidly varying (mrv) subexpressions of 'e'"
    if not has(e,x): return []
    elif e.isequal(x): return [x]
    elif isinstance(e,s.mul): 
        a,b=e.getab()
        return max(mrv(a,x),mrv(b,x),x)
    elif isinstance(e,s.add): 
        a,b=e.getab()
        return max(mrv(a,x),mrv(b,x),x)
    elif isinstance(e,s.pow) and isinstance(e.b,s.number):
        return mrv(e.a,x)
    elif isinstance(e,s.exp): 
        if limitinf_manual(e.arg,x) in [-1,1]:
            return max([e],mrv(e.arg,x),x)
        else:
            return mrv(e.arg,x)
    raise "unimplemented in mrv: %s"%e

def max(f,g,x):
    """Computes the maximum of two sets of expressions f and g, which 
    are in the same comparability class, i.e. max() compares (two elements of)
    f and g and returns the set, which is in the higher comparability class
    of the union of both, if they have the same order of variation.

    page 40 (47).
    """
#    print "max:",f,g,member(x,g)
    if f==[]: return g
    elif g==[]: return f
    elif intersect(f,g): return union(f,g)
    elif member(x,f): return g
    elif member(x,g): return f
    else:
        c=compare(f[0],g[0],x)
        if c==">": return f
        elif c=="<": return g
        else: return union(f,g)
    raise "max error",f,g

def compare(a,b,x):
    """Returns "<" if a<b (at x=infinity), "=" for a==b, ">" for a>b"""
#    print "compare:",a,b
    if a==s.exp(-x) and b==x: return ">"
    elif b==s.exp(-x) and a==x: return "<"
    elif a==s.exp(x+1/x) and b==x: return ">"
    elif a==s.exp(x) and b==x**5: return ">"
    elif a==s.exp(x**2) and b==s.exp(x)**2: return ">"
    elif a==s.exp(x) and b==s.exp(x+s.exp(-x)): return "="
    elif b==s.exp(x) and a==s.exp(x+s.exp(-x)): return "="
    elif a==s.exp(s.exp(x)) and b==s.exp(x+s.exp(-s.exp(x))): return ">"
    elif b==s.exp(-x) and a==s.exp(x+s.exp(-x)): return "="
    elif a==s.exp(-s.exp(x)) and b==s.exp(x): return ">"
    elif a==s.exp(s.exp(-s.exp(x))+x) and b==s.exp(-s.exp(x)): return "<"
    elif a==s.exp(s.exp(-x**2)+x) and b==s.exp(-x**2): return "<"
    return mapping(a-b,[])
    c=mrvleadterm(s.ln(a)/s.ln(b),x)
    d=signum(c[2])
    if d==-1: return ">"
    elif d==1: return "<"
    elif d==0: return "="
    raise "compare error"
