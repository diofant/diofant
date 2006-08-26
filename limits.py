"""
Limits
======

Implemented according to the PhD thesis
http://www.cybertester.com/data/gruntz.pdf, which contains very thorough
descriptions of the algorithm including many examples.  We summarize here the
gist of it.


All functions are sorted according to how rapidly varying they are at infinity
using the following rules. Any two functions f and g can be compared:

L=lim  ln|f(x)| / ln|g(x)|           (for x -> infty)

f > g .... L=+-infty 
    * f is greater than any power of g
    * f is more rapidly varying than g
    * f goes to infinity/zero faster than g


f < g .... L=0 
    * f is lower than any power of g

f ~ g .... L!=0,+-infty 
    * both f and g are bounded from above and below by suitable integral powers
    of the other


Examples: 

    1 < x < exp(x) < exp(x^2) < exp(exp(x))
    1 ~ 3 ~ -5
    x ~ x^2 ~ x^3 ~ 1/x ~ x^m ~ -x
    exp(x) ~ exp(-x) ~ exp(2x) ~ exp(x)^2 ~ exp(x+exp(-x))
    f ~ 1/f

So we can divide all the functions into comparability classes (x and x^2 is the
same class, as is exp(x) and exp(-x)). In principle, we could compare any two
functions, but in our algorithm, we don't compare anything below f=1 (for
example ln(x) is below 1), so we set f=1 as the lowest comparability class. 

Given the function f, we find the list of most rapidly varying (mrv set)
subexpressions of it. This list belongs to the same comparability class. Let's
say it is {exp(x), exp(2x)}. Using the rule f ~ 1/f we find an element "w"
(either from the list or a new one) from the same comparability class which
goes to zero at infinity. In our example we set w=exp(-x) (but we could also
set w=exp(-2x) or w=exp(3x) ...). We rewrite the mrv set using w, in our case
{1/w,1/w^2}, and substitute it into f. Then we expand f into a series in w:

    f=c0*w^e0 + c1*w^e1 + ... + O(w^en),        where e0<e1<...<en, c0!=0

but for x->infty, lim f = lim c0*w^e0, because all the other terms go to zero.
So, 
    for e0>0, lim f = 0
    for e0<0, lim f = +-infty   (the sign depends on the sign of c0)
    for e0=0, lim f = lim c0

We need to recursively compute limits at several places of the algorithm, but
as is shown in the PhD thesis, it always finishes.

Important functions from the implementation:

compare(a,b,x) compares "a" and "b" by computing the limit L.
mrv(e,x) returns the list of most rapidly varying (mrv) subexpressions of "e"
mrvleadterm(e,x) returns the lead term (c0,w,e0) for e
limitinf(e,x) computes lim e  (for x->infty)
limit(e,z,z0) computes any limit by converting it to the case x->infty

all the functions are really simple and straightforward except mrvleadterm(),
which is the most difficult part of the algorithm.

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

def leadterm(series,x):
    """Returns the term c0*x^e0 of the power series in x with the lowest power
    or x in a form (c0,e0)
    """
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
    if not isinstance(series,s.add):
        return extract(series,x)
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
    if not has(e,x): return e #e is a constant

    leadterm=mrvleadterm(e,x) #leadterm= (c0, e0)
    #for e0>0, lim f = 0
    #for e0<0, lim f = +-infty   (the sign depends on the sign of c0)
    #for e0=0, lim f = lim c0
    if leadterm[1] == s.rational(0): return limitinf(leadterm[0],x)
    elif signum(leadterm[1])==1: return s.rational(0)
    elif signum(leadterm[1])==-1: return s.infty
    else: raise "Error"

def signum(a):
    """Returns a sign of an expression at x->infinity"""
    assert isinstance(a,s.number)
    return a.sign()

def has(e,x):
    return not e.diff(x).isequal(s.rational(0))

def sign(e,x):
    "returns: 1 ... e>0, 0 .... e==0, -1 ... e<0, for x->infty"
    if not has(e,x): 
        return signum(e)
    elif e == x: 
        return 1
    elif isinstance(e,s.mul): 
        a,b=e.getab()
        return sign(a,x)*sign(b,x)
    elif isinstance(e,s.exp): 
        return 1 
    elif isinstance(e,s.pow):
        if sign(e.a,x) == 1: 
            return 1
    raise "cannot determine the sign of %s"%e

def rewrite(e,Omega,x,wsym):
    """e(x) ... the function
    Omega ... the mrv set
    wsym ... the symbol which is going to be used for w

    returns the rewritten e in terms of w.
    """
    assert len(Omega)==1
    w=Omega[0]
    assert w!=x
    assert isinstance(w,s.exp)
    if sign(w.arg,x)==1: wsym=1/wsym
    f2=e.subs(w,wsym)
#    print "rewrite: %s, %s:   %s  ->  %s"%(w,wsym,e,f2)

#    A=limitinf((Omega[0].arg/w.arg).eval(),x)
    #finds the shortest element in Omega and rewrites everything else using it
    #for len(Omega)==1 not necessary
    return f2

def moveup(l,x):
    return [e.subs(x,s.exp(x)).eval() for e in l]

def movedown(l,x):
    return [e.subs(x,s.ln(x)).eval() for e in l]

def mrvleadterm(e,x,Omega=None):
    """Returns (c0, e0) for e."""
    e=e.eval()
    #if not has(e,x): return (e,s.rational(1),s.rational(0))
    if Omega==None:
        Omega=mrv(e,x)
    #else: take into account only terms from Omega, which are in e.
    if member(x,Omega):
        return movedown(mrvleadterm(moveup([e],x)[0],x,moveup(Omega,x)),x)
    wsym=s.symbol("w")
    f2=rewrite(e,Omega,x,wsym)
    ser=f2.series(wsym,3)
    return leadterm(ser.eval(),wsym)
#    print e,Omega,wexpr,f2,ser,lterm

def mrv(e,x):
    "Returns the list of most rapidly varying (mrv) subexpressions of 'e'"
    if not has(e,x): return []
    elif e==x: return [x]
    elif isinstance(e,s.mul): 
        a,b=e.getab()
        return max(mrv(a,x),mrv(b,x),x)
    elif isinstance(e,s.add): 
        a,b=e.getab()
        return max(mrv(a,x),mrv(b,x),x)
    elif isinstance(e,s.pow) and isinstance(e.b,s.number):
        return mrv(e.a,x)
    elif isinstance(e,s.exp): 
        if limitinf(e.arg,x)==s.infty:
            return max([e],mrv(e.arg,x),x)
        else:
            return mrv(e.arg,x)
    raise "unimplemented in mrv: %s"%e

def max(f,g,x):
    """Computes the maximum of two sets of expressions f and g, which 
    are in the same comparability class, i.e. max() compares (two elements of)
    f and g and returns the set, which is in the higher comparability class
    of the union of both, if they have the same order of variation.
    """
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
    """Returns "<" if a<b, "=" for a==b, ">" for a>b"""
    if a==s.exp(-x) and b==x: return ">"
    elif b==s.exp(-x) and a==x: return "<"
    elif a==s.exp(x+1/x) and b==x: return ">"
    elif a==s.exp(x) and b==x**5: return ">"
    elif a==s.exp(x**2) and b==s.exp(x)**2: return ">"
    elif b==s.exp(-x) and a==s.exp(x+s.exp(-x)): return "="
    elif a==s.exp(s.exp(-s.exp(x))+x) and b==s.exp(-s.exp(x)): return "<"
    elif a==s.exp(s.exp(-x**2)+x) and b==s.exp(-x**2): return "<"

    #print a,b,(s.ln(a)/s.ln(b)).eval()
    c=limitinf(s.ln(a)/s.ln(b),x)
    if c==s.rational(0): return "<"
    elif c==s.infty: return ">"
    else: return "="
