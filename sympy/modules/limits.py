"""
Limits
======

Implemented according to the PhD thesis
http://www.cybertester.com/data/gruntz.pdf, which contains very thorough
descriptions of the algorithm including many examples.  We summarize here the
gist of it.


All functions are sorted according to how rapidly varying they are at infinity
using the following rules. Any two functions f and g can be compared using the
properties of L:

L=lim  log|f(x)| / log|g(x)|           (for x -> infty) 

We define >, < ~ according to::
    
    1. f > g .... L=+-infty 
    
        - f is greater than any power of g
        - f is more rapidly varying than g
        - f goes to infinity/zero faster than g
    
    
    2. f < g .... L=0 
    
        - f is lower than any power of g
    
    3. f ~ g .... L!=0,+-infty 
    
        - both f and g are bounded from above and below by suitable integral powers
        of the other


Examples
========
::
    1 < x < exp(x) < exp(x^2) < exp(exp(x))
    1 ~ 3 ~ -5
    x ~ x^2 ~ x^3 ~ 1/x ~ x^m ~ -x
    exp(x) ~ exp(-x) ~ exp(2x) ~ exp(x)^2 ~ exp(x+exp(-x))
    f ~ 1/f

So we can divide all the functions into comparability classes (x and x^2 is the
same class, exp(x) and exp(-x) is some other class). In principle, we could
compare any two functions, but in our algorithm, we don't compare anything
below f=1 (for example log(x) is below 1), so we set f=1 as the lowest
comparability class. 

Given the function f, we find the list of most rapidly varying (mrv set)
subexpressions of it. This list belongs to the same comparability class. Let's
say it is {exp(x), exp(2x)}. Using the rule f ~ 1/f we find an element "w"
(either from the list or a new one) from the same comparability class which
goes to zero at infinity. In our example we set w=exp(-x) (but we could also
set w=exp(-2x) or w=exp(-3x) ...). We rewrite the mrv set using w, in our case
{1/w,1/w^2}, and substitute it into f. Then we expand f into a series in w::

    f = c0*w^e0 + c1*w^e1 + ... + O(w^en),        where e0<e1<...<en, c0!=0

but for x->infty, lim f = lim c0*w^e0, because all the other terms go to zero,
because w goes to zero faster than the ci and ei. So::

    for e0>0, lim f = 0
    for e0<0, lim f = +-infty   (the sign depends on the sign of c0)
    for e0=0, lim f = lim c0

We need to recursively compute limits at several places of the algorithm, but
as is shown in the PhD thesis, it always finishes.

Important functions from the implementation:

compare(a,b,x) compares "a" and "b" by computing the limit L.
mrv(e,x) returns the list of most rapidly varying (mrv) subexpressions of "e"
rewrite(e,Omega,x,wsym) rewrites "e" in terms of w
leadterm(f,x) returns the lowest power term in the series of f
mrvleadterm(e,x) returns the lead term (c0,e0) for e
limitinf(e,x) computes lim e  (for x->infty)
limit(e,z,z0) computes any limit by converting it to the case x->infty

all the functions are really simple and straightforward except rewrite(),
which is the most difficult part of the algorithm.

"""

import sympy as s
from sympy import Basic, mhash
from sympy.core.prettyprint import StringPict

from decorator import decorator

#Debugging:
#import the limits.py in your code and set limits.debug=True. 
#this will print a nice tree of recursive calls to all methods here, which
#are decorated with @decorator(maketree)
#you can apply this decorator to any method here and it will be included
#in the tree.
debug = False

def tree(subtrees):
    def indent(s,type=1):
        x = s.split("\n")
        r = "+-%s\n"%x[0]
        for a in x[1:]:
            if a=="": continue
            if type==1:
                r += "| %s\n"%a
            else:
                r += "  %s\n"%a
        return r
    if len(subtrees)==0: return ""
    f="";
    for a in subtrees[:-1]:
        f += indent(a)
    f += indent(subtrees[-1],2)
    return f

tmp=[]
iter=0
def maketree(f,*args,**kw):
    global tmp
    global iter
    if debug:
        oldtmp=tmp
        tmp=[]
        iter+=1

    r = f(*args,**kw)

    if debug:
        iter-=1
        s = "%s%s = %s\n" % (f.func_name,args,r)
        if tmp!=[]: s += tree(tmp)
        tmp=oldtmp
        tmp.append(s)
        if iter == 0: 
            print tmp[0]
            tmp=[]
    return r

def getattr_(obj, name, default_thunk):
    "Similar to .setdefault in dictionaries."
    try:
        return getattr(obj, name)
    except AttributeError:
        default = default_thunk()
        setattr(obj, name, default)
        return default

@decorator
def memoize(func, *args):
    dic = getattr_(func, "memoize_dic", dict)
    # memoize_dic is created at the first call
    argshash=tuple([x.hash() for x in args])
    if argshash in dic:
        return dic[argshash]
    else:
        result = func(*args)
        dic[argshash] = result
        return result

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

def limit(e,z,z0):
    """Currently only limit z->z0+"""
    x=s.Symbol("x",True)
    e0=e.subs(z,z0+1/x)
    return limitinf(e0,x)

#@decorator(maketree)
@memoize
def limitinf(e,x):
    """Limit e(x) for x-> infty"""
    if not e.has(x): return e #e is a constant
    c0,e0 = mrv_leadterm(e,x) 
    sig=sign(e0,x)
    if sig==1: return s.Rational(0) # e0>0: lim f = 0
    elif sig==-1: 
        s.infty.sig=sign(c0,x)
        return s.infty #e0<0: lim f = +-infty   (the sign depends on the sign of c0)
    elif sig==0: return limitinf(c0,x) #e0=0: lim f = lim c0

@memoize
def sign(e,x):
    """Returns a sign of an expression at x->infty.
    
        e>0 ... 1
        e==0 .. 0
        e<0 ... -1
    """
    #print "sign:",e
    if isinstance(e,s.Number):
        return e.sign()
    elif not e.has(x):
        return e.evalf() > 0
    elif e == x: 
        return 1
    elif isinstance(e,s.Mul): 
        a,b=e.getab()
        return sign(a,x)*sign(b,x)
#    elif isinstance(e,s.add): 
#        a,b=e.getab()
#        return sign(a,x)*sign(b,x)
    elif isinstance(e,s.exp): 
        return 1 
    elif isinstance(e,s.Pow):
        if sign(e.base,x) == 1: 
            return 1
    elif isinstance(e,s.log): 
        return sign(e.arg-1,x)
    elif isinstance(e,s.Add):
        #print limitinf(e,x) 
        #print sign(limitinf(e,x),x) 
        return sign(limitinf(e,x),x) #FIXME this is wrong for -infty
    raise "cannot determine the sign of %s"%e

def tryexpand(a):
    if isinstance(a,s.Mul) or isinstance(a,s.Pow) or isinstance(a,s.Add):
        return a.expand()
    else:
        return a

#@decorator(maketree)
def rewrite(e,Omega,x,wsym):
    """e(x) ... the function
    Omega ... the mrv set
    wsym ... the symbol which is going to be used for w

    returns the rewritten e in terms of w. and log(w)
    """
    for t in Omega: assert isinstance(t,s.exp)
    assert len(Omega)!=0
    def cmpfunc(a,b):
        #FIXME: this is really, really slow...
        return -cmp(len(mrv(a,x)), len(mrv(b,x)))
    #sort Omega (mrv set) from the most complicated to the simplest ones
    #the complexity of "a" from Omega: the length of the mrv set of "a"
    Omega.sort(cmp=cmpfunc)
    g=Omega[-1] #g is going to be the "w" - the simplest one in the mrv set
    assert isinstance(g,s.exp) #all items in Omega should be exponencials
    sig= (sign(g.arg,x)==1) 
    if sig: wsym=1/wsym #if g goes to infty, substitute 1/w
    #O2 is a list, which results by rewriting each item in Omega using "w"
    O2=[]
    for f in Omega: 
        assert isinstance(f,s.exp) #all items in Omega should be exponencials
        c=mrv_leadterm(f.arg/g.arg,x)
        assert c[1]==0
        O2.append(s.exp(tryexpand(f.arg-c[0]*g.arg))*wsym**c[0])
    #Remember that Omega contains subexpressions of "e". So now we find
    #them in "e" and substitute them for our rewriting, stored in O2
    f=e 
    for a,b in zip(Omega,O2):
        f=f.subs(a,b)

    #tmp.append("Omega=%s; O2=%s; w=%s; wsym=%s\n"%(Omega,O2,g,wsym))

    #finally compute the logarithm of w (logw). 
    logw=g.arg
    if sig: logw=-logw     #log(w)->log(1/w)=-log(w)
    return f,logw

def moveup(l,x):
    return [e.subs(x,s.exp(x)) for e in l]

def movedown(l,x):
    return [e.subs(x,s.log(x)) for e in l]

def subexp(e,sub):
    """Is "sub" a subexpression of "e"? """
    n=s.Symbol("x",True)
    #we substitute some symbol for the "sub", and if the 
    #expression changes, the substitution was successful, thus the answer
    #is yes.
    return e.subs(sub,n)!=e

#@decorator(maketree)
def mrv_leadterm(e,x,Omega=[]):
    """Returns (c0, e0) for e."""
    if not e.has(x): return (e,s.Rational(0))
    Omega = [t for t in Omega if subexp(e,t)]
    if Omega == []:
        Omega = mrv(e,x)
    if member(x,Omega):
        return movedown(mrv_leadterm(moveup([e],x)[0],x,moveup(Omega,x)),x)
    wsym = s.Symbol("w",True)
    f,logw=rewrite(e,Omega,x,wsym)
    series=f.series(wsym,1)
    n = 2
    while series==0 and n<10:
        series=f.series(wsym,n)
        n += 1
    assert series!=0
    series=series.subs(s.log(wsym),logw)
    return series.leadterm(wsym)

#@decorator(maketree)
@memoize
def mrv(e,x):
    "Returns the list of most rapidly varying (mrv) subexpressions of 'e'"
    if not e.has(x): return []
    elif e == x: return [x]
    elif isinstance(e, s.Mul): 
        a,b = e.getab()
        return max(mrv(a,x),mrv(b,x),x)
    elif isinstance(e, s.Add): 
        a,b = e.getab()
        return max(mrv(a,x),mrv(b,x),x)
    elif isinstance(e, s.Pow):
        if e.exp.has(x):
            return mrv(s.exp(e.exp * s.log(e.base)),x)
        else:
            return mrv(e.base,x)
    elif isinstance(e, s.log): 
        return mrv(e.arg, x)
    elif isinstance(e, s.exp): 
        if limitinf(e.arg,x) == s.infty:
            return max([e],mrv(e.arg, x), x)
        else:
            return mrv(e.arg,x)
    elif isinstance(e, s.Function): 
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
    c=limitinf(s.log(a)/s.log(b),x)
    if c==s.Rational(0): return "<"
    elif c==s.infty: return ">"
    else: return "="

class Limit(Basic):

    def __init__(self,e,x,x0):
        Basic.__init__(self)
        self.e=self.sympify(e)
        self.x=self.sympify(x)
        self.x0=self.sympify(x0)

    def doit(self):
        return limit(self.e,self.x,self.x0)

    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.e.hash())
        self.mhash.addint(self.x.hash())
        self.mhash.addint(self.x0.hash())
        return self.mhash.value

    def print_pretty(self):
         e, x, t = [a.print_pretty() for a in (self.e,self.x,self.x0)]
         return StringPict('lim').below(StringPict.next(x, '->', t)) \
                 .right(' ', e)

def limit(e,z,z0):
    """Currently only limit z->z0+"""
    x=s.Symbol("x",True)
    e0=e.subs(z,z0+1/x)
    return limitinf(e0,x)
