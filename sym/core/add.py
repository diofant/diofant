import hashing
from basic import Basic
from numbers import Number,Rational
from power import Pow,pole_error

class Pair(Basic):
    
    def __init__(self,*args):
        Basic.__init__(self)
        if len(args) == 2:
            self.args = [args[0],args[1]]
        elif len(args) == 1:
            self.args = args[0]
            assert len(self.args) > 1
        else:
            raise Exception("accept only 1 or 2 arguments")
            
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        for i in self.args:
            self.mhash.add(i.hash())
        return self.mhash.value
        
    def tryexpand(self,a):
        if isinstance(a,Mul) or isinstance(a,Pow):
            return a.expand()
        else:
            return a
            
    def evalargs(self,a):
        b=[]
        for t in a:
            b.append(t.eval())
        return b
        
    def flatten(self,a):
        """flatten([add(x,4),Mul(a,5),add(x,b),x]) ->
                [x,4,Mul(a,5),x,b,x] if self is add
                [add(x,4),a,5,add(x,b),x] if self is Mul

        returns a copy of "a", where the the classes of the same type as this
        (e.g. if we are Mul, then all the Muls, if we are add, then all the
        adds) are substituted for their arguments. all other classes are
        left intact.
        """
        b=[]
        for x in a:
            if isinstance(x,type(self)):
                b.extend(x.args)
            else:
                b.append(x)
        return b
        
    def coerce(self,a,action):
        """coerce([x,y,z],action) -> action(action(action([],x),y),z)"""
        #equivalent code:
        #exp=[]
        #for x in a:
        #    exp=action(exp,x)
        #return exp
        return reduce(action,a,[])
        
    def coerce_numbers(self,a,action,default):
        """coercenumbers([x,4,a,10],action,Rational(1)) ->
                (action(action(Rational(1),4),10),[x,a])

        picks out the numbers of the list "a" and applies the action on them
        (add or Mul).
        """
        n=default
        b=[]
        for x in a:
            if isinstance(x,Number):
                n=action(n,x)
            else:
                b.append(x)
        return (n,b)
        
    def print_tree(self):
        def indent(s,type=1):
            x = s.split("\n")
            r = "+--%s\n"%x[0]
            for a in x[1:]:
                if a=="": continue
                if type==1:
                    r += "|  %s\n"%a
                else:
                    r += "   %s\n"%a
            return r
        if isinstance(self,Mul):
            f="Mul\n"
        else:
            assert(isinstance(self,Add))
            f="Add\n"
        for a in self.args[:-1]:
            f += indent(a.printtree())
        f += indent(self.args[-1].printtree(),2)
        return f

class Mul(Pair):
    
    def print_normal(self):
        f = ""
        a = self.args
        if isinstance(a[0],Rational):
            if a[0].isminusone():
                f = "-"
                a = self.args[1:]
            elif a[0].isone():
                f = ""
                a = self.args[1:]
        for x in a:
            if isinstance(x,Pair):
                f += "(%s)*"
            else:
                f += "%s*"
        f = f[:-1]
        return f % tuple([str(x) for x in a])
        
    def print_prog(self):
        f = "Mul(%s"+",%s"*(len(self.args)-1)+")"
        return f % tuple([str(x) for x in self.args])
        
    def __str__(self):
        return self.print_normal()
        
    def extractnumericandnonnumeric(self):
        "extract numeric and non numeric part"
        if isinstance(self.args[0],Number):
            return self.getab()
        else:
            return (Rational(1),self)
            
    def get_baseandexp(self,a):
        if isinstance(a,Pow):
            return a.get_baseandexp()
        else:
            return (a,Rational(1))
            
    def eval(self):
        "Flatten, put all Rationals in the front, sort arguments"
        def _mul(exp,x):
            a,aexp = self.get_baseandexp(x)
            e = []
            ok = False
            for y in exp:
                b,bexp = self.get_baseandexp(y)
                if (not ok) and a.isequal(b):
                    e.append(Pow(a,Add(aexp,bexp)).eval())
                    ok = True
                else:
                    e.append(y)
            if not ok: e.append(x)
            return e
        
        if self.evaluated: return self
        a = self.evalargs(self.args)
        a = self.flatten(a)
        #create Powers: a*b*a -> a^2*b
        a = self.coerce(a,_mul)
        #coerce and multiply through the numbers
        n,a = self.coerce_numbers(a,Number.mulnumber,Rational(1))
        if n.iszero(): return Rational(0)
        a.sort(Basic.cmphash)
        #put the number in front of all the other args
        if not n.isone(): a=[n]+a
        if len(a)>1:
            return Mul(a).hold()
        elif len(a)==1:
            return a[0].hold()
        else:
            return Rational(1)
            
    def evalf(self):
        a,b=self.getab()
        return a.evalf()*b.evalf()
        
    def getab(self):
        """Pretend that self=a*b and return a,b
        
        in general, self=a*b*c*d*..., but in many algorithms, we 
        want to have just 2 arguments to Mul. Use this function to 
        simulate this interface. (the returned b = b*c*d.... )
        """
        a=self.args[0]
        if len(self.args)==2:
            b=self.args[1]
        else:
            assert len(self.args) > 2
            b=Mul(self.args[1:])
        return (a,b)
        
    def diff(self,sym):
        r=Rational(0)
        for i in range(len(self.args)):
            d=self.args[i].diff(sym)
            for j in range(len(self.args)):
                if i!=j:
                    d*=self.args[j]
            r+=d
        return r.eval()
        
    def series(self,sym,n):
        """expansion for Mul
        need to handle this correctly:
        (e^x-1)/x  -> 1+x/2+x^2/6
        the first term is (e^x-1)/x evaluated at x=0. normally, we would use
        a limit. but in cas, the limit is computed using series. so we must
        calculate series differently - using the bottom up approach:
        first expand x, then e^x, then e^x-1, and finally (e^x-1)/x
        """
        a,b=self.getab()
        return (a.series(sym,n)*b.series(sym,n)).expand()
        
    def expand(self):
        a,b = self.getab()
        a = self.tryexpand(a)
        b = self.tryexpand(b)
        if isinstance(a,Add):
            d = Rational(0)
            for t in a.args:
                d += (t*b).expand()
            return d.eval()
        elif isinstance(b,Add):
            d=Rational(0)
            for t in b.args:
                d += (a*t).expand()
            return d.eval()
        else:
            return a*b
    def subs(self,old,new):
        a,b=self.getab()
        e=a.subs(old,new)*b.subs(old,new)
        return e.eval()

class Add(Pair):
    
    def print_prog(self):
        f = "Add(%s"+",%s"*(len(self.args)-1)+")"
        return f % tuple([str(x) for x in self.args])
    
    def print_normal(self):
        """Returns a string representation of the expression in self."""
        #old primitive way of printing:
        #f="%s"+"+%s"*(len(self.args)-1)
        #return f%tuple([str(x) for x in self.args])

        #sophisticated way of printing (correctly handles 2+(-x)):
        _n = len(self.args)
        from copy import deepcopy
        _temp_args = deepcopy(self.args)
        f = "%s"
        for i in range(1,_n):
            try:
                if _temp_args[i].args[0] < 0:
                    f += "-%s"
                    _temp_args[i].args[0] = (Rational(-1)*_temp_args[i].args[0]).eval()
                else:
                    f += "+%s"
            except (AttributeError, NotImplementedError):
                if isinstance(_temp_args[1], Number) and _temp_args[i] < 0:
                        # case the second member is a negative Rational
                        f += "-%s"
                        _temp_args[i] = (Rational(-1)*_temp_args[i]).eval()
                else:
                    f += "+%s"
        return f % tuple([str(x) for x in _temp_args])
                
    def __str__(self):
        return self.print_normal()
    
    def getab(self):
        """Pretend that self = a+b and return a,b
        
        in general, self=a+b+c+d+..., but in many algorithms, we 
        want to have just 2 arguments to add. Use this function to 
        simulate this interface. (the returned b = b+c+d.... )
        """
        a=self.args[0]
        if len(self.args)==2:
            b=self.args[1]
        else:
            assert len(self.args) > 2
            b = Add(self.args[1:])
        return (a,b)
    
    def extractnumericandnonnumeric(self,a):
        "extract numeric and non numeric part of 'a'"
        if isinstance(a,Mul):
            return a.extractnumericandnonnumeric()
        elif isinstance(a,Number):
            return (a,Rational(1))
        else:
            return (Rational(1),a)
        
    def eval(self):
        "Flatten, put all Rationals in the back, coerce, sort"
        def _add(exp,x):
            an,a=self.extractnumericandnonnumeric(x)
            e=[]
            ok=False
            for y in exp:
                bn,b=self.extractnumericandnonnumeric(y)
                if (not ok) and a.isequal(b):
                    e.append(Mul(an.addnumber(bn),a).eval())
                    ok=True
                else:
                    e.append(y)
            if not ok: e.append(x)
            return e
        if self.evaluated: return self
        a = self.evalargs(self.args)
        a = self.flatten(a)
        a = self.coerce(a,_add)
        n,a = self.coerce_numbers(a,Number.addnumber,Rational(0))
        a.sort(Basic.cmphash)
        if not n.iszero(): a=[n]+a
        if len(a)>1:
            return Add(a).hold()
        elif len(a)==1:
            return a[0].hold()
        else:
            return Rational(0)
        
    def evalf(self):
        a,b = self.getab()
        return a.evalf() + b.evalf()
    def diff(self,sym):
        d = Rational(0)
        for x in self.args:
            d += x.diff(sym)
        return d.eval()
    
    def expand(self):
        """Tries to expands all the terms in the sum."""
        d = Rational(0)
        for x in self.args:
            d += self.tryexpand(x)
        return d.eval()
    
    def subs(self,old,new):
        d = Rational(0)
        for x in self.args:
            d += x.subs(old,new)
        return d.eval()
    
    def series(self,sym,n):
        """expansion for add
        need to handle this correctly:
        x+1/x
        tries to use Basic.series (which substitutes x->0), if it fails,
        expands term by term
        """
        try:
            return Basic.series(self,sym,n)
        except pole_error:
            a,b = self.getab()
            #there is a cancelation problem here:
            return (a.series(sym,n)+b.series(sym,n)).eval()

class NCMul(Mul):
    
    def print_normal(self):
        f = ""
        a = self.args
        for x in a:
            if isinstance(x,Pair):
                f += "(%s)*"
            else:
                f += "%s*"
        f = f[:-1]
        return f % tuple([str(x) for x in a])
        
    def __str__(self):
        return self.print_normal()
    
    def eval(self):
        "Flatten, put all Rationals in the front, sort arguments"
        def _mul(exp,x):
            a,aexp=self.get_baseandexp(x)
            e=[]
            ok=False
            for y in exp:
                b,bexp=self.get_baseandexp(y)
                if (not ok) and a.isequal(b):
                    e.append(Pow(a,Add(aexp,bexp)).eval())
                    ok=True
                else:
                    e.append(y)
            if not ok: e.append(x)
            return e
        if self.evaluated: return self
        a = self.evalargs(self.args)
        a = self.flatten(a)
        a = self.coerce(a,_mul)
        n,a = self.coerce_numbers(a,Number.mulnumber,Rational(1))
        if n.iszero(): return Rational(0)
        if not n.isone(): a=[n]+a
        if len(a)>1:
            return Mul(a).hold()
        elif len(a)==1:
            return a[0].hold()
        else:
            return Rational(1)
