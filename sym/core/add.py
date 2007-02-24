import hashing
from basic import basic
from numbers import number,rational
from power import pow,pole_error

class pair(basic):
    def __init__(self,*args):
        basic.__init__(self)
        if len(args)==2:
            self.args=[args[0],args[1]]
        elif len(args)==1:
            self.args=args[0]
            assert len(self.args)>1
        else:
            raise "accept only 1 or 2 arguments"
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        for i in self.args:
            self.mhash.add(i.hash())
        return self.mhash.value
    def tryexpand(self,a):
        if isinstance(a,mul) or isinstance(a,pow):
            return a.expand()
        else:
            return a
    def evalargs(self,a):
        b=[]
        for t in a:
            b.append(t.eval())
        return b
    def flatten(self,a):
        """flatten([add(x,4),mul(a,5),add(x,b),x]) ->
                [x,4,mul(a,5),x,b,x] if self is add
                [add(x,4),a,5,add(x,b),x] if self is mul

        returns a copy of "a", where the the classes of the same type as this
        (e.g. if we are mul, then all the muls, if we are add, then all the
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
    def coercenumbers(self,a,action,default):
        """coercenumbers([x,4,a,10],action,rational(1)) ->
                (action(action(rational(1),4),10),[x,a])

        picks out the numbers of the list "a" and applies the action on them
        (add or mul).
        """
        n=default
        b=[]
        for x in a:
            if isinstance(x,number):
                n=action(n,x)
            else:
                b.append(x)
        return (n,b)

class mul(pair):
    def printnormal(self):
        f=""
        a=self.args
        if isinstance(a[0],rational):
            if a[0].isminusone():
                f="-"
                a=self.args[1:]
            elif a[0].isone():
                f = ""
                a=self.args[1:]
        for x in a:
            if isinstance(x,pair):
                f+="(%s)*"
            else:
                f+="%s*"
        f=f[:-1]
        return f%tuple([str(x) for x in a])
    def printprog(self):
        f="mul(%s"+",%s"*(len(self.args)-1)+")"
        return f%tuple([str(x) for x in self.args])
    def __str__(self):
        return self.printnormal()
    def extractnumericandnonnumeric(self):
        "extract numeric and non numeric part"
        if isinstance(self.args[0],number):
            return self.getab()
        else:
            return (rational(1),self)
    def getbaseandexp(self,a):
        if isinstance(a,pow):
            return a.getbaseandexp()
        else:
            return (a,rational(1))
    def eval(self):
        "Flatten, put all rationals in the front, sort arguments"
        def mul2(exp,x):
            a,aexp=self.getbaseandexp(x)
            e=[]
            ok=False
            for y in exp:
                b,bexp=self.getbaseandexp(y)
                if (not ok) and a.isequal(b):
                    e.append(pow(a,add(aexp,bexp)).eval())
                    ok=True
                else:
                    e.append(y)
            if not ok: e.append(x)
            return e
        if self.evaluated: return self
        a=self.evalargs(self.args)
        a=self.flatten(a)
        #create powers: a*b*a -> a^2*b
        a=self.coerce(a,mul2)
        #coerce and multiply through the numbers
        n,a=self.coercenumbers(a,number.mulnumber,rational(1))
        if n.iszero(): return rational(0)
        a.sort(basic.cmphash)
        #put the number in front of all the other args
        if not n.isone(): a=[n]+a
        if len(a)>1:
            return mul(a).hold()
        elif len(a)==1:
            return a[0].hold()
        else:
            return rational(1)
    def evalf(self):
        a,b=self.getab()
        return a.evalf()*b.evalf()
    def getab(self):
        """Pretend that self=a*b and return a,b
        
        in general, self=a*b*c*d*..., but in many algorithms, we 
        want to have just 2 arguments to mul. Use this function to 
        simulate this interface. (the returned b = b*c*d.... )
        """
        a=self.args[0]
        if len(self.args)==2:
            b=self.args[1]
        else:
            assert len(self.args) > 2
            b=mul(self.args[1:])
        return (a,b)
    def diff(self,sym):
        r=rational(0)
        for i in range(len(self.args)):
            d=self.args[i].diff(sym)
            for j in range(len(self.args)):
                if i!=j:
                    d*=self.args[j]
            r+=d
        return r.eval()
    def series(self,sym,n):
        """expansion for mul
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
        a,b=self.getab()
        a=self.tryexpand(a)
        b=self.tryexpand(b)
        if isinstance(a,add):
            d=rational(0)
            for t in a.args:
                d+=(t*b).expand()
            return d.eval()
        elif isinstance(b,add):
            d=rational(0)
            for t in b.args:
                d+=(a*t).expand()
            return d.eval()
        else:
            return a*b
    def subs(self,old,new):
        a,b=self.getab()
        e=a.subs(old,new)*b.subs(old,new)
        return e.eval()

class add(pair):
    def printprog(self):
        f="add(%s"+",%s"*(len(self.args)-1)+")"
        return f%tuple([str(x) for x in self.args])
    def printnormal(self):
        """Returns a string representation of the expression in self."""
        if (len(self.args) == 2 and isinstance(self.args[1], mul) 
                and isinstance(self.args[1].args[0],number) 
                and self.args[1].args[0] < 0):
            # if the second member is -something
            f = "%s-%s"*(len(self.args)-1)
            tmp_args = self.args
            tmp_args[1].args[0] = (rational(-1)*tmp_args[1].args[0]).eval()
        elif (len(self.args) == 2 and isinstance(self.args[1], number) and self.args[1] < 0):
            # if the second member is a negative number
            f = "%s-%s"*(len(self.args)-1)
            tmp_args = self.args
            tmp_args[1] = (rational(-1)*tmp_args[1]).eval()
        else:
            f = "%s+%s"*(len(self.args)-1)
            tmp_args = self.args
        return f%tuple([str(x) for x in tmp_args])

    def __str__(self):
        return self.printnormal()
        return self.printprog()
    def getab(self):
        """Pretend that self=a+b and return a,b
        
        in general, self=a+b+c+d+..., but in many algorithms, we 
        want to have just 2 arguments to add. Use this function to 
        simulate this interface. (the returned b = b+c+d.... )
        """
        a=self.args[0]
        if len(self.args)==2:
            b=self.args[1]
        else:
            assert len(self.args) > 2
            b=add(self.args[1:])
        return (a,b)
    def extractnumericandnonnumeric(self,a):
        "extract numeric and non numeric part of 'a'"
        if isinstance(a,mul):
            return a.extractnumericandnonnumeric()
        elif isinstance(a,number):
            return (a,rational(1))
        else:
            return (rational(1),a)
    def eval(self):
        "Flatten, put all rationals in the back, coerce, sort"
        def add2(exp,x):
            an,a=self.extractnumericandnonnumeric(x)
            e=[]
            ok=False
            for y in exp:
                bn,b=self.extractnumericandnonnumeric(y)
                if (not ok) and a.isequal(b):
                    e.append(mul(an.addnumber(bn),a).eval())
                    ok=True
                else:
                    e.append(y)
            if not ok: e.append(x)
            return e
        if self.evaluated: return self
        a=self.evalargs(self.args)
        a=self.flatten(a)
        a=self.coerce(a,add2)
        n,a=self.coercenumbers(a,number.addnumber,rational(0))
        a.sort(basic.cmphash)
        if not n.iszero(): a=[n]+a
        if len(a)>1:
            return add(a).hold()
        elif len(a)==1:
            return a[0].hold()
        else:
            return rational(0)
    def evalf(self):
        a,b=self.getab()
        return a.evalf()+b.evalf()
    def diff(self,sym):
        d=rational(0)
        for x in self.args:
            d+=x.diff(sym)
        return d.eval()
    def expand(self):
        """Tries to expands all the terms in the sum."""
        d=rational(0)
        for x in self.args:
            d+=self.tryexpand(x)
        return d.eval()
    def subs(self,old,new):
        d=rational(0)
        for x in self.args:
            d+=x.subs(old,new)
        return d.eval()
    def series(self,sym,n):
        """expansion for add
        need to handle this correctly:
        x+1/x
        tries to use basic.series (which substitutes x->0), if it fails,
        expands term by term
        """
        try:
            return basic.series(self,sym,n)
        except pole_error:
            a,b=self.getab()
            #there is a cancelation problem here:
            return (a.series(sym,n)+b.series(sym,n)).eval()

class ncmul(mul):
    def printnormal(self):
        f=""
        a=self.args
        for x in a:
            if isinstance(x,pair):
                f+="(%s)*"
            else:
                f+="%s*"
        f=f[:-1]
        return f%tuple([str(x) for x in a])
    def __str__(self):
        return self.printnormal()
    def eval(self):
        "Flatten, put all rationals in the front, sort arguments"
        def mul2(exp,x):
            a,aexp=self.getbaseandexp(x)
            e=[]
            ok=False
            for y in exp:
                b,bexp=self.getbaseandexp(y)
                if (not ok) and a.isequal(b):
                    e.append(pow(a,add(aexp,bexp)).eval())
                    ok=True
                else:
                    e.append(y)
            if not ok: e.append(x)
            return e
        if self.evaluated: return self
        a=self.evalargs(self.args)
        a=self.flatten(a)
        a=self.coerce(a,mul2)
        n,a=self.coercenumbers(a,number.mulnumber,rational(1))
        if n.iszero(): return rational(0)
        if not n.isone(): a=[n]+a
        if len(a)>1:
            return mul(a).hold()
        elif len(a)==1:
            return a[0].hold()
        else:
            return rational(1)
