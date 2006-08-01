import hashing

def c(a):
    if isinstance(a,int):
        return rational(a)
    elif isinstance(a,float):
        return real(a)
    else:
        assert isinstance(a,basic)
        return a
def doadd(a,b):
    return add(c(a),c(b))
def domul(a,b):
    return mul(c(a),c(b))
def dopow(a,b):
    return pow(c(a),c(b))

class pole_error(Exception):
    pass

def sign(x):
    if x < 0: return -1
    elif x==0: return 0
    else: return 1

class basic(object):
    def __init__(self):
        self.evaluated=False;
        self.mhash=0
    def __neg__(self):
        return domul(rational(-1),self)
    def __pos__(self):
        return self
    def __add__(self,a):
        return doadd(self,a)
    def __radd__(self,a):
        return doadd(a,self)
    def __sub__(self,a):
        return doadd(self,-a)
    def __rsub__(self,a):
        return doadd(a,-self)
    def __mul__(self,a):
        return domul(self,a)
    def __rmul__(self,a):
        return domul(a,self)
    def __div__(self,a):
        return domul(self,dopow(a,rational(-1)))
    def __rdiv__(self,a):
        return domul(a,dopow(self,rational(-1)))
    def __pow__(self,a):
        return dopow(self,a)
    def __rpow__(self,a):
        return dopow(a,self)

    def eval(self):
        return self
    def hold(self):
        self.evaluated=True
        return self
    def isequal(self,a):
        return self.hash()==a.hash()
    def cmphash(a,b):
        return sign(a.hash()-b.hash())
    def diffn(self,sym,n):
        while n:
            self=self.diff(sym)
            n-=1
        return self
    def series(self,sym,n):
        f=self
        e=f.subs(sym,rational(0))
        fact=rational(1)
        for i in range(1,n+1):
            fact*=rational(i)
            f=f.diff(sym)
            e+=f.subs(sym,rational(0))*sym**i/fact
        return e.eval()
    def subs(self,old,new):
        if self.isequal(old):
            return new
        else:
            return self

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
        b=[]
        for x in a:
            if isinstance(x,type(self)):
                b.extend(x.args)
            else:
                b.append(x)
        return b
    def coerce(self,a,action):
        "mul the same terms together"
        exp=[]
        for x in a:
            exp=action(exp,x)
        return exp
    def coercenumbers(self,a,action,default):
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
        if isinstance(a[0],number):
            if a[0].isminusone():
                f="-"
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
        a=self.coerce(a,mul2)
        n,a=self.coercenumbers(a,number.mulnumber,rational(1))
        if n.iszero(): return rational(0)
        a.sort(basic.cmphash)
        if not n.isone(): a=[n]+a
        if len(a)>1:
            return mul(a).hold()
        elif len(a)==1:
            return a[0].hold()
        else:
            return rational(1)
    def getab(self):
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
        f="%s"+"+%s"*(len(self.args)-1)
        return f%tuple([str(x) for x in self.args])
    def __str__(self):
        return self.printnormal()
        return self.printprog()
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

    def diff(self,sym):
        d=rational(0)
        for x in self.args:
            d+=x.diff(sym)
        return d.eval()
    def expandterms(self):
        d=rational(0)
        for x in self.args:
            d+=self.tryexpand(x)
        return d.eval()
    def subs(self,old,new):
        d=rational(0)
        for x in self.args:
            d+=x.subs(old,new)
        return d.eval()

class pow(basic):
    def __init__(self,a,b):
        basic.__init__(self)
        self.a=a
        self.b=b
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.add(self.a.hash())
        self.mhash.add(self.b.hash())
        return self.mhash.value
    def printnormal(self):
        f=""
        if isinstance(self.a,pair) or isinstance(self.a,pow):
            f+="(%s)"
        else:
            f+="%s"
        f+="^"
        if isinstance(self.b,pair) or isinstance(self.b,pow):
            f+="(%s)"
        else:
            f+="%s"
        return f%(self.a,self.b)
    def __str__(self):
        return self.printnormal()
    def getbaseandexp(self):
        return (self.a,self.b)
    def eval(self):
        if self.evaluated: return self
        self.a=self.a.eval()
        self.b=self.b.eval()
        if isinstance(self.b,rational) and self.b.iszero():
            return rational(1)
        if isinstance(self.b,rational) and self.b.isone():
            return self.a
        if isinstance(self.a,rational) and self.a.iszero():
            if isinstance(self.b,rational) and self.b.isinteger():
                if self.b.iszero():
                    raise pole_error("pow::eval(): 0^0.")
                elif self.b.getinteger()<0:
                    raise pole_error("pow::eval(): Division by 0.")
            return rational(0)
        if isinstance(self.a,rational) and self.a.isone():
            return rational(1)
        if isinstance(self.a,real) and isinstance(self.b,real):
            return self.a.pownumber(self.b)
        if isinstance(self.a,rational) and isinstance(self.b,rational) and \
                self.b.isinteger():
            return self.a.pownumber(self.b)
        if isinstance(self.a,pow): 
            return pow(self.a.a,self.a.b*self.b).eval()
        return pow(self.a,self.b).hold()
    def diff(self,sym):
        f=self.a
        g=self.b
        return (self*(g*ln(f)).diff(sym)).eval()
    def series(self,sym,n):
        if isinstance(self.a,symbol) and isinstance(self.b,rational):
            if self.b.isinteger():
                return self
        return basic.series(self,sym,n)
    def expand(self):
        if isinstance(self.b,number):
            if self.b.isinteger():
                n=self.b.getinteger()
                if n>1:
                    a=self.a
                    while n>1:
                        a*=self.a
                        n-=1
                    return a.expand()
        return self
    def subs(self,old,new):
        return (self.a.subs(old,new)**self.b.subs(old,new)).eval()

class symbol(basic):
    def __init__(self,name):
        basic.__init__(self)
        self.name=name
    def __str__(self):
        return str(self.name)
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addstr(self.name)
        return self.mhash.value
    def diff(self,sym):
        if self==sym:
            return rational(1)
        else:
            return rational(0)

class number(basic):
    def addnumber(self,a):
        if isinstance(a,real):
            self=self.evalf()
        return self.add(a)
    def mulnumber(self,a):
        if isinstance(a,real):
            self=self.evalf()
        return self.mul(a)
    def pownumber(self,a):
        if isinstance(a,real):
            self=self.evalf()
        return self.pow(a)
    def diff(self,sym):
        return rational(0)

class real(number):
    def __init__(self,num):
        basic.__init__(self)
        if isinstance(num,str):
            num=float(num)
        assert isinstance(num,float)
        self.num=num
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addfloat(self.num)
        return self.mhash.value
    def __str__(self):
        if self.num < 0:
            f="(%d)"
        else:
            f="%d"
        return f%(self.num)
    def add(self,a):
        return real(self.num+a.evalf().num)
    def mul(self,a):
        return real(self.num*a.evalf().num)
    def pow(self,a):
        return real(self.num**a.evalf().num)
    def iszero(self):
        return False
    def isone(self):
        return False
    def isinteger(self):
        return False
    def evalf(self):
        return self

class rational(number):
    def __init__(self,*args):
        basic.__init__(self)
        if len(args)==1:
            p=args[0]
            q=1
        elif len(args)==2:
            p=args[0]
            q=args[1]
        else:
            raise "invalid number of arguments"
        assert (isinstance(p,int) or isinstance(p,long)) and \
                (isinstance(q,int) or isinstance(q,long))
        assert q!=0
        s=sign(p)*sign(q)
        p=abs(p)
        q=abs(q)
        c=self.gcd(p,q)
        self.p=p/c*s
        self.q=q/c
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.p)
        self.mhash.addint(self.q)
        return self.mhash.value
    def gcd(self,a,b):
        while b!=0:
            c=a % b
            a=b
            b=c
        return a
    def __str__(self):
        if self.q == 1:
            if self.p < 0:
                f="(%d)"
            else:
                f="%d"
            return f%(self.p)
        else:
            if self.p < 0:
                f="(%d/%d)"
            else:
                f="%d/%d"
            return f%(self.p,self.q)
    def mul(self,a):
        return rational(self.p*a.p,self.q*a.q)
    def add(self,a):
        return rational(self.p*a.q+self.q*a.p,self.q*a.q)
    def pow(self,a):
        assert a.q==1
        if a.p > 0:
            return rational(self.p**a.p,self.q**a.p)
        else:
            return rational(self.q**(-a.p),self.p**(-a.p))
    def iszero(self):
        return self.p==0 
    def isone(self):
        return self.p==1 and self.q==1
    def isminusone(self):
        return self.p==-1 and self.q==1
    def isinteger(self):
        return self.q==1
    def getinteger(self):
        assert self.isinteger()
        return self.p
    def evalf(self):
        return real(float(self.p)/self.q)
    def diff(self,sym):
        return rational(0)


class function(basic):
    def __init__(self,arg):
        basic.__init__(self)
        self.arg=arg
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.arg.hash())
        return self.mhash.value
    def diff(self,sym):
        return (self.derivative()*self.arg.diff(sym)).eval()
    def subs(self,old,new):
        e=basic.subs(self,old,new)
        if e==self:
            return (type(self)(self.arg.subs(old,new))).eval()
        else:
            return e
    def __str__(self):
        f="%s(%s)"
        return f%(self.getname(),str(self.arg))

class sin(function):
    def getname(self):
        return "sin"
    def derivative(self):
        return cos(self.arg)
    def eval(self):
        if self.evaluated: return self
        if isinstance(self.arg,rational) and self.arg.iszero():
            return rational(0)
        return self.hold()

class cos(function):
    def getname(self):
        return "cos"
    def derivative(self):
        return -sin(self.arg)
    def eval(self):
        if self.evaluated: return self
        if isinstance(self.arg,rational) and self.arg.iszero():
            return rational(1)
        return self.hold()

class tan(function):
    def getname(self):
        return "tan"
    def derivative(self):
        return rational(1)/cos(self.arg)**rational(2)
    def eval(self):
        return (sin(self.arg)/cos(self.arg)).eval()

class ln(function):
    def getname(self):
        return "ln"
    def derivative(self):
        return rational(1)/self.arg
    def eval(self):
        if self.evaluated: return self
        self.arg=self.arg.eval()
        if isinstance(self.arg,rational) and self.arg.isone():
            return rational(0)
        return self.hold()

class exp(function):
    def getname(self):
        return "exp"
    def derivative(self):
        return exp(self.arg)
    def eval(self):
        if self.evaluated: return self
        self.arg=self.arg.eval()
        if isinstance(self.arg,rational) and self.arg.iszero():
            return rational(1)
        return self.hold()

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
