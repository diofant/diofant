#! /usr/bin/python

import unittest
import sym as g

class test_sym(unittest.TestCase):
    def e(self,a,b):
        self.failUnless(a.eval().isequal(b.eval()))
    def ne(self,a,b):
        self.failIf(a.eval().isequal(b.eval()))

    def testsymbol(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.symbol("c")
        self.assertEqual(str(a),"a")
        self.assertEqual(str(b),"b")
        e=a*b
        self.e(e,a*b)
        self.e(a*b*b,a*b**2)
        self.e(a*b*b+c,c+a*b**2)
        self.e(a*b*b-c,-c+a*b**2)
    def testarit(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.symbol("c")
        p=g.rational(5)
        e=a*b
        self.e(e,a*b)
        e=a*b+b*a
        self.e(e,2*a*b)
        e=a*b+b*a+a*b+p*b*a
        self.e(e,8*a*b)
        e=a*b+b*a+a*b+p*b*a+a
        self.e(e,a+8*a*b)
        e=a+a
        self.e(e,2*a)
        e=a+b+a
        self.e(e,b+2*a)
        e=a+b*b+a+b*b
        self.e(e,2*a+2*b**2)
        e=a+g.rational(2)+b*b+a+b*b+p
        self.e(e,7+2*a+2*b**2)
        e=(a+b*b+a+b*b)*p
        self.e(e,5*(2*a+2*b**2))
        e=(a*b*c+c*b*a+b*a*c)*p
        self.e(e,15*a*b*c)
        e=(a*b*c+c*b*a+b*a*c)*p-g.rational(15)*a*b*c
        self.e(e,g.rational(0))
        e=g.rational(50)*(a-a)
        self.e(e,g.rational(0))
        e=b*a-b-a*b+b
        self.e(e,g.rational(0))
        e=a*b+c**p
        self.e(e,a*b+c**5)
        e=a/b
        self.e(e,a*b**(-1))
        e=a*2*2
        self.e(e,4*a)
        e=2+a*2/2
        self.e(e,2+a)
        e=2-a-2
        self.e(e,-a)
        e=2*a*2
        self.e(e,4*a)
        e=2/a/2
        self.e(e,a**(-1))
        e=2**a**2
        self.e(e,2**(a**2))
    def testeval(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.rational(1)
        p=g.rational(5)
        e=(a*b+c+p).eval()
        self.e(e,a*b+6)
        e=(c+a+p).eval()
        self.e(e,a+6)
        e=(c+a-p).eval()
        self.e(e,a+(-4))
        e=(c+a+b*c+a-p).eval()
        self.e(e,2*a+b+(-4))
        e=(a+a).eval()
        self.e(e,2*a)
        e=(a+p+a).eval()
        self.e(e,2*a+5)
        e=(a*g.rational(2)+p+a)
        self.e(e,a*2+5+a)
        e=(a*g.rational(2)+p+a).eval()
        self.e(e,3*a+5)
        e=(c+p).eval()
        self.e(e,g.rational(6))
        e=(a*g.rational(2)+a).eval()
        self.e(e,3*a)
        e=b+a-b
        self.e(e,a)
    def xtestcmp(self):
        b=g.symbol("b")
        a=g.symbol("a")
        c=g.symbol("c")
        e1=a+b
        e2=a*b
        self.assertEqual(e1.isequal(e2),False)

        self.assertEqual(a.cmp(b),-1)
        self.assertEqual(b.cmp(a),1)
        self.assertEqual(b.cmp(b),0)
        self.assertEqual(a.cmp(a),0)
        l=[c,a,b]
        l.sort(g.cmp_expr)
        self.assertEqual(l,[a,b,c])
        self.assertNotEqual(l,[a,c,b])
        self.assertNotEqual(l,[c,a,b])

        n2=g.rational(2)
        n1=g.rational(1)
        n1p=g.rational(1)

        l=[c,n2,n1p,a,b,n1]
        l.sort(g.cmp_expr)
        self.failUnless(l==[n1,n1p,n2,a,b,c] or l==[n1p,n1,n2,a,b,c])

        p1=a**a
        p2=a**b
        p3=b**b
        l=[c,p3,n2,p1,a,b,p2,n1]
        l.sort(g.cmp_expr)
        #self.assertEqual(l,[n1,n2,a,b,c,p1,p2,p3])

        p1=a**2
        p2=a**3
        p3=b**3
        e1=a*b**2
        e2=a**3*b**2
        e3=a**2*b**2
        e4=a**2*b**2*c**2
        e4=e4.eval()
        self.assertEqual(str(e1.gethighestpower()),"2")
        self.assertEqual(str(e2.gethighestpower()),"3")
        l=[e2,c,p3,e4,e3,n2,p1,e1,a,b,p2,n1]
        l.sort(g.cmp_expr)
        #print [str(x) for x in l]
        #self.assertEqual(l,[n1,n2,a,b,c,p1,e1,e3,e4,p2,p3,e2])

        e1=6*a*g.cos(a)
        e2=6*g.sin(a)
        self.assertEqual(e1.cmp(e2),1)
        self.assertEqual(e2.cmp(e1),-1)

        e1=6*g.sin(a)
        e2=6*g.sin(a)
        self.assertEqual(e1.cmp(e2),0)
        self.assertEqual(e2.cmp(e1),0)

        e1=6*g.sin(a)
        e2=6*g.sin(2*a)
        self.assertEqual(e1.cmp(e2),-1)
        self.assertEqual(e2.cmp(e1),1)
    def testdiv(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.symbol("c")
        e=a/b
        self.e(e,a*b**(-1))
        e=a/b+c/2
        self.e(e,a*b**(-1)+g.rational(1)/2*c)
        e=(1-b)/(b-1)
        self.e(e,(1+-b)*((-1)+b)**(-1))
    def testpow(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.symbol("c")
        n1=g.rational(1)
        n2=g.rational(2)
        n5=g.rational(5)
        e=a*a
        self.e(e,a**2)
        e=a*a*a
        self.e(e,a**3)
        e=a*a*a*a**g.rational(6)
        self.e(e,a**9)
        e=a*a*a*a**g.rational(6)-a**g.rational(9)
        self.e(e,g.rational(0))
        e=a**(b+c)*a**(-b)
        self.e(e,a**c)
        e=a**(b+c)*a*a**(-b)*a**(-c)/a
        self.e(e,g.rational(1))
        e=a**(b-b)
        self.e(e,g.rational(1))
        e=(a-a)**b
        self.e(e,g.rational(0))
        e=(a+g.rational(1)-a)**b
        self.e(e,g.rational(1))

        e=(a+b+c)**n2
        self.e(e,(a+b+c)**2)
        self.e(e.expand(),2*b*c+2*a*c+2*a*b+a**2+c**2+b**2)

        e=(a+b)**n2
        self.e(e,(a+b)**2)
        self.e(e.expand(),2*a*b+a**2+b**2)

        e=(a+b)**(n1/n2)
        self.e(e,(a+b)**(g.rational(1)/2))
        self.e(e.expand(),(a+b)**(g.rational(1)/2))

        n=n5**(n1/n2)
        self.e(n,g.rational(5)**(g.rational(1)/2))
        e=n*a*b-n*b*a
        self.e(e,g.rational(0))
        e=n*a*b+n*b*a
        self.e(e,2*a*b*5**(g.rational(1)/2))
        self.e(e.eval().diff(a),2*b*5**(g.rational(1)/2))
        self.e(e.diff(a),2*b*5**(g.rational(1)/2))
        e=a/b**2
        self.e(e,a*b**(-2))
        e=g.ln(1/a).eval()
        self.assertNotEqual(str(e),"ln(1*a^(-1))")
        self.assertEqual(str(e),"ln(a^(-1))")
    def testexpand(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.symbol("c")
        p=g.rational(5)
        e=(a+b)*c
        self.e(e.eval(),c*(a+b))
        self.e((e.eval().expand()-a*c-b*c).eval(),g.rational(0))
        e=(a+b)*(a+b)
        self.e(e.eval(),(a+b)**2)
        self.e(e.eval().expand(),2*a*b+a**2+b**2)
        e=(a+b)*(a+b)**g.rational(2)
        self.e(e.eval(),(a+b)**3)
        self.e(e.expand(),3*b*a**2+3*a*b**2+a**3+b**3)
        self.e(e.eval().expand(),3*b*a**2+3*a*b**2+a**3+b**3)
        e=(a+b)*(a+c)*(b+c)
        self.e(e.eval(),(a+c)*(a+b)*(b+c))
        self.e(e.expand(),2*a*b*c+b*a**2+c*a**2+b*c**2+a*c**2+c*b**2+a*b**2)
        e=(a+g.rational(1))**p
        self.e(e.eval(),(1+a)**5)
        self.e(e.expand(),1+5*a+10*a**2+10*a**3+5*a**4+a**5)
        e=(a+b+c)*(a+c+p)
        self.e(e.eval(),(5+a+c)*(a+b+c))
        self.e(e.expand(),5*a+5*b+5*c+2*a*c+b*c+a*b+a**2+c**2)
        x=g.symbol("x")
        s=g.exp(x*x)-1
        e=s.series(x,4)/x**2
        e=e.eval()
        self.e(e,(x**2+g.rational(1)/2*x**4)*x**(-2))
        self.e(e.expand(), 1+g.rational(1)/2*x**2)
    def testeval(self):
        a=g.symbol("a")
        b=g.symbol("b")
        e=a+b+a+b
        s1=str(e)
        e.eval()
        s2=str(e)
        self.assertEqual(s1,s2)
    def xtestfind(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.symbol("c")
        p=g.rational(5)
        e=a*b+b**p
        self.failUnless(e.find(b))
        self.failIf(e.find(c))
    def testdiff(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.symbol("c")
        p=g.rational(5)
        e=a*b+b**p
        self.e(e.diff(a),b)
        self.e(e.diff(b),a+5*b**4)
        self.e(e.diff(b).diff(a),g.rational(1))
        e=a*(b+c)
        self.e(e.diff(a),b+c)
        self.e(e.diff(b),a)
        self.e(e.diff(b).diff(a),g.rational(1))
        e=c**p
        self.e(e.diffn(c,6),g.rational(0))
        self.e(e.diffn(c,5),g.rational(120))
        e=c**g.rational(2)
        self.e(e.diff(c),2*c)
        e=(a*b*c).eval()
        self.e(e.diff(c),a*b)
    def testfunc(self):
        a=g.symbol("a")
        b=g.symbol("b")
        c=g.symbol("c")
        p=g.rational(5)
        e=a*b+g.sin(b**p)
        self.e(e,a*b+g.sin(b**5))
        self.e(e.diff(a),b)
        self.e(e.diff(b),a+5*b**4*g.cos(b**5))
        e=g.tan(c)
        self.e(e,g.tan(c))
        self.e(e.diff(c),g.cos(c)**(-2))
        e=c*g.ln(c)-c
        self.e(e,-c+c*g.ln(c))
        self.e(e.diff(c),g.ln(c))
        e=g.ln(g.sin(c))
        self.e(e,g.ln(g.sin(c)))
        self.e(e.diff(c),g.sin(c)**(-1)*g.cos(c))
        self.ne(e.diff(c),g.cos(c)**(-1)*g.sin(c))
        self.ne(e.diff(c),g.sin(c)**(-2)*g.cos(c))
        self.ne(e.diff(c),g.sin(c)**(-3)*g.cos(c))
        t=g.rational(2)
        e=(t**a/g.ln(t)).eval()
        self.e(e,2**a*g.ln(g.rational(2))**(-1))
        self.e(e.diff(a),2**a)
    def testhash(self):
        n0=g.rational(-0)
        n1=g.rational(-1)
        n2=g.rational(-2)
        n3=g.rational(-3)
        self.failIf(n1.hash()==n2.hash())
        self.failIf(n1.hash()==n3.hash())
        self.failIf(n0.hash()==n1.hash())
        x=g.symbol("x")
        y=g.symbol("y")
        z=g.symbol("y1")
        z2=g.symbol("y1")
        self.failIf(x.hash()==y.hash())
        self.failIf(x.hash()==z.hash())
        self.failIf(y.hash()==z.hash())
        self.failUnless(z.hash()==z2.hash())
    def testdiff2(self):
        n3=g.rational(3)
        n2=g.rational(2)
        n6=g.rational(6)
        x=g.symbol("x")
        c=g.symbol("c")
        e=n3*(-n2+x**n2)*g.cos(x)+x*(-n6+x**n2)*g.sin(x)
        self.e(e,3*((-2)+x**2)*g.cos(x)+x*((-6)+x**2)*g.sin(x))
        self.e(e.diff(x).expandterms(),x**3*g.cos(x))

        e=(x+1)**3
        self.e(e.diff(x),3*(x+1)**2)
        e=x*(x+1)**3
        self.e(e.diff(x),(x+1)**3+3*x*(x+1)**2)
        e=(2*g.exp(x*x)*x).eval()
        self.e(e.diff(x),2*g.exp(x*x)+4*x**2*g.exp(x*x))

        #e=g.ln(x)/(c**3)-g.ln(-c**3+x**3)/(3*c**3)
        #self.assertEqual(str(e.eval()),"tan(x*2^(-1))")
        #self.assertEqual(str(e.diff(x)),"tan(x*2^(-1))")
    def testsubs(self):
        n3=g.rational(3)
        n2=g.rational(2)
        n6=g.rational(6)
        x=g.symbol("x")
        c=g.symbol("c")
        e=x
        self.assertEqual(str(e),"x")
        e=e.subs(x,n3)
        self.assertEqual(str(e),"3")

        e=2*x
        self.e(e,2*x)
        e=e.subs(x,n3)
        self.assertEqual(str(e),"6")

        e=(g.sin(x)**2).diff(x)
        self.e(e,2*g.sin(x)*g.cos(x))
        e=e.subs(x,n3)
        self.e(e,2*g.cos(n3)*g.sin(n3))

        e=(g.sin(x)**2).diff(x)
        self.e(e,2*g.sin(x)*g.cos(x))
        e=e.subs(g.sin(x),g.cos(x))
        self.e(e,2*g.cos(x)**2)

    def testseries(self):
        n3=g.rational(3)
        n2=g.rational(2)
        n6=g.rational(6)
        x=g.symbol("x")
        c=g.symbol("c")
        e=g.sin(x)
        self.assertEqual(str(e),"sin(x)")
        self.assertEqual(str(e.series(x,0)),"0")
        self.assertEqual(str(e.series(x,1)),"x")
        self.assertEqual(str(e.series(x,2)),"x")
        self.e(e.series(x,3),x+(-g.rational(1)/6)*x**3)
        self.e(e.series(x,4),x+(-g.rational(1)/6)*x**3)

        e=((g.exp(x)-1)/x).eval()
        self.e(e.series(x,1),g.rational(1))
        self.assertRaises(g.pole_error,g.basic.series,e,x,0)

        #e=2*g.sin(x)*g.cos(x)
        #print
        #print e.series(x,5)
        #e=g.sin(2*x)
        #e=g.tan(2*x)
        #e=1/g.cos(x)
        #print e.series(x,8)
    def testrational(self):
        n1=g.rational(1,4)
        n2=g.rational(1,3)
        n3=g.rational(2,4)
        n4=g.rational(2,-4)
        n5=g.rational(0)
        n6=g.rational(1)
        n7=g.rational(3)
        n8=g.rational(-3)
        self.assertEqual(str(n1.mul(n2)),"1/12")
        self.assertEqual(str(n1.mul(n2)),"1/12")
        self.assertEqual(str(n3),"1/2")
        self.assertEqual(str(n1.mul(n3)),"1/8")
        self.assertEqual(str(n1.add(n3)),"3/4")
        self.assertEqual(str(n1.add(n2)),"7/12")
        self.assertEqual(str(n1.add(n4)),"(-1/4)")
        self.assertEqual(str(n4.mul(n4)),"1/4")
        self.assertEqual(str(n4.add(n2)),"(-1/6)")
        self.assertEqual(str(n4.add(n5)),"(-1/2)")
        self.assertEqual(str(n4.mul(n5)),"0")
        self.assertEqual(str(n3.add(n4)),"0")
        self.assertEqual(str(n1.pow(n7)),"1/64")
        self.assertEqual(str(n2.pow(n7)),"1/27")
        self.assertEqual(str(n2.pow(n8)),"27")
        self.assertEqual(str(n7.pow(n8)),"1/27")
    def testncmul(self):
        x=g.symbol("x")
        c=g.symbol("c")
        e=g.mul(x,c)+g.mul(c,x)
        self.assertEqual(str(e),"x*c+c*x")
        self.e(e,2*c*x)
        e=g.ncmul(c,x)+g.ncmul(x,c)
        self.assertEqual(str(e),"c*x+x*c")
        self.e(e,g.ncmul(c,x)+g.ncmul(x,c))
        
class test_arithmetic(unittest.TestCase):

    def dotest(self,s):
        def t(a,b):
            s(a,b)
            s(b,a)
        a=g.rational(2)
        b=g.real("1.3") 
        c=g.symbol("x")
        d=g.symbol("y")
        e=pow(c,d)*d
        f=5
        h=5.5

        t(a,a)
        t(a,b)
        t(a,c)
        t(a,d)
        t(a,e)
        t(a,f)
        t(a,h)

        t(b,b)
        t(b,c)
        t(b,d)
        t(b,e)
        t(b,f)
        t(b,h)

        t(c,c)
        t(c,d)
        t(c,e)
        t(c,f)
        t(c,h)

        t(d,d)
        t(d,e)
        t(d,f)
        t(d,h)

        t(e,e)
        t(e,f)
        t(e,h)

        #t(f,f)
        #t(f,h)

        #t(h,h)

    def testbasic(self):
        def s(a,b):
            x= a
            x= +a
            x= -a
            x= a+b;x.eval()
            x= a-b;x.eval()
            x= a*b;x.eval()
            x= a/b;x.eval()
            x= a**b;x.eval()
        self.dotest(s)

    def testibasic(self):
        def s(a,b):
            x= a
            x+=b;x.eval()
            x= a
            x-=b;x.eval()
            x= a
            x*=b;x.eval()
            x= a
            x/=b;x.eval()
        self.dotest(s)

if __name__ == "__main__":
    #unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
    unittest.main()
