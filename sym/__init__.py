"""Computer Algebra System"""
from basic import basic
from symbol import symbol
from functions import sin,cos,tan,exp,ln
from numbers import rational,real,number
from power import pow,pole_error
from add import add,mul,ncmul

"""2006 Ondrej Certik

Basic usage:
>>> import sym
>>> x=sym.symbol("x")
>>> e=1/sym.cos(x)
>>> print e.series(x,10)
1+1/2*x^2+5/24*x^4+61/720*x^6+277/8064*x^8+50521/3628800*x^10


Sympy
=====

Computer algebra system (CAS) in python.

Sympy is a symbolic manipulation package, written in python.
Currently, Sympy capabilities include:
    - basic aritmetics *,/,+,- 
    - basic symplification (like a*b*b+2*b*a*b  -> 3*a*b^2)
    - expansion (like (a+b)^2 -> a^2 + 2*a*b+ b^2)
    - functions (exp, ln, sin, cos, tan, ...)
    - differentiation
    - taylor series
    - basic substitution (like x-> ln(x))
    - arbitrary precision integers and rationals
    - standard (python) floats

Working on:
    - limits, according to: http://www.cybertester.com/data/gruntz.pdf

More distant plans:
    - reach capabilities of Ginac. 

Why another CAS:
----------------

Why not to use CAS which has its own language (like Maple/Mathematica/Maxima)
is best described by: http://www.ginac.de/FAQ.html#whynotmaple

We want to use CAS from a "normal" language, like C, C++, python (maybe
also Java, Ruby, C#). There are currently only 2 libraries, Ginac and Giac
(both for C++), which satisfy this need. Even them are unfortunately too
complicated and difficult to extend. Sympy tries to be as simple as possible,
the primary aim is: easy to learn, easy to extend. Yet Simpy wants to be a full
feature CAS.

Implementation details:
-----------------------

All symbolic things are implemented using subclasses of the "basic" class.
First, you need to create symbols using symbol("x") or numbers using
rational(5) or real(34.3). Then you construct the expression using any class
from the CAS system.  For example add(symbol("a"),symbol("b")) gives an
instance of the add class.  You can call all methods, which the particular
class supports. 

For easier use, there is a syntactic sugar for expressions like:
sym.cos(x)+1 ... basic.__add__(1)  ... add(sym.cos(x),rational(1))
1/sym.cos(x) ... basic.__rdiv__(1) ... 
           ... mul(rational(1),pow(sym.cos(x),rational(-1)))
So, you can write normal expressions using python arithmetics, but from the CAS
point of view, we just need the classes add,mul,pow.

During the construction of the expression, the result is unevaluated, so this
phase is really fast. For computation, the expression needs to be in a
canonical form, this is achieved using the method eval(), which  only performs
unexpensive operations necessary to put the expression in the canonical form.
So the canonical form doesn't mean the simplest possible expresion. The exact
list of operations performed by eval() depends on the implementation.
Obviously, the definition of the canonical form is arbitrary, the only
requirement is that all equivalent expressions must have the same canonical
form. We tried to achieve a canonical, standard form as fast as possible and
also in a way so that the result is what you would write by hand - so "ba + -4
+ b + ab + 4 + (a+b)^2" becomes "2ab + b + (a+b)^2".  The order of terms in the
sum is sorted according to their hash values, so they don't have to be in the
alphabetical order (depends on the hash implementation). 

There is no given requiremens on classes in the library. For example, if they
don't implement the diff() method and you construct an expression using such a
class, then trying to use basic.series() methods will raise an exception of not
founding the diff() method of your class. This "duck typing" has an advantage
that you just implement the functionality which you need. You can define the
function cos as this
class cos(basic):
    def __init__(self,arg):
        basic.__init__(self)
        self.arg=arg
and use it like "1+cos(x)", but if you don't implement the diff() method, you
will not be able to call (1+cos(x)).series().

useful things to implement in new classes: diff, subs, series

the symbolic object is characterized (defined) by the things which it can do.
so implementing more methods like diff, subs etc., you are creating a "shape" of
the symbolic object.

the issue of comparisons - you basically compare evaluated forms. 
so, if expansion operation is not performed upon evaluation (which is
reasonable), then (a+b)^2 != (a^2+2ab+b^2). on the other hand
((a+b)^2).expand == (a^2+2ab+b^2). so, we should define, what we mean by eval,
and take this into account during comparisons.


All objects in the sympy should be immutable - in the sense, that any
operation (like eval() for example) should just return a new instance (it can
return the same instance only if we didn't change). This is a common mistake
to change the current instance, like "self.arg=self.arg.eval()" (wrong!). Use 
"arg=self.arg.eval();return exp(arg)" instead. The object in immutable in the
sense of the symbolic expression they represent. They can modify itself to keep
track of for example the evaluated variable. Or hash. Or they can precalculate
anything regarding the expression they contain. But the expression cannot be
changed. So you can pass any instance to other objects, because you don't have
to worry that it will change, or that this would break anything.

a=symbol("x") and another b=symbol("x") is the same thing, i.e a==b is True.
Note that this is different from Ginac. We chose a==b, because it seems more
natural. In the future, we should implement a=symbol(), to get a unique symbol,
to use for internal computations. 

"""
