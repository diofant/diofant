.. _tutorial-manipulation:

==================================
 Advanced Expression Manipulation
==================================

In this section, we discuss some ways that we can perform advanced
manipulation of expressions.

Understanding Expression Trees
==============================

Before we can do this, we need to understand how expressions are represented
in Diofant.  A mathematical expression is represented as a tree.  Let us take
the expression `x^2 + xy`, i.e., ``x**2 + x*y``.  We can see what this
expression looks like internally by using ``repr``

    >>> from diofant import *
    >>> x, y, z = symbols('x y z')

    >>> expr = x**2 + x*y
    >>> repr(expr)
    "Add(Pow(Symbol('x'), Integer(2)), Mul(Symbol('x'), Symbol('y')))"

The easiest way to tear this apart is to look at a diagram of the expression
tree:

.. This comes from dotprint(x**2 + x*y, labelfunc=repr)

.. graphviz::

    digraph{

    # Graph style
    "bgcolor"="transparent"
    "ordering"="out"
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Add(Pow(Symbol('x'), Integer(2)), Mul(Symbol('x'), Symbol('y')))_()" ["color"="black", "label"="Add", "shape"="ellipse"];
    "Pow(Symbol('x'), Integer(2))_(0,)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Symbol('x')_(0, 0)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Integer(2)_(0, 1)" ["color"="black", "label"="Integer(2)", "shape"="ellipse"];
    "Mul(Symbol('x'), Symbol('y'))_(1,)" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Symbol('x')_(1, 0)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Symbol('y')_(1, 1)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Add(Pow(Symbol('x'), Integer(2)), Mul(Symbol('x'), Symbol('y')))_()" -> "Pow(Symbol('x'), Integer(2))_(0,)";
    "Add(Pow(Symbol('x'), Integer(2)), Mul(Symbol('x'), Symbol('y')))_()" -> "Mul(Symbol('x'), Symbol('y'))_(1,)";
    "Pow(Symbol('x'), Integer(2))_(0,)" -> "Symbol('x')_(0, 0)";
    "Pow(Symbol('x'), Integer(2))_(0,)" -> "Integer(2)_(0, 1)";
    "Mul(Symbol('x'), Symbol('y'))_(1,)" -> "Symbol('x')_(1, 0)";
    "Mul(Symbol('x'), Symbol('y'))_(1,)" -> "Symbol('y')_(1, 1)";
    }

.. note::

   The above diagram was made using `Graphviz <http://www.graphviz.org/>`_ and
   the :py:meth:`dotprint <diofant.printing.dot.dotprint>` function.

First, let's look at the leaves of this tree.  Symbols are instances of the
class Symbol.  While we have been doing

    >>> x = symbols('x')

we could have also done

    >>> x = Symbol('x')

Either way, we get a Symbol with the name "x" [#symbols-fn]_.  For the number
in the expression, 2, we got ``Integer(2)``.  ``Integer`` is the Diofant class
for integers.  It is similar to the Python built-in type ``int``, except that
``Integer`` plays nicely with other Diofant types.

When we write ``x**2``, this creates a ``Pow`` object.  ``Pow`` is short for
"power".

    >>> repr(x**2)
    "Pow(Symbol('x'), Integer(2))"

We could have created the same object by calling ``Pow(x, 2)``

    >>> Pow(x, 2)
    x**2

Note that in the ``repr`` output, we see ``Integer(2)``, the Diofant version of
integers, even though technically, we input ``2``, a Python int.  In general,
whenever you combine a Diofant object with a non-Diofant object via some function
or operation, the non-Diofant object will be converted into a Diofant object.  The
function that does this is ``sympify`` [#sympify-fn]_.

    >>> type(2)
    <... 'int'>
    >>> type(sympify(2))
    <class 'diofant.core.numbers.Integer'>

We have seen that ``x**2`` is represented as ``Pow(x, 2)``.  What about
``x*y``?  As we might expect, this is the multiplication of ``x`` and ``y``.
The Diofant class for multiplication is ``Mul``.

    >>> repr(x*y)
    "Mul(Symbol('x'), Symbol('y'))"

Thus, we could have created the same object by writing ``Mul(x, y)``.

    >>> Mul(x, y)
    x*y

Now we get to our final expression, ``x**2 + x*y``.  This is the addition of
our last two objects, ``Pow(x, 2)``, and ``Mul(x, y)``.  The Diofant class for
addition is ``Add``, so, as you might expect, to create this object, we use
``Add(Pow(x, 2), Mul(x, y))``.

    >>> Add(Pow(x, 2), Mul(x, y))
    x**2 + x*y

There is no subtraction class in Diofant.  ``x - y`` is represented as ``x +
-y``, or, more completely, ``x + -1*y``, i.e., ``Add(x, Mul(-1, y))``.

    >>> expr = x - y
    >>> repr(x - y)
    "Add(Symbol('x'), Mul(Integer(-1), Symbol('y')))"

.. dotprint(x - y, labelfunc=repr)

.. graphviz::

    digraph{

    # Graph style
    "bgcolor"="transparent"
    "ordering"="out"
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Add(Symbol('x'), Mul(Integer(-1), Symbol('y')))_()" ["color"="black", "label"="Add", "shape"="ellipse"];
    "Symbol('x')_(0,)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Mul(Integer(-1), Symbol('y'))_(1,)" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Integer(-1)_(1, 0)" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];
    "Symbol('y')_(1, 1)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Add(Symbol('x'), Mul(Integer(-1), Symbol('y')))_()" -> "Symbol('x')_(0,)";
    "Add(Symbol('x'), Mul(Integer(-1), Symbol('y')))_()" -> "Mul(Integer(-1), Symbol('y'))_(1,)";
    "Mul(Integer(-1), Symbol('y'))_(1,)" -> "Integer(-1)_(1, 0)";
    "Mul(Integer(-1), Symbol('y'))_(1,)" -> "Symbol('y')_(1, 1)";
    }

Similarly to subtraction, there is no class in Diofant for division.  Rather,
division is represented by a power of -1.  Hence, we have ``Pow(y, -1)``.  What
if we had divided something other than 1 by ``y``, like ``x/y``?  Let's see.

    >>> expr = x/y
    >>> repr(expr)
    "Mul(Symbol('x'), Pow(Symbol('y'), Integer(-1)))"

.. dotprint(x/y, labelfunc=repr)

.. graphviz::

    digraph{

    # Graph style
    "bgcolor"="transparent"
    "ordering"="out"
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Mul(Symbol('x'), Pow(Symbol('y'), Integer(-1)))_()" ["color"="black", "label"="Mul", "shape"="ellipse"];
    "Symbol('x')_(0,)" ["color"="black", "label"="Symbol('x')", "shape"="ellipse"];
    "Pow(Symbol('y'), Integer(-1))_(1,)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Symbol('y')_(1, 0)" ["color"="black", "label"="Symbol('y')", "shape"="ellipse"];
    "Integer(-1)_(1, 1)" ["color"="black", "label"="Integer(-1)", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Mul(Symbol('x'), Pow(Symbol('y'), Integer(-1)))_()" -> "Symbol('x')_(0,)";
    "Mul(Symbol('x'), Pow(Symbol('y'), Integer(-1)))_()" -> "Pow(Symbol('y'), Integer(-1))_(1,)";
    "Pow(Symbol('y'), Integer(-1))_(1,)" -> "Symbol('y')_(1, 0)";
    "Pow(Symbol('y'), Integer(-1))_(1,)" -> "Integer(-1)_(1, 1)";
    }

We see that ``x/y`` is represented as ``x*y**-1``, i.e., ``Mul(x, Pow(y,
-1))``.

But what about ``x/2``?   Following the pattern of the previous example, we
might expect to see ``Mul(x, Pow(Integer(2), -1))``.  But instead, we have
``Mul(Rational(1, 2), x)``.  Rational numbers are always combined into a single
term in a multiplication, so that when we divide by 2, it is represented as
multiplying by 1/2.

Finally, one last note.  You may have noticed that the order we entered our
expression and the order that it came out from ``repr`` or in the graph were
different.  You may have also noticed this phenonemon earlier in the
tutorial.  For example

     >>> 1 + x
     x + 1

This because in Diofant, the arguments of the commutative operations ``Add`` and
``Mul`` are stored in an arbitrary (but consistent!) order, which is
independent of the order inputted (if you're worried about noncommutative
multiplication, don't be.  In Diofant, you can create noncommutative Symbols
using ``Symbol('A', commutative=False)``, and the order of multiplication for
noncommutative Symbols is kept the same as the input).  Furthermore, as we
shall see in the next section, the printing order and the order in which
things are stored internally need not be the same either.

.. tip::

   The way an expression is represented internally and the way it is printed
   are often not the same.

In general, an important thing to keep in mind when working with Diofant expression
trees is this:  the internal representation of an expression and the way it is
printed need not be the same.  The same is true for the input form.   If some
expression manipulation algorithm is not working in the way you expected it
to, chances are, the internal representation of the object is different from
what you thought it was.

Recursing through an Expression Tree
====================================

Now that you know how expression trees work in Diofant, let's look at how to dig
our way through an expression tree.  Every object in Diofant has two very
important attributes, ``func``, and ``args``.


func
----

``func`` is the head of the object. For example, ``(x*y).func`` is ``Mul``.
Usually it is the same as the class of the object (though there are exceptions
to this rule).

Two notes about ``func``.  First, the class of an object need not be the same
as the one used to create it.  For example

    >>> expr = Add(x, x)
    >>> expr.func
    <class 'diofant.core.mul.Mul'>

We created ``Add(x, x)``, so we might expect ``expr.func`` to be ``Add``, but
instead we got ``Mul``.  Why is that?  Let's take a closer look at ``expr``.

    >>> expr
    2*x

``Add(x, x)``, i.e., ``x + x``, was automatically converted into ``Mul(2,
x)``, i.e., ``2*x``, which is a ``Mul``.   Diofant classes make heavy use of the
``__new__`` class constructor, which, unlike ``__init__``, allows a different
class to be returned from the constructor.

Second, some classes are special-cased, usually for efficiency reasons
[#singleton-fn]_.

    >>> Integer(2).func
    <class 'diofant.core.numbers.Integer'>
    >>> Integer(0).func
    <class 'diofant.core.numbers.Zero'>
    >>> Integer(-1).func
    <class 'diofant.core.numbers.NegativeOne'>

For the most part, these issues will not bother us.  The special classes
``Zero``, ``One``, ``NegativeOne``, and so on are subclasses of ``Integer``,
so as long as you use ``isinstance``, it will not be an issue.

args
----

``args`` are the top-level arguments of the object.  ``(x*y).args`` would be
``(x, y)``.  Let's look at some examples

    >>> expr = 3*y**2*x
    >>> expr.func
    <class 'diofant.core.mul.Mul'>
    >>> expr.args
    (3, x, y**2)

From this, we can see that ``expr == Mul(3, y**2, x)``.  In fact, we can see
that we can completely reconstruct ``expr`` from its ``func`` and its
``args``.

    >>> expr.func(*expr.args)
    3*x*y**2
    >>> expr == expr.func(*expr.args)
    True

Note that although we entered ``3*y**2*x``, the ``args`` are ``(3, x, y**2)``.
In a ``Mul``, the Rational coefficient will come first in the ``args``, but
other than that, the order of everything else follows no special pattern.  To
be sure, though, there is an order.

    >>> expr = y**2*3*x
    >>> expr.args
    (3, x, y**2)

Mul's ``args`` are sorted, so that the same ``Mul`` will have the same
``args``.  But the sorting is based on some criteria designed to make the
sorting unique and efficient that has no mathematical significance.

The ``repr`` form of our ``expr`` is ``Mul(3, x, Pow(y, 2))``.  What if we
want to get at the ``args`` of ``Pow(y, 2)``.  Notice that the ``y**2`` is in
the third slot of ``expr.args``, i.e., ``expr.args[2]``.

    >>> expr.args[2]
    y**2

So to get the ``args`` of this, we call ``expr.args[2].args``.

    >>> expr.args[2].args
    (y, 2)

Now what if we try to go deeper.  What are the args of ``y``.  Or ``2``.
Let's see.

    >>> y.args
    ()
    >>> Integer(2).args
    ()

They both have empty ``args``.  In Diofant, empty ``args`` signal that we have
hit a leaf of the expression tree.

So there are two possibilities for a Diofant expression. Either it has empty
``args``, in which case it is a leaf node in any expression tree, or it has
``args``, in which case, it is a branch node of any expression tree.  When it
has ``args``, it can be completely rebuilt from its ``func`` and its ``args``.
This is expressed in the key invariant.

.. topic:: Key Invariant

   Every well-formed Diofant expression must either have empty ``args`` or
   satisfy ``expr == expr.func(*expr.args)``.

(Recall that in Python if ``a`` is a tuple, then ``f(*a)`` means to call ``f``
with arguments from the elements of ``a``, e.g., ``f(*(1, 2, 3))`` is the same
as ``f(1, 2, 3)``.)

This key invariant allows us to write simple algorithms that walk expression
trees, change them, and rebuild them into new expressions.

Walking the Tree
----------------

With this knowledge, let's look at how we can recurse through an expression
tree.  The nested nature of ``args`` is a perfect fit for recursive functions.
The base case will be empty ``args``.  Let's write a simple function that goes
through an expression and prints all the ``args`` at each level.

    >>> def pre(expr):
    ...     print(expr)
    ...     for arg in expr.args:
    ...         pre(arg)

See how nice it is that ``()`` signals leaves in the expression tree.  We
don't even have to write a base case for our recursion; it is handled
automatically by the for loop.

Let's test our function.

    >>> expr = x*y + 1
    >>> pre(expr)
    x*y + 1
    1
    x*y
    x
    y

Can you guess why we called our function ``pre``?  We just wrote a pre-order
traversal function for our expression tree.   See if you can write a
post-order traversal function.

Such traversals are so common in Diofant that the generator functions
``preorder_traversal`` and ``postorder_traversal`` are provided to make such
traversals easy.  We could have also written our algorithm as

    >>> for arg in preorder_traversal(expr):
    ...     print(arg)
    x*y + 1
    1
    x*y
    x
    y

.. rubric:: Footnotes

.. [#symbols-fn] We have been using ``symbols`` instead of ``Symbol`` because it
  automatically splits apart strings into multiple ``Symbol``\ s.
  ``symbols('x y z')`` returns a tuple of three ``Symbol``\ s.  ``Symbol('x y
  z')`` returns a single ``Symbol`` called ``x y z``.
.. [#sympify-fn] Technically, it is an internal function called ``_sympify``,
  which differs from ``sympify`` in that it does not convert strings.  ``x +
  '2'`` is not allowed.
.. [#singleton-fn] Classes like ``One`` and ``Zero`` are singletonized, meaning
  that only one object is ever created, no matter how many times the class is
  called.  This is done for space efficiency, as these classes are very
  common.  For example, ``Zero`` might occur very often in a sparse matrix
  represented densely.  As we have seen, ``NegativeOne`` occurs any time we
  have ``-x`` or ``1/x``.  It is also done for speed efficiency because
  singletonized objects can be compared by ``is``.  The unique objects for
  each singletonized class can be accessed from the ``S`` object.
