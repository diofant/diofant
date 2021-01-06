.. _tutorial-printing:

==========
 Printing
==========

..
    >>> init_printing(pretty_print=True, use_unicode=True)

As we have already seen, Diofant can pretty print its output using
Unicode characters.  This is a short introduction to the most common
printing options available in Diofant.  The most common ones are

- `Str`_
- `Repr`_
- `2D Pretty Printer`_ (Unicode or ASCII)
- `LaTeX`_
- `Dot`_

In addition to these, there are also "printers" that can output
Diofant objects to code, such as C, Fortran, or Mathematica.

Best printer is enabled automatically for interactive session
(i.e. `\LaTeX` in the IPython notebooks, pretty printer in the IPython
console or str printer in the Python console).  If you want manually
configure pretty printing, please use the
:func:`~diofant.interactive.printing.init_printing` function.

Lets take this simple expression

    >>> expr = Integral(sqrt(1/x))

and try several available printers.

Str
===

To get a string form of an expression, use :class:`str`.  This is also
the form that is produced by :func:`print`.  String forms are designed
to be easy to read and mostly to be in a form that is a correct Python syntax
so that it can be copied and pasted.

    >>> str(expr)
    'Integral(sqrt(1/x), x)'
    >>> print(expr)
    Integral(sqrt(1/x), x)

Repr
====

The repr form of an expression is designed to show the exact form of
an expression, it would yield an object with the same value when
passed to :func:`eval`.  To get it, use :func:`repr`.

    >>> repr(expr)
    "Integral(Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2)), Tuple(Symbol('x')))"

The repr form is mostly useful for understanding how an expression is
built internally.

.. _d-pretty-printer:

2D Pretty Printer
=================

A two-dimensional (2D) textual representation of the expression can be
obtained with :func:`~diofant.printing.pretty.pretty.pretty`.

    >>> pretty(expr)
    '⌠           \n⎮     ___   \n⎮    ╱ 1    \n⎮   ╱  ─  dx\n⎮ ╲╱   x    \n⌡           '
    >>> print(_)
    ⌠
    ⎮     ___
    ⎮    ╱ 1
    ⎮   ╱  ─  dx
    ⎮ ╲╱   x
    ⌡

.. note::

    Unicode pretty-printing is enabled by default in the `IPython`_
    terminal frontend.

You can pass ``use_unicode=False`` to use ASCII symbols.

    >>> print(pretty(expr, use_unicode=False))
       /
      |
      |     ___
      |    / 1
      |   /  -  dx
      | \/   x
      |
     /

:func:`~diofant.printing.pretty.pretty.pprint` prints the output to
the screen.

    >>> pprint(expr)
    ⌠
    ⎮     ___
    ⎮    ╱ 1
    ⎮   ╱  ─  dx
    ⎮ ╲╱   x
    ⌡

LaTeX
=====

To get the `\LaTeX` form of an expression, use
:func:`~diofant.printing.latex.latex`.

    >>> print(latex(expr))
    \int \sqrt{\frac{1}{x}}\, dx

Dot
===

:func:`~diofant.printing.dot.dotprint` function prints output to dot
format, which can be rendered with `Graphviz
<http://www.graphviz.org/>`_:

.. graphviz::

    digraph{

    # Graph style
    "bgcolor"="transparent"
    "ordering"="out"
    "rankdir"="TD"

    #########
    # Nodes #
    #########

    "Integral(Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2)), Tuple(Symbol('x')))_()" ["color"="black", "label"="Integral", "shape"="ellipse"];
    "Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2))_(0,)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Pow(Symbol('x'), Integer(-1))_(0, 0)" ["color"="black", "label"="Pow", "shape"="ellipse"];
    "Symbol('x')_(0, 0, 0)" ["color"="black", "label"="x", "shape"="ellipse"];
    "Integer(-1)_(0, 0, 1)" ["color"="black", "label"="-1", "shape"="ellipse"];
    "Rational(1, 2)_(0, 1)" ["color"="black", "label"="1/2", "shape"="ellipse"];
    "Tuple(Symbol('x'))_(1,)" ["color"="blue", "label"="Tuple", "shape"="ellipse"];
    "Symbol('x')_(1, 0)" ["color"="black", "label"="x", "shape"="ellipse"];

    #########
    # Edges #
    #########

    "Integral(Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2)), Tuple(Symbol('x')))_()" -> "Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2))_(0,)";
    "Integral(Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2)), Tuple(Symbol('x')))_()" -> "Tuple(Symbol('x'))_(1,)";
    "Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2))_(0,)" -> "Pow(Symbol('x'), Integer(-1))_(0, 0)";
    "Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2))_(0,)" -> "Rational(1, 2)_(0, 1)";
    "Pow(Symbol('x'), Integer(-1))_(0, 0)" -> "Symbol('x')_(0, 0, 0)";
    "Pow(Symbol('x'), Integer(-1))_(0, 0)" -> "Integer(-1)_(0, 0, 1)";
    "Tuple(Symbol('x'))_(1,)" -> "Symbol('x')_(1, 0)";
    }

.. _IPython: https://ipython.readthedocs.io/en/stable/
