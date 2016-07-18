.. _tutorial-printing:

==========
 Printing
==========

As we have already seen, Diofant can pretty print its output using Unicode
characters.  This is a short introduction to the most common printing options
available in Diofant.

Printers
========

There are several printers available in Diofant.  The most common ones are

- str
- repr
- ASCII pretty printer
- Unicode pretty printer
- LaTeX
- MathML
- Dot

In addition to these, there are also "printers" that can output Diofant objects
to code, such as C, Fortran, Javascript, Theano, and Python.  These are not
discussed in this tutorial.

Setting up Pretty Printing
==========================

Pretty printing is enabled automatically for IPython's sessions, best
printer is used by default (i.e. `\LaTeX` in the IPython QTConsole or
notebooks and Unicode pretty printer in the IPython console).

  .. image:: ../pics/ipythonconsole.png

If you want manually configure pretty printing, please use the
``init_printing()`` function.  Without arguments it will just automatically
enable the best printer available in your environment.

To explicitly not use `\LaTeX`, pass ``use_latex=False`` to
``init_printing()`` or ``init_session()``.  To explicitly not use Unicode,
pass ``use_unicode=False``.

For convinience, the ``init_session()`` function is available, it will
automatically import everything in Diofant, create some common Symbols, setup
plotting, and run ``init_printing()``.

Printing Functions
==================

In addition to automatic printing, you can explicitly use any one of the
printers by calling the appropriate function.

str
---

To get a string form of an expression, use ``str(expr)``.  This is also the
form that is produced by ``print(expr)``.  String forms are designed to be
easy to read, but in a form that is correct Python syntax so that it can be
copied and pasted.  The ``str()`` form of an expression will usually look
exactly the same as the expression as you would enter it.

    >>> from diofant import *
    >>> x, y, z = symbols('x y z')
    >>> str(Integral(sqrt(1/x), x))
    'Integral(sqrt(1/x), x)'
    >>> print(Integral(sqrt(1/x), x))
    Integral(sqrt(1/x), x)

repr
----

The repr form of an expression is designed to show the exact form of an
expression.  It will be discussed more in the :ref:`tutorial-manipulation`
section.  To get it, use ``repr()``.

    >>> repr(Integral(sqrt(1/x), x))
    "Integral(Pow(Pow(Symbol('x'), Integer(-1)), Rational(1, 2)), Tuple(Symbol('x')))"

The repr form is mostly useful for understanding how an expression is built
internally.


ASCII Pretty Printer
--------------------

The ASCII pretty printer is accessed from ``pprint()``.  If the terminal does
not support Unicode, the ASCII printer is used by default.  Otherwise, you
must pass ``use_unicode=False``.

    >>> pprint(Integral(sqrt(1/x), x), use_unicode=False)
      /
     |
     |     ___
     |    / 1
     |   /  -  dx
     | \/   x
     |
    /

``pprint()`` prints the output to the screen.  If you want the string form,
use ``pretty()``.

    >>> pretty(sqrt(1/x), use_unicode=False)
    '    ___\n   / 1 \n  /  - \n\\/   x '
    >>> print(pretty(Integral(sqrt(1/x), x), use_unicode=False))
      /
     |
     |     ___
     |    / 1
     |   /  -  dx
     | \/   x
     |
    /

Unicode Pretty Printer
----------------------

The Unicode pretty printer is also accessed from ``pprint()`` and
``pretty()``.  It the terminal supports Unicode, it is used automatically.  If
``pprint()`` is not able to detect that the terminal supports unicode, you can
pass ``use_unicode=True`` to force it to use Unicode.

    >>> pprint(Integral(sqrt(1/x), x), use_unicode=True)
    ⌠
    ⎮     ___
    ⎮    ╱ 1
    ⎮   ╱  ─  dx
    ⎮ ╲╱   x
    ⌡

.. _LaTeX:

`\LaTeX`
--------

To get the `\LaTeX` form of an expression, use ``latex()``.

    >>> print(latex(Integral(sqrt(1/x), x)))
    \int \sqrt{\frac{1}{x}}\, dx

The ``latex()`` function has many options to change the formatting of
different things.  See :py:meth:`its documentation
<diofant.printing.latex.latex>` for more details.

MathML
------

There is also a printer to MathML, called ``print_mathml()``.  It must be
imported from ``diofant.printing.mathml``.

    >>> from diofant.printing.mathml import print_mathml
    >>> print_mathml(Integral(sqrt(1/x), x))
    <apply>
        <int/>
        <bvar>
            <ci>x</ci>
        </bvar>
        <apply>
            <root/>
            <apply>
                <power/>
                <ci>x</ci>
                <cn>-1</cn>
            </apply>
        </apply>
    </apply>

``print_mathml()`` prints the output.  If you want the string, use the
function ``mathml()``.

Dot
---

The ``dotprint()`` function in ``diofant.printing.dot`` prints output to dot
format, which can be rendered with Graphviz.  See the
:ref:`tutorial-manipulation` section for some examples of the output of this
printer.
