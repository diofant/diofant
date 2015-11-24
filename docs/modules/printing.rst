Printing System
===============

See the :ref:`tutorial-printing` section in Tutorial for introduction into
printing.

This guide documents the printing system in SymPy and how it works
internally.

Printer Class
-------------

.. automodule:: sympy.printing.printer

The main class responsible for printing is ``Printer`` (see also its
`source code
<https://github.com/sympy/sympy/blob/master/sympy/printing/printer.py>`_):

.. autoclass:: Printer
    :members: doprint, _print, set_global_settings, order

    .. autoattribute:: Printer.printmethod


PrettyPrinter Class
-------------------

The pretty printing subsystem is implemented in ``sympy.printing.pretty.pretty``
by the ``PrettyPrinter`` class deriving from ``Printer``. It relies on
the modules ``sympy.printing.pretty.stringPict``, and
``sympy.printing.pretty.pretty_symbology`` for rendering nice-looking
formulas.

The module ``stringPict`` provides a base class ``stringPict`` and a derived
class ``prettyForm`` that ease the creation and manipulation of formulas
that span across multiple lines.

The module ``pretty_symbology`` provides primitives to construct 2D shapes
(hline, vline, etc) together with a technique to use unicode automatically
when possible.

.. module:: sympy.printing.pretty.pretty

.. autoclass:: PrettyPrinter
   :members: _use_unicode, doprint

   .. autoattribute:: PrettyPrinter.printmethod

.. autofunction:: pretty
.. autofunction:: pretty_print

CCodePrinter
------------

.. module:: sympy.printing.ccode

This class implements C code printing (i.e. it converts Python expressions
to strings of C code).

Usage::

    >>> from sympy.printing import print_ccode
    >>> from sympy.functions import sin, cos, Abs
    >>> from sympy.abc import x
    >>> print_ccode(sin(x)**2 + cos(x)**2)
    pow(sin(x), 2) + pow(cos(x), 2)
    >>> print_ccode(2*x + cos(x), assign_to="result")
    result = 2*x + cos(x);
    >>> print_ccode(Abs(x**2))
    fabs(pow(x, 2))

.. autodata:: sympy.printing.ccode.known_functions

.. autoclass:: sympy.printing.ccode.CCodePrinter
   :members:

   .. autoattribute:: CCodePrinter.printmethod


.. autofunction:: sympy.printing.ccode.ccode

.. autofunction:: sympy.printing.ccode.print_ccode

Fortran Printing
----------------

The ``fcode`` function translates a sympy expression into Fortran code. The main
purpose is to take away the burden of manually translating long mathematical
expressions. Therefore the resulting expression should also require no (or
very little) manual tweaking to make it compilable. The optional arguments
of ``fcode`` can be used to fine-tune the behavior of ``fcode`` in such a way
that manual changes in the result are no longer needed.

.. module:: sympy.printing.fcode
.. autofunction:: fcode
.. autofunction:: print_fcode
.. autoclass:: FCodePrinter
   :members:

   .. autoattribute:: FCodePrinter.printmethod


Two basic examples:

    >>> from sympy import *
    >>> x = symbols("x")
    >>> fcode(sqrt(1-x**2))
    '      sqrt(-x**2 + 1)'
    >>> fcode((3 + 4*I)/(1 - conjugate(x)))
    '      (cmplx(3,4))/(-conjg(x) + 1)'

An example where line wrapping is required:

    >>> expr = sqrt(1-x**2).series(x,n=20).removeO()
    >>> print(fcode(expr))
          -715.0d0/65536.0d0*x**18 - 429.0d0/32768.0d0*x**16 - 33.0d0/
         @ 2048.0d0*x**14 - 21.0d0/1024.0d0*x**12 - 7.0d0/256.0d0*x**10 -
         @ 5.0d0/128.0d0*x**8 - 1.0d0/16.0d0*x**6 - 1.0d0/8.0d0*x**4 - 1.0d0
         @ /2.0d0*x**2 + 1

In case of line wrapping, it is handy to include the assignment so that lines
are wrapped properly when the assignment part is added.

    >>> print(fcode(expr, assign_to="var"))
          var = -715.0d0/65536.0d0*x**18 - 429.0d0/32768.0d0*x**16 - 33.0d0/
         @ 2048.0d0*x**14 - 21.0d0/1024.0d0*x**12 - 7.0d0/256.0d0*x**10 -
         @ 5.0d0/128.0d0*x**8 - 1.0d0/16.0d0*x**6 - 1.0d0/8.0d0*x**4 - 1.0d0
         @ /2.0d0*x**2 + 1

For piecewise functions, the ``assign_to`` option is mandatory:

    >>> print(fcode(Piecewise((x,x<1),(x**2,True)), assign_to="var"))
          if (x < 1) then
            var = x
          else
            var = x**2
          end if

Note that by default only top-level piecewise functions are supported due to
the lack of a conditional operator in Fortran 77. Inline conditionals can be
supported using the ``merge`` function introduced in Fortran 95 by setting of
the kwarg ``standard=95``:

    >>> print(fcode(Piecewise((x,x<1),(x**2,True)), standard=95))
          merge(x, x**2, x < 1)

Loops are generated if there are Indexed objects in the expression. This
also requires use of the assign_to option.

    >>> A, B = map(IndexedBase, ['A', 'B'])
    >>> m = Symbol('m', integer=True)
    >>> i = Idx('i', m)
    >>> print(fcode(2*B[i], assign_to=A[i]))
        do i = 1, m
            A(i) = 2*B(i)
        end do

Repeated indices in an expression with Indexed objects are interpreted as
summation. For instance, code for the trace of a matrix can be generated
with

    >>> print(fcode(A[i, i], assign_to=x))
          x = 0
          do i = 1, m
              x = x + A(i, i)
          end do

By default, number symbols such as ``pi`` and ``E`` are detected and defined as
Fortran parameters. The precision of the constants can be tuned with the
precision argument. Parameter definitions are easily avoided using the ``N``
function.

    >>> print(fcode(x - pi**2 - E))
          parameter (E = 2.71828182845905d0)
          parameter (pi = 3.14159265358979d0)
          x - pi**2 - E
    >>> print(fcode(x - pi**2 - E, precision=25))
          parameter (E = 2.718281828459045235360287d0)
          parameter (pi = 3.141592653589793238462643d0)
          x - pi**2 - E
    >>> print(fcode(N(x - pi**2, 25)))
          x - 9.869604401089358618834491d0

When some functions are not part of the Fortran standard, it might be desirable
to introduce the names of user-defined functions in the Fortran expression.

    >>> print(fcode(1 - gamma(x)**2, user_functions={'gamma': 'mygamma'}))
          -mygamma(x)**2 + 1

However, when the user_functions argument is not provided, ``fcode`` attempts to
use a reasonable default and adds a comment to inform the user of the issue.

    >>> print(fcode(1 - gamma(x)**2))
    C     Not supported in Fortran:
    C     gamma
          -gamma(x)**2 + 1

By default the output is human readable code, ready for copy and paste. With the
option ``human=False``, the return value is suitable for post-processing with
source code generators that write routines with multiple instructions. The
return value is a three-tuple containing: (i) a set of number symbols that must
be defined as 'Fortran parameters', (ii) a list functions that can not be
translated in pure Fortran and (iii) a string of Fortran code. A few examples:

    >>> fcode(1 - gamma(x)**2, human=False) == (set(), {gamma(x)}, '      -gamma(x)**2 + 1')
    True
    >>> fcode(1 - sin(x)**2, human=False) == (set(), set(), '      -sin(x)**2 + 1')
    True
    >>> fcode(x - pi**2, human=False) == ({(pi, '3.14159265358979d0')}, set(), '      x - pi**2')
    True

Mathematica code printing
-------------------------

.. module:: sympy.printing.mathematica

.. autodata:: sympy.printing.mathematica.known_functions

.. autoclass:: sympy.printing.mathematica.MCodePrinter
   :members:

   .. autoattribute:: MCodePrinter.printmethod

.. autofunction:: sympy.printing.mathematica.mathematica_code

LambdaPrinter
-------------

.. module:: sympy.printing.lambdarepr

This classes implements printing to strings that can be used by the
:py:func:`sympy.utilities.lambdify.lambdify` function.

.. autoclass:: LambdaPrinter

   .. autoattribute:: LambdaPrinter.printmethod


.. autofunction:: lambdarepr

LatexPrinter
------------

.. module:: sympy.printing.latex

This class implements LaTeX printing. See ``sympy.printing.latex``.

.. autodata:: accepted_latex_functions

.. autoclass:: LatexPrinter
   :members:

   .. autoattribute:: LatexPrinter.printmethod

.. autofunction:: latex

.. autofunction:: print_latex

MathMLPrinter
-------------

.. module:: sympy.printing.mathml

This class is responsible for MathML printing. See ``sympy.printing.mathml``.

More info on mathml content: http://www.w3.org/TR/MathML2/chapter4.html

.. autoclass:: MathMLPrinter
   :members:

   .. autoattribute:: MathMLPrinter.printmethod

.. autofunction:: mathml

.. autofunction:: print_mathml

PythonPrinter
-------------

.. module:: sympy.printing.python

This class implements Python printing. Usage::

    >>> from sympy import print_python, sin
    >>> from sympy.abc import x

    >>> print_python(5*x**3 + sin(x))
    x = Symbol('x')
    e = 5*x**3 + sin(x)

ReprPrinter
-----------

.. module:: sympy.printing.repr

This printer generates executable code. This code satisfies the identity
``eval(srepr(expr)) == expr``.

.. autoclass:: ReprPrinter
   :members:

   .. autoattribute:: ReprPrinter.printmethod

.. autofunction:: srepr

StrPrinter
----------

.. module:: sympy.printing.str

This module generates readable representations of SymPy expressions.

.. autoclass:: StrPrinter
   :members: parenthesize, stringify, emptyPrinter

   .. autoattribute:: StrPrinter.printmethod

.. autofunction:: sstrrepr

Tree Printing
-------------

.. module:: sympy.printing.tree

The functions in this module create a representation of an expression as a
tree.

.. autofunction:: pprint_nodes

.. autofunction:: print_node

.. autofunction:: tree

.. autofunction:: print_tree

Implementation - Helper Classes/Functions
-----------------------------------------

.. module:: sympy.printing.conventions

.. autofunction:: split_super_sub

CodePrinter
+++++++++++

.. module:: sympy.printing.codeprinter

This class is a base class for other classes that implement code-printing
functionality, and additionally lists a number of functions that cannot be
easily translated to C or Fortran.

.. autoclass:: sympy.printing.codeprinter.CodePrinter

   .. autoattribute:: CodePrinter.printmethod

.. autoexception:: sympy.printing.codeprinter.AssignmentError

Precedence
++++++++++

.. module:: sympy.printing.precedence

.. autodata:: PRECEDENCE

   Default precedence values for some basic types.

.. autodata:: PRECEDENCE_VALUES

   A dictionary assigning precedence values to certain classes. These values
   are treated like they were inherited, so not every single class has to be
   named here.

.. autodata:: PRECEDENCE_FUNCTIONS

   Sometimes it's not enough to assign a fixed precedence value to a
   class. Then a function can be inserted in this dictionary that takes an
   instance of this class as argument and returns the appropriate precedence
   value.

.. autofunction:: precedence

Pretty-Printing Implementation Helpers
--------------------------------------

.. module:: sympy.printing.pretty.pretty_symbology

.. autofunction:: U
.. autofunction:: pretty_use_unicode
.. autofunction:: pretty_try_use_unicode
.. autofunction:: xstr

The following two functions return the Unicode version of the inputted Greek
letter.

.. autofunction:: g
.. autofunction:: G
.. autodata:: greek_letters
.. autodata:: digit_2txt
.. autodata:: symb_2txt

The following functions return the Unicode subscript/superscript version of
the character.

.. autodata:: sub
.. autodata:: sup

The following functions return Unicode vertical objects.

.. autofunction:: xobj
.. autofunction:: vobj
.. autofunction:: hobj

The following constants are for rendering roots and fractions.

.. autodata:: root
.. autofunction:: VF
.. autodata:: frac

The following constants/functions are for rendering atoms and symbols.

.. autofunction:: xsym
.. autodata:: atoms_table
.. autofunction:: pretty_atom
.. autofunction:: pretty_symbol
.. autofunction:: annotated

.. automodule:: sympy.printing.pretty.stringpict

.. autoclass:: stringPict
   :members:

.. autoclass:: prettyForm
   :members:

dotprint
--------

.. autofunction:: sympy.printing.dot.dotprint
