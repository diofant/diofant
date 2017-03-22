============
SymPy 0.5.10
============

4 Jan 2008

* ``view`` renamed to ``preview``, ``pngview``, ``pdfview``, ``dviview`` added.
* Latex printer was rewritten, ``preview`` uses builtin ``pyglet``.
* Square root denesting implemented.
* Parser of simple Mathematica expressions added.
* TeXmacs interface written.
* Some integration fixes.
* Line width in 2D plotting can be specified.
* README was updated.
* ``pyglet`` and ``mpmath`` were updated and moved to ``sympy/thirdparty``
* All ``sys.path`` hacks were moved to just 2 places.
* SymPy objects should work in ``numpy`` arrays now.
* Hand written ``sympify()`` parser was rewritten and simplified using Python AST.
