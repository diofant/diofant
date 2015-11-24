Parsing input
=============

Parsing Functions Reference
---------------------------

.. autofunction:: sympy.parsing.sympy_parser.parse_expr

.. autofunction:: sympy.parsing.sympy_parser.stringify_expr

.. autofunction:: sympy.parsing.sympy_parser.eval_expr

.. autofunction:: sympy.parsing.maxima.parse_maxima

.. autofunction:: sympy.parsing.mathematica.mathematica

Parsing Transformations Reference
---------------------------------

A transformation is a function that accepts the arguments ``tokens,
local_dict, global_dict`` and returns a list of transformed tokens. They can
be used by passing a list of functions to
:py:func:`~sympy.parsing.sympy_parser.parse_expr` and are
applied in the order given.

.. autodata:: sympy.parsing.sympy_parser.standard_transformations

.. autofunction:: sympy.parsing.sympy_parser.split_symbols

.. autofunction:: sympy.parsing.sympy_parser.split_symbols_custom

.. autofunction:: sympy.parsing.sympy_parser.implicit_multiplication

.. autofunction:: sympy.parsing.sympy_parser.implicit_application

.. autofunction:: sympy.parsing.sympy_parser.function_exponentiation

.. autofunction:: sympy.parsing.sympy_parser.implicit_multiplication_application

.. autofunction:: sympy.parsing.sympy_parser.rationalize

.. autofunction:: sympy.parsing.sympy_parser.convert_xor

These are included in
:data:``sympy.parsing.sympy_parser.standard_transformations`` and generally
don't need to be manually added by the user.

.. autofunction:: sympy.parsing.sympy_parser.auto_symbol

.. autofunction:: sympy.parsing.sympy_parser.auto_number
