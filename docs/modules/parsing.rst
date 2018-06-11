Parsing
=======

.. module:: diofant.parsing

Parsing Functions Reference
---------------------------

.. autofunction:: diofant.parsing.sympy_parser.parse_expr

.. autofunction:: diofant.parsing.sympy_parser.stringify_expr

.. autofunction:: diofant.parsing.sympy_parser.eval_expr

.. autofunction:: diofant.parsing.maxima.parse_maxima

.. autofunction:: diofant.parsing.mathematica.mathematica

Parsing Transformations Reference
---------------------------------

A transformation is a function that accepts the arguments ``tokens,
local_dict, global_dict`` and returns a list of transformed tokens. They can
be used by passing a list of functions to
:py:func:`~diofant.parsing.sympy_parser.parse_expr` and are
applied in the order given.

.. autodata:: diofant.parsing.sympy_parser.standard_transformations

.. autofunction:: diofant.parsing.sympy_parser.split_symbols

.. autofunction:: diofant.parsing.sympy_parser.split_symbols_custom

.. autofunction:: diofant.parsing.sympy_parser.implicit_multiplication

.. autofunction:: diofant.parsing.sympy_parser.implicit_application

.. autofunction:: diofant.parsing.sympy_parser.function_exponentiation

.. autofunction:: diofant.parsing.sympy_parser.implicit_multiplication_application

.. autofunction:: diofant.parsing.sympy_parser.rationalize

.. autofunction:: diofant.parsing.sympy_parser.convert_xor

These are included in
:data:``diofant.parsing.sympy_parser.standard_transformations`` and generally
don't need to be manually added by the user.

.. autofunction:: diofant.parsing.sympy_parser.auto_symbol

.. autofunction:: diofant.parsing.sympy_parser.auto_number
