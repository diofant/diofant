"""
Helper module for setting up interactive Diofant sessions.

AST transformations, provided here, `could be used
<http://ipython.readthedocs.io/en/stable/config/inputtransforms.html>`_
in IPython to reduce boilerplate while interacting with Diofant
due to the Python language syntax.
"""

from .printing import init_printing

init_printing()
