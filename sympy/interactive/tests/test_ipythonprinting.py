"""Tests that the IPython printing module is properly loaded. """

import pytest

from sympy.core.compatibility import u
from sympy.interactive.session import init_ipython_session
from sympy.external import import_module

ipython = import_module("IPython", min_module_version="2.3.0")

# disable tests if ipython is not present
if not ipython:
    disabled = True


def test_ipythonprinting():
    # Initialize and setup IPython session
    app = init_ipython_session()
    app.run_cell("ip = get_ipython()")
    app.run_cell("inst = ip.instance()")
    app.run_cell("format = inst.display_formatter.format")
    app.run_cell("from sympy import Symbol")

    # Printing by default
    app.run_cell("a = format(Symbol('pi'))")
    app.run_cell("a2 = format(Symbol('pi')**2)")
    assert app.user_ns['a'][0]['text/plain'] in (u('\N{GREEK SMALL LETTER PI}'), 'pi')
    assert app.user_ns['a2'][0]['text/plain'] in (u(' 2\n\N{GREEK SMALL LETTER PI} '), '  2\npi ')

    # Use different printing setup
    app.run_cell("from sympy import init_printing")
    app.run_cell("init_printing(use_unicode=False, pretty_print=False)")
    app.run_cell("a = format(Symbol('pi'))")
    app.run_cell("a2 = format(Symbol('pi')**2)")
    assert app.user_ns['a'][0]['text/plain'] == "pi"
    assert app.user_ns['a2'][0]['text/plain'] == "pi**2"


def test_print_builtin_option():
    # Initialize and setup IPython session
    app = init_ipython_session()
    app.run_cell("ip = get_ipython()")
    app.run_cell("inst = ip.instance()")
    app.run_cell("format = inst.display_formatter.format")
    app.run_cell("from sympy import Symbol")
    app.run_cell("from sympy import init_printing")

    app.run_cell("a = format({Symbol('pi'): 3.14, Symbol('n_i'): 3})")
    text = app.user_ns['a'][0]['text/plain']
    pytest.raises(KeyError, lambda: app.user_ns['a'][0]['text/latex'])
    # Note : In Python 3 the text is unicode, but in 2 it is a string.
    # XXX: How can we make this ignore the terminal width? This test fails if
    # the terminal is too narrow.
    assert text in ("{pi: 3.14, n_i: 3}",
                    u('{n\N{LATIN SUBSCRIPT SMALL LETTER I}: 3, \N{GREEK SMALL LETTER PI}: 3.14}'),
                    "{n_i: 3, pi: 3.14}",
                    u('{\N{GREEK SMALL LETTER PI}: 3.14, n\N{LATIN SUBSCRIPT SMALL LETTER I}: 3}'))

    # If we enable the default printing, then the dictionary's should render
    # as a LaTeX version of the whole dict: ${\pi: 3.14, n_i: 3}$
    app.run_cell("inst.display_formatter.formatters['text/latex'].enabled = True")
    app.run_cell("init_printing(use_latex=True)")
    app.run_cell("a = format({Symbol('pi'): 3.14, Symbol('n_i'): 3})")
    text = app.user_ns['a'][0]['text/plain']
    latex = app.user_ns['a'][0]['text/latex']
    assert text in ("{pi: 3.14, n_i: 3}",
                    u('{n\N{LATIN SUBSCRIPT SMALL LETTER I}: 3, \N{GREEK SMALL LETTER PI}: 3.14}'),
                    "{n_i: 3, pi: 3.14}",
                    u('{\N{GREEK SMALL LETTER PI}: 3.14, n\N{LATIN SUBSCRIPT SMALL LETTER I}: 3}'))
    assert latex == r'$$\left \{ n_{i} : 3, \quad \pi : 3.14\right \}$$'

    app.run_cell("inst.display_formatter.formatters['text/latex'].enabled = True")
    app.run_cell("init_printing(use_latex=True, print_builtin=False)")
    app.run_cell("a = format({Symbol('pi'): 3.14, Symbol('n_i'): 3})")
    text = app.user_ns['a'][0]['text/plain']
    pytest.raises(KeyError, lambda: app.user_ns['a'][0]['text/latex'])
    # Note : In Python 3 the text is unicode, but in 2 it is a string.
    assert text in ("{pi: 3.14, n_i: 3}", "{n_i: 3, pi: 3.14}")


def test_matplotlib_bad_latex():
    # Initialize and setup IPython session
    app = init_ipython_session()
    app.run_cell("import IPython")
    app.run_cell("ip = get_ipython()")
    app.run_cell("inst = ip.instance()")
    app.run_cell("format = inst.display_formatter.format")
    app.run_cell("from sympy import init_printing, Matrix")
    app.run_cell("init_printing(use_latex='matplotlib')")

    # The png formatter is not enabled by default in this context
    app.run_cell("inst.display_formatter.formatters['image/png'].enabled = True")

    # Make sure no warnings are raised by IPython
    app.run_cell("import warnings")
    app.run_cell("warnings.simplefilter('error', IPython.core.formatters.FormatterWarning)")

    # This should not raise an exception
    app.run_cell("a = format(Matrix([1, 2, 3]))")
