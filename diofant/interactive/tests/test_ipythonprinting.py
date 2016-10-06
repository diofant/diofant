"""Tests that the IPython printing module is properly loaded. """

import pytest

from diofant.core import Symbol, Integer, Rational
from diofant.external import import_module

ipython = import_module("IPython", min_module_version="2.3.0")


def init_ipython_session():
    ip = ipython.get_ipython()

    if not ip:
        app = ipython.terminal.ipapp.TerminalIPythonApp()
        app.display_banner = False
        app.initialize([])
        ip = app.shell

    return ip


@pytest.mark.skipif(ipython is None, reason="no IPython")
def test_int_to_Integer():
    app = init_ipython_session()
    app.run_cell("ip = get_ipython()")
    app.run_cell("from diofant import Integer")
    app.run_cell("from diofant.interactive.session import IntegerWrapper")
    app.run_cell("ip.ast_transformers.append(IntegerWrapper())")

    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], Integer)

    app.run_cell("a = 1/2")
    assert isinstance(app.user_ns['a'], Rational)
    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], Integer)
    app.run_cell("a = int(1)")
    assert isinstance(app.user_ns['a'], int)
    app.run_cell("a = (1/\n2)")
    assert app.user_ns['a'] == Rational(1, 2)


@pytest.mark.skipif(ipython is None, reason="no IPython")
def test_automatic_symbols():
    app = init_ipython_session()
    app.run_cell("ip = get_ipython()")
    app.run_cell("from diofant.interactive.session import AutomaticSymbols")
    app.run_cell("ip.ast_transformers.append(AutomaticSymbols())")
    app.run_cell("from diofant import Symbol, factorial")

    symbol = "verylongsymbolname"
    assert symbol not in app.user_ns
    app.run_cell("a = %s" % symbol)
    app.run_cell("a = type(%s)" % symbol)
    assert app.user_ns['a'] == Symbol
    app.run_cell("%s = Symbol('%s')" % (symbol, symbol))
    assert symbol in app.user_ns

    # Check that built-in names aren't overridden
    app.run_cell("a = all == __builtin__.all")
    assert "all" not in app.user_ns
    assert app.user_ns['a'] is True

    # Check that diofant names aren't overridden
    app.run_cell("import diofant")
    app.run_cell("a = factorial == diofant.factorial")
    assert app.user_ns['a'] is True


@pytest.mark.skipif(ipython is None, reason="no IPython")
def test_printing():
    # Initialize and setup IPython session
    app = init_ipython_session()
    app.run_cell("ip = get_ipython()")
    app.run_cell("inst = ip.instance()")
    app.run_cell("format = inst.display_formatter.format")
    app.run_cell("from diofant import Symbol")
    app.run_cell("from diofant import init_printing")

    # Printing by default
    app.run_cell("a = format(Symbol('pi'))")
    app.run_cell("a2 = format(Symbol('pi')**2)")
    app.run_cell("init_printing()")
    assert app.user_ns['a'][0]['text/plain'] in ('\N{GREEK SMALL LETTER PI}', 'pi')
    assert app.user_ns['a2'][0]['text/plain'] in (' 2\n\N{GREEK SMALL LETTER PI} ', '  2\npi ')

    # Use different printing setup
    app.run_cell("init_printing(use_unicode=False)")
    app.run_cell("a = format(Symbol('pi'))")
    app.run_cell("a2 = format(Symbol('pi')**2)")
    assert app.user_ns['a'][0]['text/plain'] == "pi"
    assert app.user_ns['a2'][0]['text/plain'] == "  2\npi "


@pytest.mark.skipif(ipython is None, reason="no IPython")
def test_dumbterm(monkeypatch):
    monkeypatch.setenv('TERM', '')
    app = init_ipython_session()
    app.run_cell("ip = get_ipython()")
    app.run_cell("inst = ip.instance()")
    app.run_cell("format = inst.display_formatter.format")
    app.run_cell("from diofant import Symbol")
    app.run_cell("from diofant import init_printing")
    app.run_cell("init_printing()")
    app.run_cell("a = format(Symbol('pi'))")
    assert app.user_ns['a'][0]['text/plain'] == "pi"


@pytest.mark.skipif(ipython is None, reason="no IPython")
def test_print_builtins():
    # Initialize and setup IPython session
    app = init_ipython_session()
    app.run_cell("ip = get_ipython()")
    app.run_cell("inst = ip.instance()")
    app.run_cell("format = inst.display_formatter.format")
    app.run_cell("from diofant import Symbol")
    app.run_cell("from diofant import init_printing, pretty, sstrrepr")

    app.run_cell("a = format(int(1))")
    pytest.raises(KeyError, lambda: app.user_ns['a'][0]['text/latex'])

    app.run_cell("init_printing()")
    app.run_cell("a = format({Symbol('pi'): 3.14, Symbol('n_i'): 3})")
    text = app.user_ns['a'][0]['text/plain']
    pytest.raises(KeyError, lambda: app.user_ns['a'][0]['text/latex'])
    # XXX: How can we make this ignore the terminal width? This test fails if
    # the terminal is too narrow.
    assert text in ('{n\N{LATIN SUBSCRIPT SMALL LETTER I}: 3, \N{GREEK SMALL LETTER PI}: 3.14}',
                    '{\N{GREEK SMALL LETTER PI}: 3.14, n\N{LATIN SUBSCRIPT SMALL LETTER I}: 3}')

    # If we enable the default printing, then the dictionary's should render
    # as a LaTeX version of the whole dict: ${\pi: 3.14, n_i: 3}$
    app.run_cell("init_printing()")
    app.run_cell("inst.display_formatter.formatters['text/latex'].enabled = True")
    app.run_cell("a = format({Symbol('pi'): 3.14, Symbol('n_i'): 3})")
    text = app.user_ns['a'][0]['text/plain']
    latex = app.user_ns['a'][0]['text/latex']
    assert text in ('{n\N{LATIN SUBSCRIPT SMALL LETTER I}: 3, \N{GREEK SMALL LETTER PI}: 3.14}',
                    '{\N{GREEK SMALL LETTER PI}: 3.14, n\N{LATIN SUBSCRIPT SMALL LETTER I}: 3}')
    assert latex == r'\begin{equation}\left \{ n_{i} : 3, \quad \pi : 3.14\right \}\end{equation}'

    app.run_cell("a = format([Symbol('x'), Symbol('y')])")
    latex = app.user_ns['a'][0]['text/latex']
    assert latex in (r'\begin{equation}\left [ x, \quad y\right ]\end{equation}',
                     r'\begin{equation}left [ y, \quad x\right ]\end{equation}')

    app.run_cell("a = format(False)")
    assert app.user_ns['a'][0]['text/plain'] == r'False'

    app.run_cell("init_printing(no_global=True)")
    app.run_cell("a = format(set([Symbol('pi')]))")
    assert app.user_ns['a'][0]['text/plain'] == '{\N{GREEK SMALL LETTER PI}}'

    app.run_cell("init_printing(use_unicode=False)")
    app.run_cell("a = format(Symbol('pi'))")
    assert app.user_ns['a'][0]['text/plain'] == 'pi'

    app.run_cell("init_printing()")
    app.run_cell("inst.display_formatter.formatters['text/latex'].enabled = False")
    app.run_cell("a = format({Symbol('pi'): 3.14, Symbol('n_i'): 3})")
    text = app.user_ns['a'][0]['text/plain']
    pytest.raises(KeyError, lambda: app.user_ns['a'][0]['text/latex'])


@pytest.mark.skipif(ipython is None, reason="no IPython")
def test_matplotlib_bad_latex():
    # Initialize and setup IPython session
    app = init_ipython_session()
    app.run_cell("import IPython")
    app.run_cell("ip = get_ipython()")
    app.run_cell("inst = ip.instance()")
    app.run_cell("format = inst.display_formatter.format")
    app.run_cell("from diofant import init_printing, Matrix")
    app.run_cell("init_printing()")

    # Make sure no warnings are raised by IPython
    app.run_cell("import warnings")
    app.run_cell("warnings.simplefilter('error', IPython.core.formatters.FormatterWarning)")

    # This should not raise an exception
    app.run_cell("a = format(Matrix([1, 2, 3]))")
