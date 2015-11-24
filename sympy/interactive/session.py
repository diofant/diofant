"""Tools for setting up interactive sessions. """

from __future__ import print_function, division

from distutils.version import LooseVersion as V

from sympy.core.compatibility import range
from sympy.external import import_module
from sympy.interactive.printing import init_printing

preexec_source = """\
from __future__ import division
from sympy import *
x, y, z, t = symbols('x y z t')
k, m, n = symbols('k m n', integer=True)
f, g, h = symbols('f g h', cls=Function)
init_printing()
"""

verbose_message = """\
These commands were executed:
%(source)s
Documentation can be found at http://docs.sympy.org/%(version)s
"""

no_ipython = """\
Couldn't locate IPython. Having IPython installed is greatly recommended.
See http://ipython.scipy.org for more details. If you use Debian/Ubuntu,
just install first the 'ipython' package.
"""


def _make_message(ipython=True, quiet=False, source=None):
    """Create a banner for an interactive session. """
    from sympy import __version__ as sympy_version
    from sympy.polys.domains import GROUND_TYPES
    from sympy.utilities.misc import ARCH
    from sympy import SYMPY_DEBUG

    import sys
    import os

    python_version = "%d.%d.%d" % sys.version_info[:3]

    if ipython:
        shell_name = "IPython"
    else:
        shell_name = "Python"

    info = ['ground types: %s' % GROUND_TYPES]

    cache = os.getenv('SYMPY_USE_CACHE')

    if cache is not None and cache.lower() == 'no':
        info.append('cache: off')

    if SYMPY_DEBUG:
        info.append('debugging: on')

    args = shell_name, sympy_version, python_version, ARCH, ', '.join(info)
    message = "%s console for SymPy %s (Python %s-%s) (%s)\n" % args

    if not quiet:
        if source is None:
            source = preexec_source

        _source = ""

        for line in source.split('\n')[:-1]:
            if not line:
                _source += '\n'
            else:
                _source += '>>> ' + line + '\n'

        doc_version = sympy_version
        if 'dev' in doc_version:
            doc_version = "dev"
        else:
            doc_version = "%s.%s.%s/" % tuple(doc_version.split('.')[:3])

        message += '\n' + verbose_message % {'source': _source,
                                             'version': doc_version}

    return message


def int_to_Integer(s):
    """
    Wrap integer literals with Integer.

    This is based on the decistmt example from
    http://docs.python.org/library/tokenize.html.

    Only integer literals are converted.  Float literals are left alone.

    Examples
    ========

    >>> from __future__ import division
    >>> from sympy.interactive.session import int_to_Integer
    >>> from sympy import Integer
    >>> s = '1.2 + 1/2 - 0x12 + a1'
    >>> int_to_Integer(s)
    '1.2 +Integer (1 )/Integer (2 )-Integer (0x12 )+a1 '
    >>> s = 'print (1/2)'
    >>> int_to_Integer(s)
    'print (Integer (1 )/Integer (2 ))'
    >>> exec(s)
    0.5
    >>> exec(int_to_Integer(s))
    1/2
    """
    from tokenize import generate_tokens, untokenize, NUMBER, NAME, OP
    from sympy.core.compatibility import StringIO

    def _is_int(num):
        """
        Returns true if string value num (with token NUMBER) represents an integer.
        """
        # XXX: Is there something in the standard library that will do this?
        if '.' in num or 'j' in num.lower() or 'e' in num.lower():
            return False
        return True

    result = []
    g = generate_tokens(StringIO(s).readline)   # tokenize the string
    for toknum, tokval, _, _, _ in g:
        if toknum == NUMBER and _is_int(tokval):  # replace NUMBER tokens
            result.extend([
                (NAME, 'Integer'),
                (OP, '('),
                (NUMBER, tokval),
                (OP, ')')
            ])
        else:
            result.append((toknum, tokval))
    return untokenize(result)


def enable_automatic_int_sympification(app):
    """
    Allow IPython to automatically convert integer literals to Integer.
    """
    hasshell = hasattr(app, 'shell')

    import ast
    if hasshell:
        old_run_cell = app.shell.run_cell
    else:
        old_run_cell = app.run_cell

    def my_run_cell(cell, *args, **kwargs):
        try:
            # Check the cell for syntax errors.  This way, the syntax error
            # will show the original input, not the transformed input.  The
            # downside here is that IPython magic like %timeit will not work
            # with transformed input (but on the other hand, IPython magic
            # that doesn't expect transformed input will continue to work).
            ast.parse(cell)
        except SyntaxError:
            pass
        else:
            cell = int_to_Integer(cell)
        old_run_cell(cell, *args, **kwargs)

    if hasshell:
        app.shell.run_cell = my_run_cell
    else:
        app.run_cell = my_run_cell


def enable_automatic_symbols(app):
    """Allow IPython to automatially create symbols. """
    # XXX: This should perhaps use tokenize, like int_to_Integer() above.
    # This would avoid re-executing the code, which can lead to subtle
    # issues.  For example:
    #
    # In [1]: a = 1
    #
    # In [2]: for i in range(10):
    #    ...:     a += 1
    #    ...:
    #
    # In [3]: a
    # Out[3]: 11
    #
    # In [4]: a = 1
    #
    # In [5]: for i in range(10):
    #    ...:     a += 1
    #    ...:     print b
    #    ...:
    # b
    # b
    # b
    # b
    # b
    # b
    # b
    # b
    # b
    # b
    #
    # In [6]: a
    # Out[6]: 12
    #
    # Note how the for loop is executed again because `b` was not defined, but `a`
    # was already incremented once, so the result is that it is incremented
    # multiple times.

    import re
    re_nameerror = re.compile(
        "name '(?P<symbol>[A-Za-z_][A-Za-z0-9_]*)' is not defined")

    def _handler(self, etype, value, tb, tb_offset=None):
        """Handle :exc:`NameError` exception and allow injection of missing symbols. """
        if etype is NameError and tb.tb_next and not tb.tb_next.tb_next:
            match = re_nameerror.match(str(value))

            if match is not None:
                # XXX: Make sure Symbol is in scope. Otherwise you'll get infinite recursion.
                self.run_cell("%(symbol)s = Symbol('%(symbol)s')" %
                    {'symbol': match.group("symbol")}, store_history=False)

                try:
                    code = self.user_ns['In'][-1]
                except (KeyError, IndexError):
                    pass
                else:
                    self.run_cell(code, store_history=False)
                    return
                finally:
                    self.run_cell("del %s" % match.group("symbol"),
                                  store_history=False)

        stb = self.InteractiveTB.structured_traceback(
            etype, value, tb, tb_offset=tb_offset)
        self._showtraceback(etype, value, stb)

    if hasattr(app, 'shell'):
        app.shell.set_custom_exc((NameError,), _handler)
    else:
        # This was restructured in IPython 0.13
        app.set_custom_exc((NameError,), _handler)


def init_ipython_session(argv=[], auto_symbols=False, auto_int_to_Integer=False):
    """Construct new IPython session. """
    import IPython

    if V(IPython.__version__) >= '0.11':
        # use an app to parse the command line, and init config
        # IPython 1.0 deprecates the frontend module, so we import directly
        # from the terminal module to prevent a deprecation message from being
        # shown.
        if V(IPython.__version__) >= '1.0':
            from IPython.terminal import ipapp
        else:
            from IPython.frontend.terminal import ipapp
        app = ipapp.TerminalIPythonApp()

        # don't draw IPython banner during initialization:
        app.display_banner = False
        app.initialize(argv)

        if auto_symbols:
            readline = import_module("readline")
            if readline:
                enable_automatic_symbols(app)
        if auto_int_to_Integer:
            enable_automatic_int_sympification(app)

        return app.shell
    else:
        from IPython.Shell import make_IPython
        return make_IPython(argv)


def init_python_session():
    """Construct new Python session. """
    from code import InteractiveConsole

    class SymPyConsole(InteractiveConsole):
        """An interactive console with readline support. """

        def __init__(self):
            InteractiveConsole.__init__(self)

            try:
                import readline
            except ImportError:
                pass
            else:
                import os
                import atexit

                readline.parse_and_bind('tab: complete')

                if hasattr(readline, 'read_history_file'):
                    history = os.path.expanduser('~/.sympy-history')

                    try:
                        readline.read_history_file(history)
                    except IOError:
                        pass

                    atexit.register(readline.write_history_file, history)

    return SymPyConsole()


def init_session(ipython=None, pretty_print=True, order=None,
        use_unicode=None, use_latex=None, quiet=False, auto_symbols=False,
        auto_int_to_Integer=False, argv=[]):
    """
    Initialize an embedded IPython or Python session. The IPython session is
    initiated with the --pylab option, without the numpy imports, so that
    matplotlib plotting can be interactive.

    Parameters
    ==========

    pretty_print: boolean
        If True, use pretty_print to stringify;
        if False, use sstrrepr to stringify.
    order: string or None
        There are a few different settings for this parameter:
        lex (default), which is lexographic order;
        grlex, which is graded lexographic order;
        grevlex, which is reversed graded lexographic order;
        old, which is used for compatibility reasons and for long expressions;
        None, which sets it to lex.
    use_unicode: boolean or None
        If True, use unicode characters;
        if False, do not use unicode characters.
    use_latex: boolean or None
        If True, use latex rendering if IPython GUI's;
        if False, do not use latex rendering.
    quiet: boolean
        If True, init_session will not print messages regarding its status;
        if False, init_session will print messages regarding its status.
    auto_symbols: boolean
        If True, IPython will automatically create symbols for you.
        If False, it will not.
        The default is False.
    auto_int_to_Integer: boolean
        If True, IPython will automatically wrap int literals with Integer, so
        that things like 1/2 give Rational(1, 2).
        If False, it will not.
        The default is False.
    ipython: boolean or None
        If True, printing will initialize for an IPython console;
        if False, printing will initialize for a normal console;
        The default is None, which automatically determines whether we are in
        an ipython instance or not.
    argv: list of arguments for IPython

    See Also
    ========

    sympy.interactive.printing.init_printing: for examples and the rest of the parameters.


    Examples
    ========

    >>> from sympy import init_session, Symbol, sin, sqrt
    >>> sin(x) #doctest: +SKIP
    NameError: name 'x' is not defined
    >>> init_session() #doctest: +SKIP
    >>> sin(x) #doctest: +SKIP
    sin(x)
    >>> sqrt(5) #doctest: +SKIP
      ___
    \/ 5
    >>> init_session(pretty_print=False) #doctest: +SKIP
    >>> sqrt(5) #doctest: +SKIP
    sqrt(5)
    >>> y + x + y**2 + x**2 #doctest: +SKIP
    x**2 + x + y**2 + y
    >>> init_session(order='grlex') #doctest: +SKIP
    >>> y + x + y**2 + x**2 #doctest: +SKIP
    x**2 + y**2 + x + y
    >>> init_session(order='grevlex') #doctest: +SKIP
    >>> y * x**2 + x * y**2 #doctest: +SKIP
    x**2*y + x*y**2
    >>> init_session(order='old') #doctest: +SKIP
    >>> x**2 + y**2 + x + y #doctest: +SKIP
    x + y + x**2 + y**2
    >>> theta = Symbol('theta') #doctest: +SKIP
    >>> theta #doctest: +SKIP
    theta
    >>> init_session(use_unicode=True) #doctest: +SKIP
    >>> theta # doctest: +SKIP
    \u03b8
    """
    import sys

    in_ipython = False

    if ipython is not False:
        try:
            import IPython
        except ImportError:
            if ipython is True:
                raise RuntimeError("IPython is not available on this system")
            ip = None
        else:
            if V(IPython.__version__) >= '0.11':
                try:
                    ip = get_ipython()
                except NameError:
                    ip = None
            else:
                ip = IPython.ipapi.get()
                if ip:
                    ip = ip.IP
        in_ipython = bool(ip)
        if ipython is None:
            ipython = in_ipython

    if ipython is False:
        ip = init_python_session()
        mainloop = ip.interact
    else:
        if ip is None:
            ip = init_ipython_session(argv=argv, auto_symbols=auto_symbols,
                auto_int_to_Integer=auto_int_to_Integer)

        if V(IPython.__version__) >= '0.11':
            # runsource is gone, use run_cell instead, which doesn't
            # take a symbol arg.  The second arg is `store_history`,
            # and False means don't add the line to IPython's history.
            ip.runsource = lambda src, symbol='exec': ip.run_cell(src, False)

            # Enable interactive plotting using pylab.
            try:
                ip.enable_pylab(import_all=False)
            except Exception:
                # Causes an import error if matplotlib is not installed.
                # Causes other errors (depending on the backend) if there
                # is no display, or if there is some problem in the
                # backend, so we have a bare "except Exception" here
                pass
        if not in_ipython:
            mainloop = ip.mainloop

    readline = import_module("readline")
    if auto_symbols and (not ipython or V(IPython.__version__) < '0.11' or not readline):
        raise RuntimeError("automatic construction of symbols is possible only in IPython 0.11 or above with readline support")
    if auto_int_to_Integer and (not ipython or V(IPython.__version__) < '0.11'):
        raise RuntimeError("automatic int to Integer transformation is possible only in IPython 0.11 or above")

    _preexec_source = preexec_source

    ip.runsource(_preexec_source, symbol='exec')
    init_printing(pretty_print=pretty_print, order=order,
        use_unicode=use_unicode, use_latex=use_latex, ip=ip)

    message = _make_message(ipython, quiet, _preexec_source)

    if not in_ipython:
        mainloop(message)
        sys.exit('Exiting ...')
    else:
        ip.write(message)
        import atexit
        atexit.register(lambda ip: ip.write("Exiting ...\n"), ip)
