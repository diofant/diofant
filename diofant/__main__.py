"""
Python shell for Diofant.

This is just a normal Python shell (IPython shell if you have the
IPython package installed), that adds default imports and run
some initialization code.
"""

import argparse
import ast
import atexit
import code
import os
import readline
import rlcompleter

from diofant.interactive.session import (AutomaticSymbols,
                                         IntegerDivisionWrapper)


__all__ = ()


parser = argparse.ArgumentParser(description=__doc__,
                                 prog='python -m diofant')
parser.add_argument('--no-wrap-division',
                    help="Don't wrap integer divisions with Fraction",
                    action='store_true')
parser.add_argument('-a', '--auto-symbols',
                    help="Automatically create missing Symbol's",
                    action='store_true')
parser.add_argument('--no-ipython', help="Don't use IPython",
                    action='store_true')


def main():
    args, ipython_args = parser.parse_known_args()

    lines = ['from diofant import *',
             'init_printing()',
             "a, b, c, d, t, x, y, z = symbols('a:d t x:z')",
             "k, m, n = symbols('k m n', integer=True)",
             "f, g, h = symbols('f g h', cls=Function)",
             'init_printing(pretty_print=True, use_unicode=True)']

    try:
        import IPython
        import traitlets
    except ImportError:
        args.no_ipython = True

    if not args.no_ipython:
        config = traitlets.config.loader.Config()
        shell = config.InteractiveShell
        ast_transformers = shell.ast_transformers
        if not args.no_wrap_division:
            ast_transformers.append(IntegerDivisionWrapper())
        shell.confirm_exit = False
        config.TerminalIPythonApp.display_banner = False
        config.TerminalInteractiveShell.autoformatter = None

        app = IPython.terminal.ipapp.TerminalIPythonApp.instance(config=config)
        app.initialize(ipython_args)
        shell = app.shell
        for l in lines:
            shell.run_cell(l, silent=True)
        if args.auto_symbols:
            shell.run_cell('from diofant.interactive.session import AutomaticSymbols')
            shell.run_cell('ip = get_ipython()')
            shell.run_cell('ip.ast_transformers.append(AutomaticSymbols(ip.user_ns))')
            shell.run_cell('del ip')
        app.start()
    else:
        ast_transformers = []
        ns = {}

        if not args.no_wrap_division:
            ast_transformers.append(IntegerDivisionWrapper())
        if args.auto_symbols:
            ast_transformers.append(AutomaticSymbols(ns))

        class DiofantConsole(code.InteractiveConsole):
            """An interactive console with readline support."""

            def __init__(self, ast_transformers=[], **kwargs):
                super().__init__(**kwargs)

                readline.set_completer(rlcompleter.Completer(ns).complete)
                readline.parse_and_bind('tab: complete')

                history = os.path.expanduser('~/.python_history')
                readline.read_history_file(history)
                atexit.register(readline.write_history_file, history)
                self.ast_transformers = ast_transformers

            def runsource(self, source, filename='<input>', symbol='single'):
                try:
                    tree = ast.parse(source)
                except SyntaxError:
                    return True

                for t in self.ast_transformers:
                    tree = t.visit(tree)
                ast.fix_missing_locations(tree)

                source = ast.unparse(tree)
                source = source.split('\n')
                source = ';'.join(source)
                return super().runsource(source, filename=filename, symbol=symbol)

        c = DiofantConsole(ast_transformers=ast_transformers, locals=ns)

        for l in lines:
            c.push(l)
        c.interact('', '')


if __name__ == '__main__':  # pragma: no branch
    main()
