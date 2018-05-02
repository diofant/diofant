"""
Python shell for Diofant.

This is just a normal Python shell (IPython shell if you have the
IPython package installed), that adds default imports and run
some initialization code.
"""

import argparse
import atexit
import code
import os
import readline

import diofant


parser = argparse.ArgumentParser(description="Python shell for Diofant.")
parser.add_argument("--no-wrap-division",
                    help="Don't wrap integer divisions",
                    action="store_true")
parser.add_argument("-a", "--auto-symbols",
                    help="Automatically create missing symbols",
                    action="store_true")
parser.add_argument("--no-ipython",
                    help="Don't use IPython",
                    action="store_true")


def main():
    args, ipython_args = parser.parse_known_args()

    lines = ["from diofant import *",
             "init_printing()",
             "x, y, z, t = symbols('x y z t')",
             "k, m, n = symbols('k m n', integer=True)",
             "f, g, h = symbols('f g h', cls=Function)"]

    try:
        import IPython
        import traitlets
    except ImportError:
        args.no_ipython = False

    if not args.no_ipython:
        config = traitlets.config.loader.Config()
        if not args.no_wrap_division:
            config.InteractiveShell.ast_transformers.append(diofant.interactive.session.IntegerDivisionWrapper())
        if args.auto_symbols:
            config.InteractiveShell.ast_transformers.append(diofant.interactive.session.AutomaticSymbols())
        config.InteractiveShell.confirm_exit = False
        config.TerminalIPythonApp.display_banner = False

        app = IPython.terminal.ipapp.TerminalIPythonApp.instance(config=config)
        app.initialize(ipython_args)
        for l in lines:
            app.shell.run_cell(l, silent=True)
        app.start()
    else:
        class DiofantConsole(code.InteractiveConsole):
            """An interactive console with readline support. """

            def __init__(self):
                super().__init__()

                readline.parse_and_bind('tab: complete')

                history = os.path.expanduser('~/.python_history')
                readline.read_history_file(history)
                atexit.register(readline.write_history_file, history)

        c = DiofantConsole()
        for l in lines:
            c.push(l)
        c.interact("")


if __name__ == "__main__":
    main()
