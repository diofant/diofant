import argparse
import atexit
import code
import os
import readline

import diofant


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-wrap-division",
                        help="Don't wrap integer divisions",
                        action="store_true")
    parser.add_argument("--auto-symbols",
                        help="Automatically create missing symbols",
                        action="store_true")
    parser.add_argument("--no-ipython",
                        help="Don't use IPython",
                        action="store_true")

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
        config.InteractiveShellApp.exec_lines = lines

        IPython.start_ipython(argv=ipython_args, config=config)
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
        banner_python = "\n".join(">>> " + l for l in lines)
        for l in lines:
            c.push(l)
            readline.add_history(l)
        c.interact(banner_python)


if __name__ == "__main__":
    main()
