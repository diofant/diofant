import argparse

import IPython
import traitlets

import diofant


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-wrap-division",
                        help="Don't wrap integer divisions",
                        action="store_true")
    parser.add_argument("--auto-symbols",
                        help="Automatically create missing symbols",
                        action="store_true")

    args, ipython_args = parser.parse_known_args()

    config = traitlets.config.loader.Config()
    lines = ["from diofant import *",
             "init_printing()",
             "x, y, z, t = symbols('x y z t')",
             "k, m, n = symbols('k m n', integer=True)",
             "f, g, h = symbols('f g h', cls=Function)"]
    if not args.no_wrap_division:
        config.InteractiveShell.ast_transformers.append(diofant.interactive.session.IntegerDivisionWrapper())
    if args.auto_symbols:
        config.InteractiveShell.ast_transformers.append(diofant.interactive.session.AutomaticSymbols())
    config.InteractiveShellApp.exec_lines = lines

    IPython.start_ipython(argv=ipython_args, config=config)


if __name__ == "__main__":
    main()
