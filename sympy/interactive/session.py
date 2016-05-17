"""Tools for setting up interactive sessions. """

import ast


class IntegerWrapper(ast.NodeTransformer):
    """Wraps all integers in a call to Integer()"""
    def visit_Num(self, node):
        if isinstance(node.n, int):
            return ast.Call(func=ast.Name(id='Integer', ctx=ast.Load()),
                            args=[node], keywords=[],
                            starargs=None, kwargs=None)
        return node

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id != "Integer":
            node = self.generic_visit(node)
        return node


def enable_automatic_symbols(app):
    """Allow IPython to automatially create symbols. """
    # XXX: This should perhaps use ast, like IntegerWrapper above.
    # This would avoid re-executing the code, which can lead to subtle
    # issues.  For example:
    #
    # In [1]: a = 1
    #
    # In [2]: for i in range(3):
    #    ...:     a += 1
    #    ...:
    #
    # In [3]: a
    # Out[3]: 4
    #
    # In [4]: a = 1
    #
    # In [5]: for i in range(3):
    #    ...:     a += 1
    #    ...:     print(b)
    #    ...:
    # b
    # b
    # b
    #
    # In [6]: a
    # Out[6]: 5
    #
    # Note how the for loop is executed again because `b` was not defined, but `a`
    # was already incremented once, so the result is that it is incremented
    # multiple times.

    import re
    re_nameerror = re.compile(
        "name '(?P<symbol>[A-Za-z_][A-Za-z0-9_]*)' is not defined")

    def _handler(self, etype, value, tb, tb_offset=None):
        """Handle :exc:`NameError` exception and allow injection of missing symbols. """
        match = re_nameerror.match(str(value))

        if match is not None:
            # XXX: Make sure Symbol is in scope. Otherwise you'll get infinite recursion.
            self.run_cell("%(symbol)s = Symbol('%(symbol)s')" %
                          {'symbol': match.group("symbol")},
                           store_history=False)

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

    app.set_custom_exc((NameError,), _handler)


def init_ipython_session(argv=[], auto_symbols=False,
                         auto_int_to_Integer=False):
    """Construct new IPython session. """
    import IPython
    import readline

    app = IPython.terminal.ipapp.TerminalIPythonApp()

    # don't draw IPython banner during initialization:
    app.display_banner = False
    app.initialize(argv)

    if auto_symbols:
        enable_automatic_symbols(app.shell)
    if auto_int_to_Integer:
        app.shell.ast_transformers.append(IntegerWrapper())

    return app.shell
