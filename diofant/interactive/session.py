"""Tools for setting up interactive sessions. """

import ast
import builtins


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


class AutomaticSymbols(ast.NodeTransformer):
    """Add missing Symbol definitions automatically."""
    def __init__(self, app):
        super(AutomaticSymbols, self).__init__()
        self.app = app
        self.names = []

    def visit_Module(self, node):
        for s in node.body:
            self.visit(s)

        for v in self.names:
            if v in self.app.user_ns.keys() or v in builtins.__dir__():
                continue

            assign = ast.Assign(targets=[ast.Name(id=v, ctx=ast.Store())],
                                value=ast.Call(func=ast.Name(id='Symbol',
                                                             ctx=ast.Load()),
                                args=[ast.Str(s=v)], keywords=[],
                                starargs=None, kwargs=None))
            node.body.insert(0, assign)

        newnode = ast.Module(body=node.body)
        ast.copy_location(newnode, node)
        ast.fix_missing_locations(newnode)
        return newnode

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Load):
            self.names.append(node.id)
        return node


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
        app.shell.ast_transformers.append(AutomaticSymbols(app.shell))
    if auto_int_to_Integer:
        app.shell.ast_transformers.append(IntegerWrapper())

    return app.shell
