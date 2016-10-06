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
    def __init__(self):
        super(AutomaticSymbols, self).__init__()
        self.names = []

    def visit_Module(self, node):
        import IPython

        app = IPython.get_ipython()
        ignored_names = list(app.user_ns.keys()) + dir(builtins)

        for s in node.body:
            self.visit(s)

        for v in self.names:
            if v in ignored_names:
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
