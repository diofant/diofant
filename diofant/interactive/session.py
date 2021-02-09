import ast
import builtins


class IntegerDivisionWrapper(ast.NodeTransformer):
    """Wrap all int divisions in a call to Rational."""

    def visit_BinOp(self, node):
        def is_integer(x):
            if isinstance(x, ast.Num) and isinstance(x.n, int):
                return True
            elif isinstance(x, ast.UnaryOp) and isinstance(x.op, (ast.USub,
                                                                  ast.UAdd)):
                return is_integer(x.operand)
            elif isinstance(x, ast.BinOp) and isinstance(x.op, (ast.Add,
                                                                ast.Sub,
                                                                ast.Mult,
                                                                ast.Pow)):
                return is_integer(x.left) and is_integer(x.right)
            else:
                return False

        if (isinstance(node.op, ast.Div) and
                all(is_integer(_) for _ in [node.left, node.right])):
            return ast.Call(func=ast.Name(id='Rational', ctx=ast.Load()),
                            args=[node.left, node.right], keywords=[],
                            starargs=None, kwargs=None)
        return self.generic_visit(node)


class AutomaticSymbols(ast.NodeTransformer):
    """Add missing Symbol definitions automatically."""

    def __init__(self):
        super().__init__()
        self.names = []

    def visit_Module(self, node):
        import IPython

        app = IPython.get_ipython()
        ignored_names = list(app.user_ns) + dir(builtins)

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


class FloatRationalizer(ast.NodeTransformer):
    """Wraps all floats in a call to Rational."""

    def visit_Constant(self, node):
        if isinstance(node.n, float):
            return ast.Call(func=ast.Name(id='Rational', ctx=ast.Load()),
                            args=[ast.Str(s=repr(node.n))], keywords=[],
                            starargs=None, kwargs=None)
        return node

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id == 'Float':
            return node
        return self.generic_visit(node)
