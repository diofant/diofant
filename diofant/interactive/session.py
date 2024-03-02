import ast
import builtins
import io
import tokenize
import unicodedata
import uuid


class IntegerDivisionWrapper(ast.NodeTransformer):
    """Wrap all int divisions in a call to :class:`~fractions.Fraction`."""

    def visit_BinOp(self, node):
        def is_integer(x):
            if isinstance(x, ast.Constant) and isinstance(x.value, int):
                return True
            if isinstance(x, ast.UnaryOp) and isinstance(x.op, (ast.USub,
                                                                ast.UAdd)):
                return is_integer(x.operand)
            if isinstance(x, ast.BinOp) and isinstance(x.op, (ast.Add,
                                                              ast.Sub,
                                                              ast.Mult,
                                                              ast.Pow)):
                return is_integer(x.left) and is_integer(x.right)
            return False

        if isinstance(node.op, ast.Div) and all(map(is_integer,
                                                    [node.left, node.right])):
            return ast.Call(ast.Name('Fraction', ast.Load()),
                            [node.left, node.right], [])
        return self.generic_visit(node)


class AutomaticSymbols(ast.NodeTransformer):
    """Add missing :class:`~diofant.core.symbol.Symbol` definitions."""

    def __init__(self, ns={}):
        """Initialize self."""
        super().__init__()
        self.names = []
        self.ns = ns

    def visit_Module(self, node):
        ignored_names = list(self.ns) + dir(builtins)

        for s in node.body:
            self.visit(s)

        for v in self.names:
            if v in ignored_names:
                continue

            assign = ast.Assign([ast.Name(v, ast.Store())],
                                ast.Call(ast.Name('Symbol', ast.Load()),
                                         [ast.Constant(v)], []))
            node.body.insert(0, assign)

        newnode = ast.Module(node.body, [])
        ast.copy_location(newnode, node)
        ast.fix_missing_locations(newnode)
        return newnode

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Load):
            self.names.append(node.id)
        return node


class FloatRationalizer(ast.NodeTransformer):
    """Wraps all floats in a call to :class:`~fractions.Fraction`."""

    def visit_Constant(self, node):
        if isinstance(node.value, float):
            return ast.Call(ast.Name('Fraction', ast.Load()),
                            [ast.Constant(repr(node.value))], [])
        return node

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id == 'Float':
            return node
        return self.generic_visit(node)


_NAMES_MAP = {}


def unicode_identifiers(lines):
    """Transform original code to allow any unicode identifiers."""
    new_lines = []
    for line in lines:
        result = []
        g = tokenize.tokenize(io.BytesIO(line.encode()).readline)
        for toknum, tokval, _, _, _ in g:
            if toknum == tokenize.NAME:
                if unicodedata.normalize('NFKC', tokval) != tokval:
                    if tokval not in _NAMES_MAP:
                        _NAMES_MAP[tokval] = f'_{uuid.uuid4().hex!s}'
                    tokval = _NAMES_MAP[tokval]
            result.append((toknum, tokval))
        new_lines.append(tokenize.untokenize(result).decode())
    return new_lines


class _WrapFloats(ast.NodeTransformer):
    def __init__(self, lines):
        super().__init__()
        self.lines = lines

    def visit_Constant(self, node):
        if isinstance(node.value, float):
            line = self.lines[node.lineno - 1]
            value = line[node.col_offset:node.end_col_offset]
            return ast.Call(ast.Name('Float', ast.Load()),
                            [ast.Constant(value)], [])
        return node

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id == 'Float':
            return node
        return self.generic_visit(node)


def wrap_float_literals(lines):
    """Wraps all float literals with :class:`~diofant.core.numbers.Float`."""
    source = '\n'.join(lines)
    tree = ast.parse(source)
    tree = _WrapFloats(lines).visit(tree)
    ast.fix_missing_locations(tree)
    source = ast.unparse(tree)
    return source.split('\n')
