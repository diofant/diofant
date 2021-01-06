from diofant import Basic, Expr, Integer, MatrixSymbol, Symbol
from diofant.abc import x
from diofant.printing.dot import (attrprint, dotedges, dotnode, dotprint,
                                  styleof)


__all__ = ()


def test_styleof():
    styles = [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
              (Expr,  {'color': 'black'})]
    assert styleof(Basic(1), styles) == {'color': 'blue', 'shape': 'ellipse'}

    assert styleof(x + 1, styles) == {'color': 'black', 'shape': 'ellipse'}


def test_attrprint():
    assert attrprint({'color': 'blue', 'shape': 'ellipse'}) == \
        '"color"="blue", "shape"="ellipse"'


def test_dotnode():
    assert dotnode(x, repeat=False) ==\
        '"Symbol(\'x\')" ["color"="black", "label"="x", "shape"="ellipse"];'
    assert dotnode(x + 2, repeat=False) == \
        '"Add(Symbol(\'x\'), Integer(2))" ["color"="black", "label"="Add", "shape"="ellipse"];'
    assert dotnode(x + x**2, repeat=False) == \
        '"Add(Pow(Symbol(\'x\'), Integer(2)), Symbol(\'x\'))" ["color"="black", "label"="Add", "shape"="ellipse"];'
    assert dotnode(x + x**2, repeat=True) == \
        '"Add(Pow(Symbol(\'x\'), Integer(2)), Symbol(\'x\'))_()" ["color"="black", "label"="Add", "shape"="ellipse"];'


def test_dotedges():
    assert sorted(dotedges(x + 2, repeat=False)) == [
        '"Add(Symbol(\'x\'), Integer(2))" -> "Integer(2)";',
        '"Add(Symbol(\'x\'), Integer(2))" -> "Symbol(\'x\')";'
    ]
    assert sorted(dotedges(x + 2, repeat=True)) == [
        '"Add(Symbol(\'x\'), Integer(2))_()" -> "Integer(2)_(0,)";',
        '"Add(Symbol(\'x\'), Integer(2))_()" -> "Symbol(\'x\')_(1,)";'
    ]
    assert dotedges(1) == []


def test_dotprint():
    text = dotprint(x+2, repeat=False)
    assert all(e in text for e in dotedges(x+2, repeat=False))
    assert all(n in text for n in [dotnode(expr, repeat=False) for expr in (x, Integer(2), x+2)])
    assert 'digraph' in text
    text = dotprint(x+x**2, repeat=False)
    assert all(e in text for e in dotedges(x+x**2, repeat=False))
    assert all(n in text for n in [dotnode(expr, repeat=False) for expr in (x, Integer(2), x**2)])
    assert 'digraph' in text
    text = dotprint(x+x**2, repeat=True)
    assert all(e in text for e in dotedges(x+x**2, repeat=True))
    assert all(n in text for n in [dotnode(expr, pos=()) for expr in [x + x**2]])
    text = dotprint(x**x, repeat=True)
    assert all(e in text for e in dotedges(x**x, repeat=True))
    assert all(n in text for n in [dotnode(x, pos=(0,)), dotnode(x, pos=(1,))])
    assert 'digraph' in text

    assert dotprint(x**x**x) != dotprint(x**x**x, maxdepth=1)
    assert (dotprint(x**x**x, maxdepth=-1) ==
            'digraph{\n\n# Graph style\n"bgcolor"="transparent"\n"ordering"'
            '="out"\n"rankdir"="TD"\n\n#########\n# Nodes #\n#########\n\n"'
            'Pow(Symbol(\'x\'), Pow(Symbol(\'x\'), Symbol(\'x\')))_()" '
            '["color"="black", "label"="Pow", "shape"="ellipse"];\n\n'
            '#########\n# Edges #\n#########\n\n\n}')


def test_dotprint_depth():
    text = dotprint(3*x+2, depth=1)
    assert dotnode(3*x+2) in text
    assert dotnode(x) not in text


def test_Matrix_and_non_basics():
    n = Symbol('n')
    assert dotprint(MatrixSymbol('X', n, n))


def test_labelfunc():
    text = dotprint(x + 2, labelfunc=repr)
    assert "Symbol('x')" in text
    assert 'Integer(2)' in text
