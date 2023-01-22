from ..core import Basic, Expr


__all__ = 'dotprint',

default_styles = [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
                  (Expr,  {'color': 'black'})]


def styleof(expr, styles=default_styles):
    """Merge style dictionaries in order.

    >>> styles = [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
    ...           (Expr, {'color': 'black'})]

    >>> styleof(Basic(1), styles)
    {'color': 'blue', 'shape': 'ellipse'}

    >>> styleof(x + 1, styles)  # this is an Expr
    {'color': 'black', 'shape': 'ellipse'}

    """
    style = {}
    for typ, sty in styles:
        if isinstance(expr, typ):
            style.update(sty)
    return style


def attrprint(d, delimiter=', '):
    """Print a dictionary of attributes.

    >>> print(attrprint({'color': 'blue', 'shape': 'ellipse'}))
    "color"="blue", "shape"="ellipse"

    """
    return delimiter.join(f'"{k}"="{v}"' for k, v in sorted(d.items()))


def dotnode(expr, styles=default_styles, labelfunc=str, pos=(), repeat=True):
    """String defining a node.

    >>> print(dotnode(x))
    "Symbol('x')_()" ["color"="black", "label"="x", "shape"="ellipse"];

    """
    style = styleof(expr, styles)

    if isinstance(expr, Basic) and not expr.is_Atom:
        label = str(expr.__class__.__name__)
    else:
        label = labelfunc(expr)
    style['label'] = label
    expr_str = repr(expr)
    if repeat:
        expr_str += f'_{pos!s}'
    return f'"{expr_str}" [{attrprint(style)}];'


def dotedges(expr, atom=lambda x: not isinstance(x, Basic), pos=(), repeat=True):
    """List of strings for all expr->expr.arg pairs.

    >>> for e in dotedges(x+2):
    ...     print(e)
    "Add(Symbol('x'), Integer(2))_()" -> "Integer(2)_(0,)";
    "Add(Symbol('x'), Integer(2))_()" -> "Symbol('x')_(1,)";

    See Also
    ========

    dotprint

    """
    if atom(expr):
        return []
    expr_str = repr(expr)
    arg_strs = [repr(arg) for arg in expr.args]
    if repeat:
        expr_str += f'_{pos!s}'
        arg_strs = [arg_str + f'_{pos + (i,)!s}' for i, arg_str in enumerate(arg_strs)]
    return [f'"{expr_str}" -> "{arg_str}";' for arg_str in arg_strs]


template = \
    """digraph{

# Graph style
%(graphstyle)s

#########
# Nodes #
#########

%(nodes)s

#########
# Edges #
#########

%(edges)s
}"""

graphstyle = {'rankdir': 'TD', 'ordering': 'out', 'bgcolor': 'transparent'}


def dotprint(expr, styles=default_styles,
             atom=lambda x: not isinstance(x, Basic),
             maxdepth=None, repeat=True, labelfunc=str, **kwargs):
    """
    DOT description of a Diofant expression tree

    Options are

    ``styles``: Styles for different classes.  The default is::

        [(Basic, {'color': 'blue', 'shape': 'ellipse'}),
        (Expr, {'color': 'black'})]``

    ``atom``: Function used to determine if an arg is an atom.  The default is
          ``lambda x: not isinstance(x, Basic)``.  Another good choice is
          ``lambda x: not x.args``.

    ``maxdepth``: The maximum depth.  The default is None, meaning no limit.

    ``repeat``: Whether to different nodes for separate common subexpressions.
          The default is True.  For example, for ``x + x*y`` with
          ``repeat=True``, it will have two nodes for ``x`` and with
          ``repeat=False``, it will have one (warning: even if it appears
          twice in the same object, like Pow(x, x), it will still only appear
          only once.  Hence, with repeat=False, the number of arrows out of an
          object might not equal the number of args it has).

    ``labelfunc``: How to label leaf nodes.  The default is ``str``.  Another
          good option is ``repr``. For example with ``str``, the leaf nodes
          of ``x + 1`` are labeled, ``x`` and ``1``.  With ``repr``, they
          are labeled ``Symbol('x')`` and ``Integer(1)``.

    Additional keyword arguments are included as styles for the graph.

    Examples
    ========

    >>> print(dotprint(x + 2))
    digraph{
    <BLANKLINE>
    # Graph style
    "bgcolor"="transparent"
    "ordering"="out"
    "rankdir"="TD"
    <BLANKLINE>
    #########
    # Nodes #
    #########
    <BLANKLINE>
    "Add(Symbol('x'), Integer(2))_()" ["color"="black", "label"="Add", "shape"="ellipse"];
    "Integer(2)_(0,)" ["color"="black", "label"="2", "shape"="ellipse"];
    "Symbol('x')_(1,)" ["color"="black", "label"="x", "shape"="ellipse"];
    <BLANKLINE>
    #########
    # Edges #
    #########
    <BLANKLINE>
    "Add(Symbol('x'), Integer(2))_()" -> "Integer(2)_(0,)";
    "Add(Symbol('x'), Integer(2))_()" -> "Symbol('x')_(1,)";
    }

    """
    # repeat works by adding a signature tuple to the end of each node for its
    # position in the graph. For example, for expr = Add(x, Pow(x, 2)), the x in the
    # Pow will have the tuple (1, 0), meaning it is expr.args[1].args[0].
    graphstyle.update(kwargs)

    nodes = []
    edges = []

    def traverse(e, depth, pos=()):
        nodes.append(dotnode(e, styles, labelfunc=labelfunc, pos=pos, repeat=repeat))
        if maxdepth and depth >= maxdepth:
            return
        edges.extend(dotedges(e, atom=atom, pos=pos, repeat=repeat))
        for i, arg in enumerate(e.args):
            if not atom(arg):
                traverse(arg, depth+1, pos + (i,))

    traverse(expr, 0)

    return template % {'graphstyle': attrprint(graphstyle, delimiter='\n'),
                       'nodes': '\n'.join(nodes),
                       'edges': '\n'.join(edges)}
