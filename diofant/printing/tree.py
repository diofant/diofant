def pprint_nodes(subtrees):
    r"""
    Prettyprints systems of nodes.

    Examples
    ========

    >>> print(pprint_nodes(["a", "b1\nb2", "c"]))
    +-a
    +-b1
    | b2
    +-c
    """
    def indent(s, type=1):
        x = s.split("\n")
        r = "+-%s\n" % x[0]
        for a in x[1:]:
            if a == "":
                continue
            if type == 1:
                r += "| %s\n" % a
            else:
                r += "  %s\n" % a
        return r
    if len(subtrees) == 0:
        return ""
    f = ""
    for a in subtrees[:-1]:
        f += indent(a)
    f += indent(subtrees[-1], 2)
    return f


def print_node(node):
    """
    Returns information about the "node".

    This includes class name, string representation and assumptions.
    """
    s = "%s: %s\n" % (node.__class__.__name__, str(node))
    d = {k: v for k, v in node._assumptions.items()
         if v is not None}
    for a in sorted(d):
        s += "%s: %s\n" % (a, d[a])
    return s


def tree(node):
    """
    Returns a tree representation of "node" as a string.

    It uses print_node() together with pprint_nodes() on node.args recursively.

    See also: print_tree()
    """
    subtrees = []
    for arg in node.args:
        subtrees.append(tree(arg))
    s = print_node(node) + pprint_nodes(subtrees)
    return s


def print_tree(node):
    """
    Prints a tree representation of "node".

    Examples
    ========

    >>> x = Symbol('x', odd=True)
    >>> y = Symbol('y', even=True)
    >>> print_tree(y**x)
    Pow: y**x
    commutative: True
    +-Symbol: y
    | algebraic: True
    | commutative: True
    | complex: True
    | even: True
    | extended_real: True
    | finite: True
    | hermitian: True
    | infinite: False
    | integer: True
    | irrational: False
    | noninteger: False
    | odd: False
    | rational: True
    | real: True
    | transcendental: False
    +-Symbol: x
      algebraic: True
      commutative: True
      complex: True
      even: False
      extended_real: True
      finite: True
      hermitian: True
      imaginary: False
      infinite: False
      integer: True
      irrational: False
      noninteger: False
      nonzero: True
      odd: True
      rational: True
      real: True
      transcendental: False
      zero: False

    See also: tree()
    """
    print(tree(node))
