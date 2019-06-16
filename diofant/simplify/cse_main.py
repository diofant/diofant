""" Tools for doing common subexpression elimination."""

from ..core import Add, Basic, Mul, Pow, Symbol, Tuple, factor_terms, sympify
from ..core.compatibility import iterable
from ..core.function import _coeff_isneg
from ..utilities.iterables import (filter_symbols, numbered_symbols, ordered,
                                   sift, topological_sort)
from . import cse_opts


# (preprocessor, postprocessor) pairs which are commonly useful. They should
# each take a diofant expression and return a possibly transformed expression.
# When used in the function ``cse()``, the target expressions will be transformed
# by each of the preprocessor functions in order. After the common
# subexpressions are eliminated, each resulting expression will have the
# postprocessor functions transform them in *reverse* order in order to undo the
# transformation if necessary. This allows the algorithm to operate on
# a representation of the expressions that allows for more optimization
# opportunities.
# ``None`` can be used to specify no transformation for either the preprocessor or
# postprocessor.


basic_optimizations = [(cse_opts.sub_pre, cse_opts.sub_post),
                       (factor_terms, None)]

# sometimes we want the output in a different format; non-trivial
# transformations can be put here for users
# ===============================================================


def reps_toposort(r):
    """Sort replacements `r` so (k1, v1) appears before (k2, v2)
    if k2 is in v1's free symbols. This orders items in the
    way that cse returns its results (hence, in order to use the
    replacements in a substitution option it would make sense
    to reverse the order).

    Examples
    ========

    >>> for l, r in reps_toposort([(x, y + 1), (y, 2)]):
    ...     print(Eq(l, r))
    ...
    Eq(y, 2)
    Eq(x, y + 1)

    """
    r = sympify(r)
    E = []
    for c1, (k1, v1) in enumerate(r):
        for c2, (k2, v2) in enumerate(r):
            if k1 in v2.free_symbols:
                E.append((c1, c2))
    return [r[i] for i in topological_sort((range(len(r)), E))]


def cse_separate(r, e):
    """Move expressions that are in the form (symbol, expr) out of the
    expressions and sort them into the replacements using the reps_toposort.

    Examples
    ========

    >>> x0, x1 = symbols('x:2')
    >>> eq = (x + 1 + exp((x + 1)/(y + 1)) + cos(y + 1))
    >>> cse([eq, Eq(x, z + 1), z - 2], postprocess=cse_separate) in [
    ... [[(x0, y + 1), (x, z + 1), (x1, x + 1)],
    ...  [x1 + exp(x1/x0) + cos(x0), z - 2]],
    ... [[(x1, y + 1), (x, z + 1), (x0, x + 1)],
    ...  [x0 + exp(x0/x1) + cos(x1), z - 2]]]
    ...
    True

    """
    d = sift(e, lambda w: w.is_Equality and w.lhs.is_Symbol)
    r = r + [w.args for w in d[True]]
    e = d[False]
    return [reps_toposort(r), e]

# ====end of cse postprocess idioms===========================


def preprocess_for_cse(expr, optimizations):
    """ Preprocess an expression to optimize for common subexpression
    elimination.

    Parameters
    ----------
    expr : diofant expression
        The target expression to optimize.
    optimizations : list of (callable, callable) pairs
        The (preprocessor, postprocessor) pairs.

    Returns
    -------
    expr : diofant expression
        The transformed expression.

    """
    for pre, post in optimizations:
        if pre is not None:
            expr = pre(expr)
    return expr


def postprocess_for_cse(expr, optimizations):
    """ Postprocess an expression after common subexpression elimination to
    return the expression to canonical diofant form.

    Parameters
    ----------
    expr : diofant expression
        The target expression to transform.
    optimizations : list of (callable, callable) pairs, optional
        The (preprocessor, postprocessor) pairs.  The postprocessors will be
        applied in reversed order to undo the effects of the preprocessors
        correctly.

    Returns
    -------
    expr : diofant expression
        The transformed expression.

    """
    for pre, post in reversed(optimizations):
        if post is not None:
            expr = post(expr)
    return expr


def opt_cse(exprs, order='canonical'):
    """Find optimization opportunities in Adds, Muls, Pows and negative
    coefficient Muls

    Parameters
    ==========

    exprs : list of diofant expressions
        The expressions to optimize.
    order : string, 'none' or 'canonical'
        The order by which Mul and Add arguments are processed. For large
        expressions where speed is a concern, use the setting order='none'.

    Returns
    =======

    opt_subs : dictionary of expression substitutions
        The expression substitutions which can be useful to optimize CSE.

    Examples
    ========

    >>> opt_subs = opt_cse([x**-2])
    >>> opt_subs
    {x**(-2): 1/(x**2)}

    """
    opt_subs = {}

    adds = set()
    muls = set()

    seen_subexp = set()

    def _find_opts(expr):

        if expr.is_Atom or expr.is_Order:
            return

        if iterable(expr):
            list(map(_find_opts, expr))
            return

        if expr in seen_subexp:
            return expr
        seen_subexp.add(expr)

        list(map(_find_opts, expr.args))

        if _coeff_isneg(expr):
            neg_expr = -expr
            if not neg_expr.is_Atom:
                opt_subs[expr] = Mul(-1, neg_expr, evaluate=False)
                seen_subexp.add(neg_expr)
                expr = neg_expr

        if expr.is_Mul:
            muls.add(expr)

        elif expr.is_Add:
            adds.add(expr)

        elif expr.is_Pow:
            if _coeff_isneg(expr.exp):
                opt_subs[expr] = Pow(Pow(expr.base, -expr.exp), -1, evaluate=False)

    for e in exprs:
        assert isinstance(e, Basic)
        _find_opts(e)

    # Process Adds and commutative Muls

    def _match_common_args(Func, funcs):
        if order != 'none':
            funcs = list(ordered(funcs))
        else:
            funcs = sorted(funcs, key=lambda x: len(x.args))

        func_args = [set(e.args) for e in funcs]
        for i in range(len(func_args)):
            for j in range(i + 1, len(func_args)):
                com_args = func_args[i].intersection(func_args[j])
                if len(com_args) > 1:
                    com_func = Func(*com_args)

                    # for all sets, replace the common symbols by the function
                    # over them, to allow recursive matches

                    diff_i = func_args[i].difference(com_args)
                    func_args[i] = diff_i | {com_func}
                    if diff_i:
                        opt_subs[funcs[i]] = Func(Func(*diff_i), com_func,
                                                  evaluate=False)

                    diff_j = func_args[j].difference(com_args)
                    func_args[j] = diff_j | {com_func}
                    opt_subs[funcs[j]] = Func(Func(*diff_j), com_func,
                                              evaluate=False)

                    for k in range(j + 1, len(func_args)):
                        if not com_args.difference(func_args[k]):
                            diff_k = func_args[k].difference(com_args)
                            func_args[k] = diff_k | {com_func}
                            opt_subs[funcs[k]] = Func(Func(*diff_k), com_func,
                                                      evaluate=False)

    # split muls into commutative
    comutative_muls = set()
    for m in muls:
        c, nc = m.args_cnc(cset=True)
        if c:
            c_mul = Mul(*c)
            if nc:
                opt_subs[m] = Mul(c_mul, Mul(*nc), evaluate=False)
            if len(c) > 1:
                comutative_muls.add(c_mul)

    _match_common_args(Add, adds)
    _match_common_args(Mul, comutative_muls)

    return opt_subs


def tree_cse(exprs, symbols, opt_subs={}, order='canonical', ignore=()):
    """Perform raw CSE on expression tree, taking opt_subs into account.

    Parameters
    ==========

    exprs : list of diofant expressions
        The expressions to reduce.
    symbols : infinite iterator yielding unique Symbols
        The symbols used to label the common subexpressions which are pulled
        out.
    opt_subs : dictionary of expression substitutions
        The expressions to be substituted before any CSE action is performed.
    order : string, 'none' or 'canonical'
        The order by which Mul and Add arguments are processed. For large
        expressions where speed is a concern, use the setting order='none'.
    ignore : iterable of Symbol's
        Substitutions containing these symbols will be ignored.

    """
    # Find repeated sub-expressions

    to_eliminate = set()

    seen_subexp = set()

    def _find_repeated(expr):
        if expr.is_Atom or expr.is_Order:
            return

        if iterable(expr):
            args = expr

        else:
            if expr in seen_subexp:
                if not any(ign in expr.free_symbols for ign in ignore):
                    to_eliminate.add(expr)
                    return

            seen_subexp.add(expr)

            if expr in opt_subs:
                expr = opt_subs[expr]

            args = expr.args

        list(map(_find_repeated, args))

    for e in exprs:
        assert isinstance(e, Basic)
        _find_repeated(e)

    # Rebuild tree

    replacements = []

    subs = {}

    def _rebuild(expr):

        if not expr.args:
            return expr

        if iterable(expr):
            new_args = [_rebuild(arg) for arg in expr]
            return expr.func(*new_args)

        if expr in subs:
            return subs[expr]

        orig_expr = expr
        if expr in opt_subs:
            expr = opt_subs[expr]

        # If enabled, parse Muls and Adds arguments by order to ensure
        # replacement order independent from hashes
        if order != 'none':
            if expr.is_Mul:
                c, nc = expr.args_cnc()
                args = list(ordered(c)) + nc
            elif expr.is_Add:
                args = list(ordered(expr.args))
            else:
                args = expr.args
        else:
            args = expr.args

        new_args = list(map(_rebuild, args))
        if new_args != args:
            new_expr = expr.func(*new_args)
        else:
            new_expr = expr

        if orig_expr in to_eliminate:
            try:
                sym = next(symbols)
            except StopIteration:
                raise ValueError("Symbols iterator ran out of symbols.")
            if not orig_expr.is_commutative and not orig_expr.is_Relational:
                sym = Symbol(sym.name, commutative=False)
            subs[orig_expr] = sym
            replacements.append((sym, new_expr))
            return sym

        else:
            return new_expr

    reduced_exprs = []
    for e in exprs:
        assert isinstance(e, Basic)
        reduced_e = _rebuild(e)
        reduced_exprs.append(reduced_e)

    return replacements, reduced_exprs


def cse(exprs, symbols=None, optimizations=None, postprocess=None,
        order='canonical', ignore=()):
    """ Perform common subexpression elimination on an expression.

    Parameters
    ==========

    exprs : list of diofant expressions, or a single diofant expression
        The expressions to reduce.
    symbols : infinite iterator yielding unique Symbols
        The symbols used to label the common subexpressions which are pulled
        out. The ``numbered_symbols`` generator is useful. The default is a
        stream of symbols of the form "x0", "x1", etc. This must be an
        infinite iterator.
    optimizations : list of (callable, callable) pairs
        The (preprocessor, postprocessor) pairs of external optimization
        functions. Optionally 'basic' can be passed for a set of predefined
        basic optimizations. Such 'basic' optimizations were used by default
        in old implementation, however they can be really slow on larger
        expressions. Now, no pre or post optimizations are made by default.
    postprocess : a function which accepts the two return values of cse and
        returns the desired form of output from cse, e.g. if you want the
        replacements reversed the function might be the following lambda:
        lambda r, e: return reversed(r), e
    order : string, 'none' or 'canonical'
        The order by which Mul and Add arguments are processed. If set to
        'canonical', arguments will be canonically ordered. If set to 'none',
        ordering will be faster but dependent on expressions hashes, thus
        machine dependent and variable. For large expressions where speed is a
        concern, use the setting order='none'.
    ignore : iterable of Symbol's
        Substitutions containing these symbols will be ignored.

    Returns
    =======

    replacements : list of (Symbol, expression) pairs
        All of the common subexpressions that were replaced. Subexpressions
        earlier in this list might show up in subexpressions later in this
        list.
    reduced_exprs : list of diofant expressions
        The reduced expressions with all of the replacements above.

    Examples
    ========

    >>> from diofant.abc import w
    >>> cse(((w + x + y + z)*(w + y + z))/(w + x)**3)
    ([(x0, y + z), (x1, w + x)], [(w + x0)*(x0 + x1)/x1**3])

    Note that currently, y + z will not get substituted if -y - z is used.

     >>> cse(((w + x + y + z)*(w - y - z))/(w + x)**3)
     ([(x0, w + x)], [(w - y - z)*(x0 + y + z)/x0**3])

    List of expressions with recursive substitutions:

    >>> m = SparseMatrix([x + y, x + y + z])
    >>> cse([(x+y)**2, x + y + z, y + z, x + z + y, m])
    ([(x0, x + y), (x1, x0 + z)], [x0**2, x1, y + z, x1, Matrix([
    [x0],
    [x1]])])

    Note: the type and mutability of input matrices is retained.

    >>> isinstance(_[1][-1], SparseMatrix)
    True

    The user may disallow substitutions containing certain symbols:
    >>> cse([y**2*(x + 1), 3*y**2*(x + 1)], ignore=(y,))
    ([(x0, x + 1)], [x0*y**2, 3*x0*y**2])

    """
    from ..matrices import (MatrixBase, Matrix, ImmutableMatrix,
                            SparseMatrix, ImmutableSparseMatrix)

    # Handle the case if just one expression was passed.
    if isinstance(exprs, (Basic, MatrixBase)):
        exprs = [exprs]

    copy = exprs
    temp = []
    for e in exprs:
        if isinstance(e, (Matrix, ImmutableMatrix)):
            temp.append(Tuple(*e._mat))
        elif isinstance(e, (SparseMatrix, ImmutableSparseMatrix)):
            temp.append(Tuple(*e._smat.items()))
        else:
            temp.append(e)
    exprs = temp
    del temp

    if optimizations is None:
        optimizations = []
    elif optimizations == 'basic':
        optimizations = basic_optimizations

    # Preprocess the expressions to give us better optimization opportunities.
    reduced_exprs = [preprocess_for_cse(e, optimizations) for e in exprs]

    excluded_symbols = set().union(*[expr.atoms(Symbol)
                                     for expr in reduced_exprs])

    if symbols is None:
        symbols = numbered_symbols()
    else:
        # In case we get passed an iterable with an __iter__ method instead of
        # an actual iterator.
        symbols = iter(symbols)

    symbols = filter_symbols(symbols, excluded_symbols)

    # Find other optimization opportunities.
    opt_subs = opt_cse(reduced_exprs, order)

    # Main CSE algorithm.
    replacements, reduced_exprs = tree_cse(reduced_exprs, symbols, opt_subs,
                                           order, ignore)

    # Postprocess the expressions to return the expressions to canonical form.
    exprs = copy
    for i, (sym, subtree) in enumerate(replacements):
        subtree = postprocess_for_cse(subtree, optimizations)
        replacements[i] = (sym, subtree)
    reduced_exprs = [postprocess_for_cse(e, optimizations)
                     for e in reduced_exprs]

    # Get the matrices back
    for i, e in enumerate(exprs):
        if isinstance(e, (Matrix, ImmutableMatrix)):
            reduced_exprs[i] = Matrix(e.rows, e.cols, reduced_exprs[i])
            if isinstance(e, ImmutableMatrix):
                reduced_exprs[i] = reduced_exprs[i].as_immutable()
        elif isinstance(e, (SparseMatrix, ImmutableSparseMatrix)):
            m = SparseMatrix(e.rows, e.cols, {})
            for k, v in reduced_exprs[i]:
                m[k] = v
            if isinstance(e, ImmutableSparseMatrix):
                m = m.as_immutable()
            reduced_exprs[i] = m

    if postprocess is None:
        return replacements, reduced_exprs

    return postprocess(replacements, reduced_exprs)
