from .basic import _aresame
from .cache import cacheit
from .compatibility import ordered
from .evaluate import global_evaluate
from .expr import Expr
from .sympify import sympify


class AssocOp(Expr):
    """ Associative operations, can separate noncommutative and
    commutative parts.

    (a op b) op c == a op (b op c) == a op b op c.

    Base class for Add and Mul.

    This is an abstract base class, concrete derived classes must define
    the attribute `identity`.

    """

    @cacheit
    def __new__(cls, *args, **options):
        from ..series import Order
        args = [sympify(a, strict=True) for a in args]

        if not options.pop('evaluate', global_evaluate[0]):
            return cls._from_args(args)
        else:
            args = [a for a in args if a is not cls.identity]

        if len(args) == 0:
            return cls.identity
        if len(args) == 1:
            return args[0]

        c_part, nc_part, order_symbols = cls.flatten(args)
        obj = cls._from_args(c_part + nc_part)

        if order_symbols is not None:
            return Order(obj, *order_symbols)
        return obj

    @classmethod
    def _from_args(cls, args):
        """Create new instance with already-processed args."""
        if len(args) == 0:
            return cls.identity
        elif len(args) == 1:
            return args[0]

        return super().__new__(cls, *args)

    def _new_rawargs(self, *args, **kwargs):
        """Create new instance of own class with args exactly as provided by
        caller but returning the self class identity if args is empty.

           This is handy when we want to optimize things, e.g.

               >>> e = Mul(3, x, y)
               >>> e.args
               (3, x, y)
               >>> Mul(*e.args[1:])
               x*y
               >>> e._new_rawargs(*e.args[1:])  # the same as above, but faster
               x*y

           Note: use this with caution. There is no checking of arguments at
           all. This is best used when you are rebuilding an Add or Mul after
           simply removing one or more terms. If modification which result,
           for example, in extra 1s being inserted (as when collecting an
           expression's numerators and denominators) they will not show up in
           the result but a Mul will be returned nonetheless:

               >>> m = (x*y)._new_rawargs(Integer(1), x); m
               x
               >>> m == x
               False
               >>> m.is_Mul
               True

           Another issue to be aware of is that the commutativity of the result
           is based on the commutativity of self. If you are rebuilding the
           terms that came from a commutative object then there will be no
           problem, but if self was non-commutative then what you are
           rebuilding may now be commutative.

           Although this routine tries to do as little as possible with the
           input, getting the commutativity right is important, so this level
           of safety is enforced: commutativity will always be recomputed if
           self is non-commutative and kwarg `reeval=False` has not been
           passed.

        """
        return self._from_args(args)

    @classmethod
    def flatten(cls, seq):
        """Return seq so that none of the elements are of type `cls`.

        This is the vanilla routine that will be used if a class derived
        from AssocOp does not define its own flatten routine.

        """
        # apply associativity, no commutativity property is used
        new_seq = []
        for o in seq:
            if o.__class__ is cls:  # classes must match exactly
                seq.extend(o.args)
            else:
                new_seq.append(o)
        return [], new_seq, None  # c_part, nc_part, order_symbols

    def _matches_commutative(self, expr, repl_dict={}):
        """
        Matches Add/Mul "pattern" to an expression "expr".

        repl_dict ... a dictionary of (wild: expression) pairs, that get
                      returned with the results

        This function is the main workhorse for Add/Mul.

        For instance:

        >>> a = Wild("a")
        >>> b = Wild("b")
        >>> c = Wild("c")
        >>> (a + sin(b)*c)._matches_commutative(x + sin(y)*z)
        {a_: x, b_: y, c_: z}

        In the example above, "a+sin(b)*c" is the pattern, and "x+sin(y)*z" is
        the expression.

        The repl_dict contains parts that were already matched. For example
        here:

        >>> (x + sin(b)*c)._matches_commutative(x + sin(y)*z, repl_dict={a: x})
        {a_: x, b_: y, c_: z}

        the only function of the repl_dict is to return it in the
        result, e.g. if you omit it:

        >>> (x + sin(b)*c)._matches_commutative(x + sin(y)*z)
        {b_: y, c_: z}

        the "a: x" is not returned in the result, but otherwise it is
        equivalent.

        """
        # make sure expr is Expr if pattern is Expr
        from .expr import Add, Expr
        from .mul import Mul
        if isinstance(self, Expr) and not isinstance(expr, Expr):
            return

        # handle simple patterns
        if self == expr:
            return repl_dict

        d = self._matches_simple(expr, repl_dict)
        if d is not None:
            return d

        # eliminate exact part from pattern: (2+a+w1+w2)._matches(expr) -> (w1+w2)._matches(expr-a-2)
        from .function import WildFunction
        from .symbol import Wild
        wild_part = []
        exact_part = []
        for p in ordered(self.args):
            if p.has(Wild, WildFunction) and (not expr.has(p)):
                # not all Wild should stay Wilds, for example:
                # (w2+w3)._matches(w1) -> (w1+w3)._matches(w1) -> w3._matches(0)
                wild_part.append(p)
            else:
                exact_part.append(p)

        if exact_part:
            exact = self.func(*exact_part)
            free = expr.free_symbols
            if free and (exact.free_symbols - free):
                # there are symbols in the exact part that are not
                # in the expr; but if there are no free symbols, let
                # the matching continue
                return
            newpattern = self.func(*wild_part)
            newexpr = self._combine_inverse(expr, exact)
            if all(isinstance(e, AssocOp) for e in [expr, newexpr]):
                if newexpr.count_ops() > expr.count_ops():
                    return
            return newpattern._matches(newexpr, repl_dict)

        # now to real work ;)
        i = 0
        saw = set()
        while expr not in saw:
            saw.add(expr)
            expr_list = (self.identity,) + tuple(ordered(self.make_args(expr)))
            for last_op in reversed(expr_list):
                for w in reversed(wild_part):
                    d1 = w._matches(last_op, repl_dict)
                    if d1 is not None:
                        d2 = self.xreplace(d1)._matches(expr, d1)
                        if d2 is not None:
                            return d2

            if i == 0:
                if self.is_Mul:
                    # make e**i look like Mul
                    if expr.is_Pow and expr.exp.is_Integer:
                        if expr.exp > 0:
                            expr = Mul(*[expr.base, expr.base**(expr.exp - 1)], evaluate=False)
                        else:
                            expr = Mul(*[1/expr.base, expr.base**(expr.exp + 1)], evaluate=False)
                        i += 1
                        continue

                elif self.is_Add:
                    # make i*e look like Add
                    c, e = expr.as_coeff_Mul()
                    if abs(c) > 1:
                        if c > 0:
                            expr = Add(*[e, (c - 1)*e], evaluate=False)
                        else:
                            expr = Add(*[-e, (c + 1)*e], evaluate=False)
                        i += 1
                        continue

                    # try collection on non-Wild symbols
                    from ..simplify.radsimp import collect
                    was = expr
                    did = set()
                    for w in reversed(wild_part):
                        c, w = w.as_coeff_mul(Wild)
                        free = c.free_symbols - did
                        if free:
                            did.update(free)
                            expr = collect(expr, free)
                    if expr != was:
                        i += 0
                        continue

                else:
                    raise NotImplementedError

                break  # if we didn't continue, there is nothing more to do

        return

    def _has_matcher(self):
        """Helper for .has()."""
        def _ncsplit(expr):
            # this is not the same as args_cnc because here
            # we don't assume expr is a Mul -- hence deal with args --
            # and always return a set.
            cpart, ncpart = [], []
            for arg in expr.args:
                if arg.is_commutative:
                    cpart.append(arg)
                else:
                    ncpart.append(arg)
            return set(cpart), ncpart

        c, nc = _ncsplit(self)
        cls = self.__class__

        def is_in(expr):
            if expr == self:
                return True
            elif isinstance(expr, cls):
                _c, _nc = _ncsplit(expr)
                if (c & _c) == c:
                    if not nc:
                        return True
                    elif len(nc) <= len(_nc):
                        for i in range(len(_nc) - len(nc)):
                            if _nc[i:i + len(nc)] == nc:
                                return True
            return False
        return is_in

    def _eval_evalf(self, prec):
        """
        Evaluate the parts of self that are numbers; if the whole thing
        was a number with no functions it would have been evaluated, but
        it wasn't so we must judiciously extract the numbers and reconstruct
        the object. This is *not* simply replacing numbers with evaluated
        numbers. Nunmbers should be handled in the largest pure-number
        expression as possible. So the code below separates ``self`` into
        number and non-number parts and evaluates the number parts and
        walks the args of the non-number part recursively (doing the same
        thing).

        """
        from .add import Add
        from .mul import Mul
        from .symbol import Symbol
        from .function import AppliedUndef

        if isinstance(self, (Mul, Add)):
            x, tail = self.as_independent(Symbol, AppliedUndef)
            # if x is an AssocOp Function then the _evalf below will
            # call _eval_evalf (here) so we must break the recursion
            if not (tail is self.identity or
                    isinstance(x, AssocOp) and x.is_Function):
                # here, we have a number so we just call to _evalf with prec;
                # prec is not the same as n, it is the binary precision so
                # that's why we don't call to evalf.
                x = x._evalf(prec) if x is not self.identity else self.identity
                args = []
                for a in self.func.make_args(tail):
                    # here we call to _eval_evalf since we don't know what we
                    # are dealing with and all other _eval_evalf routines should
                    # be doing the same thing (i.e. taking binary prec and
                    # finding the evalf-able args)
                    newa = a._eval_evalf(prec)
                    if newa is None:
                        args.append(a)
                    else:
                        args.append(newa)
                if not _aresame(tuple(args), self.func.make_args(tail)):
                    tail = self.func(*args)
                return self.func(x, tail)

        # this is the same as above, but there were no pure-number args to
        # deal with
        args = []
        for a in self.args:
            newa = a.evalf(prec, strict=False)
            args.append(newa)
        if not _aresame(tuple(args), self.args):
            return self.func(*args)
        return self

    @classmethod
    def make_args(cls, expr):
        """
        Return a sequence of elements `args` such that cls(*args) == expr

        >>> Mul.make_args(x*y)
        (x, y)
        >>> Add.make_args(x*y)
        (x*y,)
        >>> set(Add.make_args(x*y + y))
        {y, x*y}

        """
        if isinstance(expr, cls):
            return expr.args
        else:
            return expr,


class ShortCircuit(Exception):
    pass


class LatticeOp(AssocOp):
    """
    Join/meet operations of an algebraic lattice[1].

    These binary operations are associative (op(op(a, b), c) = op(a, op(b, c))),
    commutative (op(a, b) = op(b, a)) and idempotent (op(a, a) = op(a) = a).
    Common examples are AND, OR, Union, Intersection, max or min. They have an
    identity element (op(identity, a) = a) and an absorbing element
    conventionally called zero (op(zero, a) = zero).

    This is an abstract base class, concrete derived classes must declare
    attributes zero and identity. All defining properties are then respected.

    >>> class my_join(LatticeOp):
    ...     zero = Integer(0)
    ...     identity = Integer(1)
    >>> my_join(2, 3) == my_join(3, 2)
    True
    >>> my_join(2, my_join(3, 4)) == my_join(2, 3, 4)
    True
    >>> my_join(0, 1, 4, 2, 3, 4)
    0
    >>> my_join(1, 2)
    2

    References
    ==========

    * https://en.wikipedia.org/wiki/Lattice_%28order%29

    """

    is_commutative = True

    def __new__(cls, *args, **options):
        args = (sympify(arg, strict=True) for arg in args)

        if options.pop('evaluate', global_evaluate[0]):
            try:
                _args = frozenset(cls._new_args_filter(args))
            except ShortCircuit:
                return sympify(cls.zero)
            if not _args:
                return sympify(cls.identity)
            elif len(_args) == 1:
                return set(_args).pop()
        else:
            _args = frozenset(args)

        obj = super(AssocOp, cls).__new__(cls, _args)
        obj._argset = _args
        return obj

    @classmethod
    def _new_args_filter(cls, arg_sequence, call_cls=None):
        """Generator filtering args."""
        ncls = call_cls or cls
        for arg in arg_sequence:
            if arg == ncls.zero:
                raise ShortCircuit(arg)
            elif arg == ncls.identity:
                continue
            elif arg.func == ncls:
                for x in arg.args:
                    yield x
            else:
                yield arg

    @classmethod
    def make_args(cls, expr):
        """
        Return a sequence of elements `args` such that cls(*args) == expr

        >>> Mul.make_args(x*y)
        (x, y)
        >>> Add.make_args(x*y)
        (x*y,)
        >>> set(Add.make_args(x*y + y))
        {y, x*y}

        """
        if isinstance(expr, cls):
            return expr._argset
        else:
            return frozenset([expr])

    @property
    @cacheit
    def args(self):
        return tuple(ordered(self._argset))
