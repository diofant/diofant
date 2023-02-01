"""
Logic expressions handling

Notes
-----
At present this is mainly needed for facts.py, feel free however to improve
this stuff for general purpose.

"""

from __future__ import annotations


def _fuzzy_group(args, quick_exit=False):
    """
    Return True if all args are True, None if there is any None else False
    unless ``quick_exit`` is True (then return None as soon as a second False
    is seen.

     ``_fuzzy_group`` is like ``fuzzy_and`` except that it is more
    conservative in returning a False, waiting to make sure that all
    arguments are True or False and returning None if any arguments are
    None. It also has the capability of permiting only a single False and
    returning None if more than one is seen. For example, the presence of a
    single transcendental amongst rationals would indicate that the group is
    no longer rational; but a second transcendental in the group would make the
    determination impossible.

    Examples
    ========

    By default, multiple Falses mean the group is broken:

    >>> _fuzzy_group([False, False, True])
    False

    If multiple Falses mean the group status is unknown then set
    `quick_exit` to True so None can be returned when the 2nd False is seen:

    >>> _fuzzy_group([False, False, True], quick_exit=True)

    But if only a single False is seen then the group is known to
    be broken:

    >>> _fuzzy_group([False, True, True], quick_exit=True)
    False

    """
    saw_other = False
    for a in args:
        if a is True:
            continue
        if a is None:
            return
        if quick_exit and saw_other:
            return
        saw_other = True
    return not saw_other


def fuzzy_bool(x):
    """
    Return True, False or None according to x.

    Whereas bool(x) returns True or False, fuzzy_bool allows
    for the None value.

    """
    if x is not None:
        return bool(x)


def fuzzy_and(args):
    """
    Return True (all True), False (any False) or None.

    Examples
    ========

    If you had a list of objects to test the commutivity of
    and you want the fuzzy_and logic applied, passing an
    iterator will allow the commutativity to only be computed
    as many times as necessary. With this list, False can be
    returned after analyzing the first symbol:

    >>> syms = [Dummy(commutative=False), Dummy()]
    >>> fuzzy_and(s.is_commutative for s in syms)
    False

    That False would require less work than if a list of pre-computed
    items was sent:

    >>> fuzzy_and([s.is_commutative for s in syms])
    False

    """
    rv = True
    for ai in args:
        ai = fuzzy_bool(ai)
        if ai is False:
            return False
        if rv:  # this will stop updating if a None is ever trapped
            rv = ai
    return rv


def fuzzy_not(v):
    """
    Not in fuzzy logic

    Return None if `v` is None else `not v`.

    Examples
    ========

    >>> fuzzy_not(True)
    False
    >>> fuzzy_not(None)
    >>> fuzzy_not(False)
    True

    """
    if v is None:
        return v
    return not v


def fuzzy_or(args):
    """
    Or in fuzzy logic. Returns True (any True), False (all False), or None

    Examples
    ========

    >>> fuzzy_or([True, False])
    True
    >>> fuzzy_or([True, None])
    True
    >>> fuzzy_or([False, False])
    False
    >>> print(fuzzy_or([False, None]))
    None

    See Also
    ========

    fuzzy_and, fuzzy_not

    """
    return fuzzy_not(fuzzy_and(fuzzy_not(i) for i in args))


class Logic:
    """Logical expression."""

    op_2class: dict[str, type[Logic]] = {}

    def __new__(cls, *args):
        obj = object.__new__(cls)
        obj.args = args
        return obj

    def __getnewargs__(self):
        return self.args

    def __hash__(self):
        return hash((type(self).__name__,) + tuple(self.args))

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return self.args == other.args

    def __str__(self):
        return f"{self.__class__.__name__}({', '.join(str(a) for a in self.args)})"

    __repr__ = __str__

    @staticmethod
    def fromstring(text):
        """
        Logic from string with space around & and | but none after ~.

        e.g.

        ~a & b | c

        """
        lexpr = None  # current logical expression
        schedop = None  # scheduled operation
        for term in text.split():
            # operation symbol
            if term in '&|':
                if schedop is not None:
                    raise ValueError(
                        f'double op forbidden: "{term} {schedop}"')
                if lexpr is None:
                    raise ValueError(
                        f'{term} cannot be in the beginning of expression')
                schedop = term
                continue
            if '&' in term or '|' in term:
                raise ValueError('& and | must have space around them')
            if term[0] == '~':
                if len(term) == 1:
                    raise ValueError('do not include space after "~"')
                term = Not(term[1:])

            # already scheduled operation, e.g. '&'
            if schedop:
                lexpr = Logic.op_2class[schedop](lexpr, term)
                schedop = None
                continue

            # this should be atom
            if lexpr is not None:
                raise ValueError(
                    f'missing op between "{lexpr}" and "{term}"')

            lexpr = term

        # let's check that we ended up in correct state
        if schedop is not None:
            raise ValueError(f'premature end-of-expression in "{text}"')
        if lexpr is None:
            raise ValueError(f'"{text}" is empty')

        # everything looks good now
        return lexpr


class AndOr_Base(Logic):
    """Base class for And and Or."""

    def __new__(cls, *args):
        bargs = []
        for a in args:
            if a == cls.op_x_notx:
                return a
            if a == (not cls.op_x_notx):
                continue    # skip this argument
            bargs.append(a)

        args = sorted(set(cls.flatten(bargs)), key=hash)

        for a in args:
            if Not(a) in args:
                return cls.op_x_notx

        if len(args) == 1:
            return args.pop()
        if len(args) == 0:
            return not cls.op_x_notx

        return Logic.__new__(cls, *args)

    @classmethod
    def flatten(cls, args):
        # quick-n-dirty flattening for And and Or
        args_queue = list(args)
        res = []

        while True:
            try:
                arg = args_queue.pop(0)
            except IndexError:
                break
            if isinstance(arg, Logic):
                if isinstance(arg, cls):
                    args_queue.extend(arg.args)
                    continue
            res.append(arg)

        args = tuple(res)
        return args


class And(AndOr_Base):
    """Logical And."""

    op_x_notx = False

    def _eval_propagate_not(self):
        # ~(a&b&c ...) == ~a | ~b | ~c ...
        return Or(*[Not(a) for a in self.args])

    # (a|b|...) & c == (a&c) | (b&c) | ...
    def expand(self):

        # first locate Or
        for i, arg in enumerate(self.args):
            if isinstance(arg, Or):
                arest = self.args[:i] + self.args[i + 1:]

                orterms = [And(*(arest + (a,))) for a in arg.args]
                for j, orterm in enumerate(orterms):
                    if isinstance(orterm, Logic):
                        orterms[j] = orterm.expand()

                return Or(*orterms)

        return self


class Or(AndOr_Base):
    """Logical Or."""

    op_x_notx = True

    def _eval_propagate_not(self):
        # ~(a|b|c ...) == ~a & ~b & ~c ...
        return And(*[Not(a) for a in self.args])


class Not(Logic):
    """Logical Not."""

    def __new__(cls, arg):
        if isinstance(arg, str):
            return Logic.__new__(cls, arg)

        if isinstance(arg, bool):
            return not arg
        if isinstance(arg, Not):
            return arg.args[0]

        if isinstance(arg, Logic):
            # XXX this is a hack to expand right from the beginning
            arg = arg._eval_propagate_not()
            return arg

        raise ValueError(f'Not: unknown argument {arg!r}')

    @property
    def arg(self):
        return self.args[0]


Logic.op_2class['&'] = And
Logic.op_2class['|'] = Or
Logic.op_2class['~'] = Not
