"""Options manager for :class:`~diofant.polys.polytools.Poly` and public API functions."""

from __future__ import annotations

import re

from ..core import Basic, I
from ..core.sympify import sympify
from ..utilities import has_dups, numbered_symbols, topological_sort
from .polyerrors import FlagError, GeneratorsError, OptionError


__all__ = 'Options', 'Order'


class Option:
    """Base class for all kinds of options."""

    option: str

    is_Flag = False

    requires: list[str] = []
    excludes: list[str] = []

    after: list[str] = []
    before: list[str] = []

    @classmethod
    def default(cls):
        return

    @classmethod
    def preprocess(cls, option):
        return   # pragma: no cover

    @classmethod
    def postprocess(cls, options):
        return


class Flag(Option):
    """Base class for all kinds of flags."""

    is_Flag = True


class BooleanOption(Option):
    """An option that must have a boolean value or equivalent assigned."""

    @classmethod
    def preprocess(cls, option):
        if option in [True, False]:
            return bool(option)
        raise OptionError(f"'{cls.option}' must have a boolean value "
                          f'assigned, got {option}')


class OptionType(type):
    """Base type for all options that does registers options."""

    def __init__(cls, *args, **kwargs):
        """Initialize cls."""
        super().__init__(cls)

        @property
        def getter(a):
            try:
                return a[cls.option]
            except KeyError:
                return cls.default()

        setattr(Options, cls.option, getter)
        Options.__options__[cls.option] = cls


class Options(dict):
    """
    Options manager for polynomial manipulation module.

    Examples
    ========

    >>> Options((x, y, z), {'domain': 'ZZ'})
    {'auto': False, 'domain': ZZ, 'gens': (x, y, z)}

    >>> build_options((x, y, z), {'domain': 'ZZ'})
    {'auto': False, 'domain': ZZ, 'gens': (x, y, z)}

    **Options**

    * Expand --- boolean option
    * Gens --- option
    * Wrt --- option
    * Sort --- option
    * Order --- option
    * Field --- boolean option
    * Greedy --- boolean option
    * Domain --- option
    * Split --- boolean option
    * Gaussian --- boolean option
    * Extension --- option
    * Modulus --- option
    * Symmetric --- boolean option
    * Strict --- boolean option

    **Flags**

    * Auto --- boolean flag
    * Frac --- boolean flag
    * Formal --- boolean flag
    * Polys --- boolean flag
    * Include --- boolean flag
    * All --- boolean flag
    * Gen --- flag

    """

    __order__: list[str] | None = None
    __options__: dict[str, type[Option]] = {}

    def __init__(self, gens, args, flags=None, strict=False):
        """Initialize self."""
        dict.__init__(self)

        if gens and args.get('gens', ()):
            raise OptionError("both '*gens' and keyword "
                              "argument 'gens' supplied")
        if gens:
            args = dict(args)
            args['gens'] = gens

        defaults = args.pop('defaults', {})

        def preprocess_options(args):
            for option, value in args.items():
                try:
                    cls = self.__options__[option]
                except KeyError as exc:
                    raise OptionError(f"'{option}' is not a "
                                      'valid option') from exc

                if issubclass(cls, Flag):
                    if strict and (flags is None or option not in flags):
                        raise OptionError(f"'{option}' flag is not "
                                          'allowed in this context')

                if value is not None:
                    self[option] = cls.preprocess(value)

        preprocess_options(args)

        for key in dict(defaults):
            if key in self:
                del defaults[key]
            else:
                for option in self:
                    cls = self.__options__[option]

                    if key in cls.excludes:
                        del defaults[key]
                        break

        preprocess_options(defaults)

        for option in self:
            cls = self.__options__[option]

            for exclude_option in cls.excludes:
                if self.get(exclude_option) is not None:
                    raise OptionError(f"'{option}' option is not allowed together with '{exclude_option}'")

        for option in self.__order__:  # pylint: disable=not-an-iterable
            self.__options__[option].postprocess(self)

    @classmethod
    def _init_dependencies_order(cls):
        """Resolve the order of options' processing."""
        if cls.__order__ is None:
            vertices, edges = [], set()

            for name, option in cls.__options__.items():
                vertices.append(name)

                for _name in option.after:
                    edges.add((_name, name))

                for _name in option.before:
                    edges.add((name, _name))

            try:
                cls.__order__ = topological_sort((vertices, list(edges)))
            except ValueError as exc:
                raise RuntimeError('cycle detected in diofant.polys'
                                   ' options framework') from exc

    def clone(self, updates={}):
        """Clone ``self`` and update specified options."""
        obj = dict.__new__(self.__class__)

        for option, value in self.items():
            obj[option] = value

        for option, value in updates.items():
            obj[option] = value

        return obj

    def __setattr__(self, attr, value):
        if attr in self.__options__:
            self[attr] = value
        else:
            super().__setattr__(attr, value)

    @property
    def args(self):
        args = {}

        for option, value in self.items():
            if value is not None and option != 'gens':
                cls = self.__options__[option]

                if not issubclass(cls, Flag):
                    args[option] = value

        return args

    @property
    def options(self):
        options = {}
        for option, cls in self.__options__.items():
            if not issubclass(cls, Flag):
                options[option] = getattr(self, option)
        return options

    @property
    def flags(self):
        flags = {}

        for option, cls in self.__options__.items():
            if issubclass(cls, Flag):
                flags[option] = getattr(self, option)

        return flags


class Expand(BooleanOption, metaclass=OptionType):
    """``expand`` option to polynomial manipulation functions."""

    option = 'expand'

    @classmethod
    def default(cls):
        return True


class Gens(Option, metaclass=OptionType):
    """``gens`` option to polynomial manipulation functions."""

    option = 'gens'

    @classmethod
    def default(cls):
        return ()

    @classmethod
    def preprocess(cls, option):
        if isinstance(option, Basic):
            option = option,

        if option == (None,):
            return ()
        if has_dups(option):
            raise GeneratorsError(f'duplicated generators: {option}')
        if any(gen.is_commutative is False for gen in option):
            raise GeneratorsError(f'non-commutative generators: {option}')
        return tuple(option)


class Wrt(Option, metaclass=OptionType):
    """``wrt`` option to polynomial manipulation functions."""

    option = 'wrt'

    _re_split = re.compile(r'\s*,\s*|\s+')

    @classmethod
    def preprocess(cls, option):
        if isinstance(option, Basic):
            return [str(option)]
        if isinstance(option, str):
            option = option.strip()
            if option.endswith(','):
                raise OptionError('Bad input: missing parameter.')
            if not option:
                return []
            return list(cls._re_split.split(option))
        if hasattr(option, '__getitem__'):
            return list(map(str, option))
        raise OptionError("invalid argument for 'wrt' option")


class Sort(Option, metaclass=OptionType):
    """``sort`` option to polynomial manipulation functions."""

    option = 'sort'

    @classmethod
    def default(cls):
        return []

    @classmethod
    def preprocess(cls, option):
        if isinstance(option, str):
            return [gen.strip() for gen in option.split('>')]
        if hasattr(option, '__getitem__'):
            return list(map(str, option))
        raise OptionError("invalid argument for 'sort' option")


class Order(Option, metaclass=OptionType):
    """``order`` option to polynomial manipulation functions."""

    option = 'order'

    @classmethod
    def default(cls):
        from .orderings import lex
        return lex

    @classmethod
    def preprocess(cls, option):
        from .orderings import monomial_key
        return monomial_key(option)


class Field(BooleanOption, metaclass=OptionType):
    """``field`` option to polynomial manipulation functions."""

    option = 'field'

    excludes = ['domain', 'split', 'gaussian']


class Greedy(BooleanOption, metaclass=OptionType):
    """``greedy`` option to polynomial manipulation functions."""

    option = 'greedy'

    excludes = ['domain', 'split', 'gaussian', 'extension', 'modulus']


class Composite(BooleanOption, metaclass=OptionType):
    """``composite`` option to polynomial manipulation functions."""

    option = 'composite'

    @classmethod
    def default(cls):
        return

    excludes = ['domain', 'split', 'gaussian', 'modulus']


class Domain(Option, metaclass=OptionType):
    """``domain`` option to polynomial manipulation functions."""

    option = 'domain'

    excludes = ['field', 'greedy', 'split', 'gaussian', 'extension']

    after = ['gens']

    _re_realfield = re.compile(r'^(R|RR)(_(\d+))?$')
    _re_complexfield = re.compile(r'^(C|CC)(_(\d+))?$')
    _re_finitefield = re.compile(r'^(FF|GF)\((\d+)\)$')
    _re_polynomial = re.compile(r'^(Z|ZZ|Q|QQ)\[(.+)\]$')
    _re_fraction = re.compile(r'^(Z|ZZ|Q|QQ)\((.+)\)$')
    _re_algebraic = re.compile(r'^(Q|QQ)\<(.+)\>$')

    @classmethod
    def preprocess(cls, option):
        from .. import domains
        if isinstance(option, domains.Domain):
            return option
        if isinstance(option, str):
            if option in ['Z', 'ZZ']:
                return domains.ZZ

            if option in ['Q', 'QQ']:
                return domains.QQ

            if option == 'EX':
                return domains.EX

            r = cls._re_realfield.match(option)

            if r is not None:
                _, _, prec = r.groups()

                if prec is None:
                    return domains.RR
                return domains.RealField(int(prec))

            r = cls._re_complexfield.match(option)

            if r is not None:
                _, _, prec = r.groups()

                if prec is None:
                    return domains.CC
                return domains.ComplexField(int(prec))

            r = cls._re_finitefield.match(option)

            if r is not None:
                return domains.FF(int(r.groups()[1]))

            r = cls._re_polynomial.match(option)

            if r is not None:
                ground, gens = r.groups()

                gens = list(map(sympify, gens.split(',')))

                if ground in ['Z', 'ZZ']:
                    return domains.ZZ.inject(*gens)
                return domains.QQ.inject(*gens)

            r = cls._re_fraction.match(option)

            if r is not None:
                ground, gens = r.groups()

                gens = list(map(sympify, gens.split(',')))

                if ground in ['Z', 'ZZ']:
                    return domains.ZZ.inject(*gens).field
                return domains.QQ.inject(*gens).field

            r = cls._re_algebraic.match(option)

            if r is not None:
                gens = list(map(sympify, r.groups()[1].split(',')))
                return domains.QQ.algebraic_field(*gens)

        raise OptionError('expected a valid domain specification, '
                          f'got {option}')

    @classmethod
    def postprocess(cls, options):
        from .. import domains
        from ..domains.compositedomain import CompositeDomain

        if 'gens' in options and 'domain' in options and isinstance(options['domain'], CompositeDomain) and \
                (set(options['domain'].symbols) & set(options['gens'])):
            raise GeneratorsError('ground domain and generators '
                                  'interfere together')
        if ('gens' not in options or not options['gens']) and \
                'domain' in options and options['domain'] == domains.EX:
            raise GeneratorsError('you have to provide generators because'
                                  ' EX domain was requested')


class Split(BooleanOption, metaclass=OptionType):
    """``split`` option to polynomial manipulation functions."""

    option = 'split'

    excludes = ['field', 'greedy', 'domain', 'gaussian', 'extension', 'modulus']

    @classmethod
    def postprocess(cls, options):
        if 'split' in options:
            raise NotImplementedError("'split' option is not implemented yet")


class Gaussian(BooleanOption, metaclass=OptionType):
    """``gaussian`` option to polynomial manipulation functions."""

    option = 'gaussian'

    excludes = ['field', 'greedy', 'domain', 'split', 'extension', 'modulus']

    @classmethod
    def postprocess(cls, options):
        if 'gaussian' in options and options['gaussian'] is True:
            options['extension'] = {I}
            Extension.postprocess(options)


class Extension(Option, metaclass=OptionType):
    """``extension`` option to polynomial manipulation functions."""

    option = 'extension'

    excludes = ['greedy', 'domain', 'split', 'gaussian', 'modulus']

    @classmethod
    def preprocess(cls, option):
        if option == 1:
            return bool(option)
        if option == 0:
            return bool(option)
        if not hasattr(option, '__iter__'):
            option = {option}
        else:
            if not option:
                option = None
            else:
                option = set(option)

        return option

    @classmethod
    def postprocess(cls, options):
        from .. import domains
        if 'extension' in options and options['extension'] not in (True, False):
            options['domain'] = domains.QQ.algebraic_field(
                *options['extension'])


class Modulus(Option, metaclass=OptionType):
    """``modulus`` option to polynomial manipulation functions."""

    option = 'modulus'

    excludes = ['greedy', 'split', 'domain', 'gaussian', 'extension']

    @classmethod
    def preprocess(cls, option):
        option = sympify(option)

        if option.is_Integer and option > 0:
            return int(option)
        raise OptionError(f"'modulus' must a positive integer, got {option}")

    @classmethod
    def postprocess(cls, options):
        from .. import domains
        if 'modulus' in options:
            modulus = options['modulus']
            options['domain'] = domains.FF(modulus)


class Strict(BooleanOption, metaclass=OptionType):
    """``strict`` option to polynomial manipulation functions."""

    option = 'strict'

    @classmethod
    def default(cls):
        return True


class Auto(BooleanOption, Flag, metaclass=OptionType):
    """``auto`` flag to polynomial manipulation functions."""

    option = 'auto'

    after = ['field', 'domain', 'extension', 'gaussian']

    @classmethod
    def default(cls):
        return True

    @classmethod
    def postprocess(cls, options):
        if ('domain' in options or 'field' in options) and 'auto' not in options:
            options['auto'] = False


class Frac(BooleanOption, Flag, metaclass=OptionType):
    """``frac`` option to polynomial manipulation functions."""

    option = 'frac'

    @classmethod
    def default(cls):
        return False


class Formal(BooleanOption, Flag, metaclass=OptionType):
    """``formal`` flag to polynomial manipulation functions."""

    option = 'formal'

    @classmethod
    def default(cls):
        return False


class Polys(BooleanOption, Flag, metaclass=OptionType):
    """``polys`` flag to polynomial manipulation functions."""

    option = 'polys'


class Include(BooleanOption, Flag, metaclass=OptionType):
    """``include`` flag to polynomial manipulation functions."""

    option = 'include'

    @classmethod
    def default(cls):
        return False


class All(BooleanOption, Flag, metaclass=OptionType):
    """``all`` flag to polynomial manipulation functions."""

    option = 'all'

    @classmethod
    def default(cls):
        return False


class Gen(Flag, metaclass=OptionType):
    """``gen`` flag to polynomial manipulation functions."""

    option = 'gen'

    @classmethod
    def default(cls):
        return 0

    @classmethod
    def preprocess(cls, option):
        if isinstance(option, (Basic, int)):
            return option
        raise OptionError("invalid argument for 'gen' option")


class Symbols(Flag, metaclass=OptionType):
    """``symbols`` flag to polynomial manipulation functions."""

    option = 'symbols'

    @classmethod
    def default(cls):
        return numbered_symbols('s', start=1)

    @classmethod
    def preprocess(cls, option):
        if hasattr(option, '__iter__'):
            return iter(option)
        raise OptionError('expected an iterator or '
                          f'iterable container, got {option}')


class Method(Flag, metaclass=OptionType):
    """``method`` flag to polynomial manipulation functions."""

    option = 'method'

    @classmethod
    def preprocess(cls, option):
        if isinstance(option, str):
            return option.lower()
        raise OptionError(f'expected a string, got {option}')


def build_options(gens, args=None):
    """Construct options from keyword arguments or ... options."""
    if args is None:
        gens, args = (), gens

    if len(args) != 1 or 'opt' not in args or gens:
        return Options(gens, args)
    return args['opt']


def allowed_flags(args, flags):
    """
    Allow specified flags to be used in the given context.

    Examples
    ========

    >>> allowed_flags({'domain': ZZ}, [])

    >>> allowed_flags({'domain': ZZ, 'frac': True}, [])
    Traceback (most recent call last):
    ...
    FlagError: 'frac' flag is not allowed in this context

    >>> allowed_flags({'domain': ZZ, 'frac': True}, ['frac'])

    """
    flags = set(flags)

    for arg in args:
        try:
            if Options.__options__[arg].is_Flag and arg not in flags:
                raise FlagError(f"'{arg}' flag is not allowed "
                                'in this context')
        except KeyError as exc:
            raise OptionError(f"'{arg}' is not a valid option") from exc


def set_defaults(options, **defaults):
    """Update options with default values."""
    if 'defaults' not in options:
        options = dict(options)
        options['defaults'] = defaults

    return options


Options._init_dependencies_order()
