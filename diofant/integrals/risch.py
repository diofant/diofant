"""
The Risch Algorithm for transcendental function integration.

The core algorithms for the Risch algorithm are here.  The subproblem
algorithms are in the rde.py and prde.py files for the Risch
Differential Equation solver and the parametric problems solvers,
respectively.  All important information concerning the differential extension
for an integrand is stored in a DifferentialExtension object, which in the code
is usually called DE.  Throughout the code and Inside the DifferentialExtension
object, the conventions/attribute names are that the base domain is QQ and each
differential extension is x, t0, t1, ..., tn-1 = DE.t. DE.x is the variable of
integration (Dx == 1), DE.D is a list of the derivatives of
x, t1, t2, ..., tn-1 = t, DE.T is the list [x, t1, t2, ..., tn-1], DE.t is the
outer-most variable of the differential extension at the given level (the level
can be adjusted using DE.increment_level() and DE.decrement_level()),
k is the field C(x, t0, ..., tn-2), where C is the constant field.  The
numerator of a fraction is denoted by a and the denominator by
d.  If the fraction is named f, fa == numer(f) and fd == denom(f).
Fractions are returned as tuples (fa, fd).  DE.d and DE.t are used to
represent the topmost derivation and extension variable, respectively.
The docstring of a function signifies whether an argument is in k[t], in
which case it will just return a Poly in t, or in k(t), in which case it
will return the fraction (fa, fd). Other variable names probably come
from the names used in Bronstein's book.
"""

import functools
import math

from ..abc import z
from ..core import Dummy, E, Eq, Integer, Lambda, Mul, Pow, Symbol, oo, sympify
from ..functions import (Piecewise, acos, acot, asin, atan, cos, cosh, cot,
                         coth, exp, log, sin, sinh, tan, tanh)
from ..polys import (Poly, PolynomialError, RootSum, cancel, gcd, real_roots,
                     reduced)
from ..utilities import default_sort_key, numbered_symbols, ordered
from .heurisch import _symbols
from .integrals import Integral, integrate


def integer_powers(exprs):
    """
    Rewrites a list of expressions as integer multiples of each other.

    For example, if you have [x, x/2, x**2 + 1, 2*x/3], then you can rewrite
    this as [(x/6) * 6, (x/6) * 3, (x**2 + 1) * 1, (x/6) * 4]. This is useful
    in the Risch integration algorithm, where we must write exp(x) + exp(x/2)
    as (exp(x/2))**2 + exp(x/2), but not as exp(x) + sqrt(exp(x)) (this is
    because only the transcendental case is implemented and we therefore cannot
    integrate algebraic extensions). The integer multiples returned by this
    function for each term are the smallest possible (their content equals 1).

    Returns a list of tuples where the first element is the base term and the
    second element is a list of `(item, factor)` terms, where `factor` is the
    integer multiplicative factor that must multiply the base term to obtain
    the original item.

    The easiest way to understand this is to look at an example:

    >>> integer_powers([x, x/2, x**2 + 1, 2*x/3])
    [(x/6, [(x, 6), (x/2, 3), (2*x/3, 4)]), (x**2 + 1, [(x**2 + 1, 1)])]

    We can see how this relates to the example at the beginning of the
    docstring.  It chose x/6 as the first base term.  Then, x can be written as
    (x/2) * 2, so we get (0, 2), and so on. Now only element (x**2 + 1)
    remains, and there are no other terms that can be written as a rational
    multiple of that, so we get that it can be written as (x**2 + 1) * 1.

    """
    # Here is the strategy:

    # First, go through each term and determine if it can be rewritten as a
    # rational multiple of any of the terms gathered so far.
    # cancel(a/b).is_Rational is sufficient for this.  If it is a multiple, we
    # add its multiple to the dictionary.

    terms = {}
    for term in exprs:
        for j, t in terms.items():
            a = cancel(term/j)
            if a.is_Rational:
                t.append((term, a))
                break
        else:
            terms[term] = [(term, Integer(1))]

    # After we have done this, we have all the like terms together, so we just
    # need to find a common denominator so that we can get the base term and
    # integer multiples such that each term can be written as an integer
    # multiple of the base term, and the content of the integers is 1.

    newterms = {}
    for term, val in terms.items():
        common_denom = functools.reduce(math.lcm, [i.as_numer_denom()[1] for _, i in
                                                   val])
        newterm = term/common_denom
        newmults = [(i, j*common_denom) for i, j in val]
        newterms[newterm] = newmults

    return sorted(newterms.items(), key=lambda item: item[0].sort_key())


class DifferentialExtension:
    """
    A container for all the information relating to a differential extension.

    The attributes of this object are (see also the docstring of __init__):

    - f: The original (Expr) integrand.
    - x: The variable of integration.
    - T: List of variables in the extension.
    - D: List of derivations in the extension; corresponds to the elements of T.
    - fa: Poly of the numerator of the integrand.
    - fd: Poly of the denominator of the integrand.
    - Tfuncs: Lambda() representations of each element of T (except for x).
      For back-substitution after integration.
    - backsubs: A (possibly empty) list of further substitutions to be made on
      the final integral to make it look more like the integrand.
    - E_K: List of the positions of the exponential extensions in T.
    - E_args: The arguments of each of the exponentials in E_K.
    - L_K: List of the positions of the logarithmic extensions in T.
    - L_args: The arguments of each of the logarithms in L_K.
    (See the docstrings of is_deriv_k() and is_log_deriv_k_t_radical() for
    more information on E_K, E_args, L_K, and L_args)
    - cases: List of string representations of the cases of T.
    - t: The top level extension variable, as defined by the current level
      (see level below).
    - d: The top level extension derivation, as defined by the current
      derivation (see level below).
    - case: The string representation of the case of self.d.
    (Note that self.T and self.D will always contain the complete extension,
    regardless of the level.  Therefore, you should ALWAYS use DE.t and DE.d
    instead of DE.T[-1] and DE.D[-1].  If you want to have a list of the
    derivations or variables only up to the current level, use
    DE.D[:len(DE.D) + DE.level + 1] and DE.T[:len(DE.T) + DE.level + 1].  Note
    that, in particular, the derivation() function does this.)

    The following are also attributes, but will probably not be useful other
    than in internal use:
    - newf: Expr form of fa/fd.
    - level: The number (between -1 and -len(self.T)) such that
      self.T[self.level] == self.t and self.D[self.level] == self.d.
      Use the methods self.increment_level() and self.decrement_level() to change
      the current level.

    """

    def __init__(self, f=None, x=None, handle_first='log', dummy=True, extension=None, rewrite_complex=False):
        """
        Tries to build a transcendental extension tower from f with respect to x.

        If it is successful, creates a DifferentialExtension object with, among
        others, the attributes fa, fd, D, T, Tfuncs, and backsubs such that
        fa and fd are Polys in T[-1] with rational coefficients in T[:-1],
        fa/fd == f, and D[i] is a Poly in T[i] with rational coefficients in
        T[:i] representing the derivative of T[i] for each i from 1 to len(T).
        Tfuncs is a list of Lambda objects for back replacing the functions
        after integrating.  Lambda() is only used (instead of lambda) to make
        them easier to test and debug. Note that Tfuncs corresponds to the
        elements of T, except for T[0] == x, but they should be back-substituted
        in reverse order.  backsubs is a (possibly empty) back-substitution list
        that should be applied on the completed integral to make it look more
        like the original integrand.

        If it is unsuccessful, it raises NotImplementedError.

        You can also create an object by manually setting the attributes as a
        dictionary to the extension keyword argument.  You must include at least
        D.  Warning, any attribute that is not given will be set to None. The
        attributes T, t, d, cases, case, x, and level are set automatically and
        do not need to be given.  The functions in the Risch Algorithm will NOT
        check to see if an attribute is None before using it.  This also does not
        check to see if the extension is valid (non-algebraic) or even if it is
        self-consistent.  Therefore, this should only be used for
        testing/debugging purposes.

        """
        # XXX: If you need to debug this function, set the break point here

        if extension:
            if 'D' not in extension:
                raise ValueError('At least the key D must be included with '
                                 'the extension flag to DifferentialExtension.')
            for attr in extension:
                setattr(self, attr, extension[attr])

            self._auto_attrs()

            return
        if f is None or x is None:
            raise ValueError('Either both f and x or a manual extension must '
                             'be given.')

        from .prde import is_deriv_k

        if handle_first not in ['log', 'exp']:
            raise ValueError(f"handle_first must be 'log' or 'exp', not {handle_first!s}.")

        # f will be the original function, self.f might change if we reset
        # (e.g., we pull out a constant from an exponential)
        self.f = f
        self.x = x
        self.reset(dummy=dummy)
        exp_new_extension, log_new_extension = True, True
        if rewrite_complex:
            rewritables = {
                (sin, cos, cot, tan, sinh, cosh, coth, tanh): exp,
                (asin, acos, acot, atan): log,
            }
        # rewrite the trigonometric components
            for candidates, rule in rewritables.items():
                self.newf = self.newf.rewrite(candidates, rule)
        else:
            if any(i.has(x) for i in self.f.atoms(sin, cos, tan, atan, asin, acos)):
                raise NotImplementedError('Trigonometric extensions are not '
                                          'supported (yet!)')

        def update(seq, atoms, func):
            s = set(seq)
            new = atoms - s
            s = atoms.intersection(s)
            s.update(list(filter(func, new)))
            return list(s)

        exps = set()
        pows = set()
        numpows = set()
        sympows = set()
        logs = set()
        symlogs = set()

        while True:
            if self.newf.is_rational_function(*self.T):
                break

            if not exp_new_extension and not log_new_extension:
                # We couldn't find a new extension on the last pass, so I guess
                # we can't do it.
                raise NotImplementedError("Couldn't find an elementary "
                                          f'transcendental extension for {f!s}.  Try using a '
                                          'manual extension with the extension flag.')

            # Pre-preparsing.
            #################
            # Get all exp arguments, so we can avoid ahead of time doing
            # something like t1 = exp(x), t2 = exp(x/2) == sqrt(t1).

            # Things like sqrt(exp(x)) do not automatically simplify to
            # exp(x/2), so they will be viewed as algebraic.  The easiest way
            # to handle this is to convert all instances of (a**b)**Rational
            # to a**(Rational*b) before doing anything else.  Note that the
            # _exp_part code can generate terms of this form, so we do need to
            # do this at each pass (or else modify it to not do that).

            ratpows = [i for i in self.newf.atoms(Pow)
                       if (i.base.is_Pow and i.exp.is_Rational)]

            ratpows_repl = [
                (i, i.base.base**(i.exp*i.base.exp)) for i in ratpows]
            self.backsubs += [(j, i) for i, j in ratpows_repl]
            self.newf = self.newf.xreplace(dict(ratpows_repl))

            # To make the process deterministic, the args are sorted
            # so that functions with smaller op-counts are processed first.
            # Ties are broken with the default_sort_key.

            # XXX Although the method is deterministic no additional work
            # has been done to guarantee that the simplest solution is
            # returned and that it would be affected be using different
            # variables. Though it is possible that this is the case
            # one should know that it has not been done intentionally so
            # further improvements may possible.

            # TODO: This probably doesn't need to be completely recomputed at
            # each pass.
            exps = update(exps, {a for a in self.newf.atoms(Pow) if a.base is E},
                          lambda i: i.exp.is_rational_function(*self.T) and
                          i.exp.has(*self.T))
            pows = update(pows, {a for a in self.newf.atoms(Pow) if a.base is not E},
                          lambda i: i.exp.is_rational_function(*self.T) and
                          i.exp.has(*self.T))
            numpows = update(numpows, set(pows),
                             lambda i: not i.base.has(*self.T))
            sympows = update(sympows, set(pows) - set(numpows),
                             lambda i: i.base.is_rational_function(*self.T) and
                             not i.exp.is_Integer)

            # The easiest way to deal with non-base E powers is to convert them
            # into base E, integrate, and then convert back.
            for i in ordered(pows):
                old = i
                new = exp(i.exp*log(i.base))
                # If exp is ever changed to automatically reduce exp(x*log(2))
                # to 2**x, then this will break.  The solution is to not change
                # exp to do that :)
                if i in sympows:
                    if i.exp.is_Rational:
                        raise NotImplementedError('Algebraic extensions are '
                                                  f'not supported ({i!s}).')
                    # We can add a**b only if log(a) in the extension, because
                    # a**b == exp(b*log(a)).
                    basea, based = frac_in(i.base, self.t)
                    A = is_deriv_k(basea, based, self)
                    if A is None:
                        # Nonelementary monomial (so far)

                        # TODO: Would there ever be any benefit from just
                        # adding log(base) as a new monomial?
                        # ANSWER: Yes, otherwise we can't integrate x**x (or
                        # rather prove that it has no elementary integral)
                        # without first manually rewriting it as exp(x*log(x))
                        self.newf = self.newf.xreplace({old: new})
                        self.backsubs += [(new, old)]
                        log_new_extension = self._log_part([log(i.base)],
                                                           dummy=dummy)
                        exps = update(exps, {a for a in self.newf.atoms(Pow)
                                             if a.base is E},
                                      lambda i: (i.exp.is_rational_function(*self.T) and
                                                 i.exp.has(*self.T)))
                        continue
                    _, u, const = A
                    newterm = exp(i.exp*(log(const) + u))
                    # Under the current implementation, exp kills terms
                    # only if they are of the form a*log(x), where a is a
                    # Number.  This case should have already been killed by the
                    # above tests.  Again, if this changes to kill more than
                    # that, this will break, which maybe is a sign that you
                    # shouldn't be changing that.  Actually, if anything, this
                    # auto-simplification should be removed.  See
                    # http://groups.google.com/group/sympy/browse_thread/thread/a61d48235f16867f

                    self.newf = self.newf.xreplace({i: newterm})

                elif i not in numpows:
                    continue
                else:
                    # i in numpows
                    newterm = new
                # TODO: Just put it in self.Tfuncs
                self.backsubs.append((new, old))
                self.newf = self.newf.xreplace({old: newterm})
                exps.append(newterm)

            atoms = self.newf.atoms(log)
            logs = update(logs, atoms,
                          lambda i: i.args[0].is_rational_function(*self.T) and
                          i.args[0].has(*self.T))
            symlogs = update(symlogs, atoms,
                             lambda i: i.has(*self.T) and i.args[0].is_Pow and
                             i.args[0].base.is_rational_function(*self.T) and
                             not i.args[0].base is E and
                             not i.args[0].exp.is_Integer)

            # We can handle things like log(x**y) by converting it to y*log(x)
            # This will fix not only symbolic exponents of the argument, but any
            # non-Integer exponent, like log(sqrt(x)).  The exponent can also
            # depend on x, like log(x**x).
            for i in ordered(symlogs):
                # Unlike in the exponential case above, we do not ever
                # potentially add new monomials (above we had to add log(a)).
                # Therefore, there is no need to run any is_deriv functions
                # here.  Just convert log(a**b) to b*log(a) and let
                # log_new_extension() handle it from there.
                lbase = log(i.args[0].base)
                logs.append(lbase)
                new = i.args[0].exp*lbase
                self.newf = self.newf.xreplace({i: new})
                self.backsubs.append((new, i))

            # remove any duplicates
            logs = sorted(set(logs), key=default_sort_key)

            if handle_first == 'exp' or not log_new_extension:
                exp_new_extension = self._exp_part(exps, dummy=dummy)
                if exp_new_extension is None:
                    # reset and restart
                    self.f = self.newf
                    self.reset(dummy=dummy)
                    exp_new_extension = True
                    continue

            if handle_first == 'log' or not exp_new_extension:
                log_new_extension = self._log_part(logs, dummy=dummy)

        self.fa, self.fd = frac_in(self.newf, self.t)
        self._auto_attrs()

    def __getattr__(self, attr):
        # Avoid AttributeErrors when debugging
        if attr not in ('f', 'x', 'T', 'D', 'fa', 'fd', 'Tfuncs', 'backsubs',
                        'E_K', 'E_args', 'L_K', 'L_args', 'cases', 'case', 't',
                        'd', 'newf', 'level', 'ts'):
            raise AttributeError(f'{self!r} has no attribute {attr!r}')

    def _auto_attrs(self):
        """Set attributes that are generated automatically."""
        if not self.T:
            # i.e., when using the extension flag and T isn't given
            self.T = [i.gen for i in self.D]
        if not self.x:
            self.x = self.T[0]
        self.cases = [get_case(d, t) for d, t in zip(self.D, self.T)]
        self.level = -1
        self.t = self.T[self.level]
        self.d = self.D[self.level]
        self.case = self.cases[self.level]

    def _exp_part(self, exps, dummy=True):
        """
        Try to build an exponential extension.

        Returns True if there was a new extension, False if there was no new
        extension but it was able to rewrite the given exponentials in terms
        of the existing extension, and None if the entire extension building
        process should be restarted.  If the process fails because there is no
        way around an algebraic extension (e.g., exp(log(x)/2)), it will raise
        NotImplementedError.

        """
        from .prde import is_log_deriv_k_t_radical

        new_extension = False
        restart = False
        expargs = [i.exp for i in exps]
        ip = integer_powers(expargs)
        for arg, others in ip:
            # Minimize potential problems with algebraic substitution
            others.sort(key=lambda i: i[1])

            arga, argd = frac_in(arg, self.t)
            A = is_log_deriv_k_t_radical(arga, argd, self)

            if A is not None:
                ans, u, n, const = A
                # if n is 1 or -1, it's algebraic, but we can handle it
                if n == -1:
                    # This probably will never happen, because
                    # Rational.as_numer_denom() returns the negative term in
                    # the numerator.  But in case that changes, reduce it to
                    # n == 1.
                    n = 1
                    u **= -1
                    const *= -1
                    ans = [(i, -j) for i, j in ans]

                if n == 1:
                    # Example: exp(x + x**2) over QQ(x, exp(x), exp(x**2))
                    self.newf = self.newf.xreplace({exp(arg): exp(const)*Mul(*[
                        u**power for u, power in ans])})
                    self.newf = self.newf.xreplace({exp(p*exparg):
                                                    exp(const*p) * Mul(*[u**power for u, power in ans])
                                                    for exparg, p in others})
                    # TODO: Add something to backsubs to put exp(const*p)
                    # back together.

                    continue

                # Bad news: we have an algebraic radical.  But maybe we
                # could still avoid it by choosing a different extension.
                # For example, integer_powers() won't handle exp(x/2 + 1)
                # over QQ(x, exp(x)), but if we pull out the exp(1), it
                # will.  Or maybe we have exp(x + x**2/2), over
                # QQ(x, exp(x), exp(x**2)), which is exp(x)*sqrt(exp(x**2)),
                # but if we use QQ(x, exp(x), exp(x**2/2)), then they will
                # all work.
                #
                # So here is what we do: If there is a non-zero const, pull
                # it out and retry.  Also, if len(ans) > 1, then rewrite
                # exp(arg) as the product of exponentials from ans, and
                # retry that.  If const == 0 and len(ans) == 1, then we
                # assume that it would have been handled by either
                # integer_powers() or n == 1 above if it could be handled,
                # so we give up at that point.  For example, you can never
                # handle exp(log(x)/2) because it equals sqrt(x).

                if const or len(ans) > 1:
                    rad = Mul(*[term**(power/n) for term, power in ans])
                    self.newf = self.newf.xreplace({exp(p*exparg):
                                                    exp(const*p)*rad for exparg, p in others})
                    self.newf = self.newf.xreplace(dict(zip(reversed(self.T),
                                                            reversed([f(self.x) for f in self.Tfuncs]))))
                    restart = True
                    break
                # TODO: give algebraic dependence in error string
                raise NotImplementedError('Cannot integrate over '
                                          'algebraic extensions.')

            arga, argd = frac_in(arg, self.t)
            darga = (argd*derivation(Poly(arga, self.t), self) -
                     arga*derivation(Poly(argd, self.t), self))
            dargd = argd**2
            darga, dargd = darga.cancel(dargd, include=True)
            darg = darga.as_expr()/dargd.as_expr()
            self.t = next(self.ts)
            self.T.append(self.t)
            self.E_args.append(arg)
            self.E_K.append(len(self.T) - 1)
            self.D.append(darg.as_poly(self.t,
                                       expand=False)*Poly(self.t, self.t,
                                                          expand=False))
            if dummy:
                i = Dummy('i')
            else:
                i = Symbol('i')
            self.Tfuncs = self.Tfuncs + [Lambda(i, exp(arg.subs({self.x: i})))]
            self.newf = self.newf.xreplace({exp(exparg): self.t**p
                                            for exparg, p in others})
            new_extension = True

        if not restart:
            return new_extension

    def _log_part(self, logs, dummy=True):
        """
        Try to build a logarithmic extension.

        Returns True if there was a new extension and False if there was no new
        extension but it was able to rewrite the given logarithms in terms
        of the existing extension.  Unlike with exponential extensions, there
        is no way that a logarithm is not transcendental over and cannot be
        rewritten in terms of an already existing extension in a non-algebraic
        way, so this function does not ever return None or raise
        NotImplementedError.

        """
        from .prde import is_deriv_k

        new_extension = False
        logargs = [i.args[0] for i in logs]
        for arg in ordered(logargs):
            # The log case is easier, because whenever a logarithm is algebraic
            # over the base field, it is of the form a1*t1 + ... an*tn + c,
            # which is a polynomial, so we can just replace it with that.
            # In other words, we don't have to worry about radicals.
            arga, argd = frac_in(arg, self.t)
            A = is_deriv_k(arga, argd, self)
            if A is not None:
                _, u, const = A
                newterm = log(const) + u
                self.newf = self.newf.xreplace({log(arg): newterm})
                continue

            arga, argd = frac_in(arg, self.t)
            darga = (argd*derivation(Poly(arga, self.t), self) -
                     arga*derivation(Poly(argd, self.t), self))
            dargd = argd**2
            darg = darga.as_expr()/dargd.as_expr()
            self.t = next(self.ts)
            self.T.append(self.t)
            self.L_args.append(arg)
            self.L_K.append(len(self.T) - 1)
            self.D.append(cancel(darg.as_expr()/arg).as_poly(self.t,
                                                             expand=False))
            if dummy:
                i = Dummy('i')
            else:
                i = Symbol('i')
            self.Tfuncs = self.Tfuncs + [Lambda(i, log(arg.subs({self.x: i})))]
            self.newf = self.newf.xreplace({log(arg): self.t})
            new_extension = True

        return new_extension

    @property
    def _important_attrs(self):
        """
        Returns some of the more important attributes of self.

        Used for testing and debugging purposes.

        The attributes are (fa, fd, D, T, Tfuncs, backsubs, E_K, E_args,
        L_K, L_args).

        """
        # XXX: This might be easier to read as a dict or something
        # Maybe a named tuple.
        return (self.fa, self.fd, self.D, self.T, self.Tfuncs,
                self.backsubs, self.E_K, self.E_args, self.L_K, self.L_args)

    # TODO: Implement __repr__

    def __str__(self):
        return str(self._important_attrs)

    def reset(self, dummy=True):
        """Reset self to an initial state.  Used by __init__."""
        self.t = self.x
        self.T = [self.x]
        self.D = [Poly(1, self.x)]
        self.level = -1
        self.L_K, self.E_K, self.L_args, self.E_args = [], [], [], []
        if dummy:
            self.ts = numbered_symbols('t', cls=Dummy)
        else:
            # For testing
            self.ts = numbered_symbols('t')
        # For various things that we change to make things work that we need to
        # change back when we are done.
        self.backsubs = []
        self.Tfuncs = []
        self.newf = self.f

    def increment_level(self):
        """
        Increment the level of self.

        This makes the working differential extension larger.  self.level is
        given relative to the end of the list (-1, -2, etc.), so we don't need
        do worry about it when building the extension.

        """
        # pylint: disable=attribute-defined-outside-init
        if self.level >= -1:
            raise ValueError('The level of the differential extension cannot '
                             'be incremented any further.')

        self.level += 1
        self.t = self.T[self.level]
        self.d = self.D[self.level]
        self.case = self.cases[self.level]

    def decrement_level(self):
        """
        Decrease the level of self.

        This makes the working differential extension smaller.  self.level is
        given relative to the end of the list (-1, -2, etc.), so we don't need
        do worry about it when building the extension.

        """
        # pylint: disable=attribute-defined-outside-init
        if self.level <= -len(self.T):
            raise ValueError('The level of the differential extension cannot '
                             'be decremented any further.')

        self.level -= 1
        self.t = self.T[self.level]
        self.d = self.D[self.level]
        self.case = self.cases[self.level]


class DecrementLevel:
    """A context manager for decrementing the level of a DifferentialExtension."""

    def __init__(self, DE):
        """Initialize self."""
        self.DE = DE

    def __enter__(self):
        self.DE.decrement_level()

    def __exit__(self, exc_type, exc_value, traceback):
        self.DE.increment_level()


class NonElementaryIntegralExceptionError(Exception):
    """
    Exception used by subroutines within the Risch algorithm to indicate to one
    another that the function being integrated does not have an elementary
    integral in the given differential field.

    """

    # TODO: Rewrite algorithms below to use this (?)

    # TODO: Pass through information about why the integral was nonelementary,
    # and store that in the resulting NonElementaryIntegral somehow.


def gcdex_diophantine(a, b, c):
    """
    Extended Euclidean Algorithm, Diophantine version.

    Given a, b in K[x] and c in (a, b), the ideal generated by a and b,
    return (s, t) such that s*a + t*b == c and either s == 0 or s.degree()
    < b.degree().

    """
    # Extended Euclidean Algorithm (Diophantine Version) pg. 13
    # TODO: This should go in densetools.py.
    # XXX: Bettter name?

    s, g = a.half_gcdex(b)
    q = c.exquo(g)  # Inexact division means c is not in (a, b)
    s = q*s

    if not s.is_zero and s.degree() >= b.degree():
        q, s = s.div(b)

    t = (c - s*a).exquo(b)

    return s, t


def frac_in(f, t, **kwargs):
    """
    Returns the tuple (fa, fd), where fa and fd are Polys in t.

    This is a common idiom in the Risch Algorithm functions, so we abstract
    it out here.  f should be a basic expresion, a Poly, or a tuple (fa, fd),
    where fa and fd are either basic expressions or Polys, and f == fa/fd.
    **kwargs are applied to Poly.

    """
    cancel = kwargs.pop('cancel', False)
    if type(f) is tuple:
        fa, fd = f
        f = fa.as_expr()/fd.as_expr()
    fa, fd = f.as_expr().as_numer_denom()
    fa, fd = fa.as_poly(t, extension=False, **kwargs), fd.as_poly(t, extension=False, **kwargs)
    if cancel:
        fa, fd = fa.cancel(fd, include=True)
    if fa is None or fd is None:
        raise ValueError(f'Could not turn {f} into a fraction in {t}.')
    return fa, fd


def as_poly_1t(p, t, z):
    """
    (Hackish) way to convert an element p of K[t, 1/t] to K[t, z].

    In other words, z == 1/t will be a dummy variable that Poly can handle
    better.

    See issue sympy/sympy#5131.

    Examples
    ========

    >>> p1 = random_poly(x, 10, -10, 10)
    >>> p2 = random_poly(x, 10, -10, 10)
    >>> p = p1 + p2.subs({x: 1/x})
    >>> as_poly_1t(p, x, z).as_expr().subs({z: 1/x}) == p
    True

    """
    # TODO: Use this on the final result.  That way, we can avoid answers like
    # (...)*exp(-x).
    pa, pd = frac_in(p, t, cancel=True)
    if not pd.is_term:
        # XXX: Is there a better Poly exception that we could raise here?
        # Either way, if you see this (from the Risch Algorithm) it indicates
        # a bug.
        raise PolynomialError(f'{p} is not an element of K[{t}, 1/{t}].')
    d = pd.degree(t)
    one_t_part = pa.slice(0, d + 1)
    r = pd.degree() - pa.degree()
    t_part = pa - one_t_part
    t_part = t_part.to_field().exquo(pd)
    # Compute the negative degree parts.
    od = max(-r - one_t_part.degree() if r < 0 < d else 0, 0)
    one_t_part = Poly([0]*od + list(reversed(one_t_part.rep.all_coeffs())),
                      *one_t_part.gens, domain=one_t_part.domain)
    if 0 < r < oo:
        one_t_part *= Poly(t**r, t)

    one_t_part = one_t_part.replace(t, z)  # z will be 1/t
    if pd.coeff_monomial((d,)):
        one_t_part *= Poly(1/pd.coeff_monomial((d,)), z, expand=False)
    ans = t_part.as_poly(t, z, expand=False) + one_t_part.as_poly(t, z,
                                                                  expand=False)

    return ans


def derivation(p, DE, coefficientD=False, basic=False):
    """
    Computes Dp.

    Given the derivation D with D = d/dx and p is a polynomial in t over
    K(x), return Dp.

    If coefficientD is True, it computes the derivation kD
    (kappaD), which is defined as kD(sum(ai*Xi**i, (i, 0, n))) ==
    sum(Dai*Xi**i, (i, 1, n)) (Definition 3.2.2, page 80).  X in this case is
    T[-1], so coefficientD computes the derivative just with respect to T[:-1],
    with T[-1] treated as a constant.

    If basic=True, the returns a Basic expression.  Elements of D can still be
    instances of Poly.

    """
    if basic:
        r = 0
    else:
        r = Poly(0, DE.t, field=True)

    t = DE.t
    if coefficientD:
        if DE.level <= -len(DE.T):
            # 'base' case, the answer is 0.
            return r
        DE.decrement_level()

    D = DE.D[:len(DE.D) + DE.level + 1]
    T = DE.T[:len(DE.T) + DE.level + 1]

    for d, v in zip(D, T):
        pv = p.as_poly(v)
        if pv is None or basic:
            pv = p.as_expr()

        if basic:
            r += d.as_expr()*pv.diff(v)
        else:
            r += (d*pv.diff(v)).as_poly(t, field=True)

    if basic:
        r = cancel(r)
    if coefficientD:
        DE.increment_level()

    return r


def get_case(d, t):
    """
    Returns the type of the derivation d.

    Returns one of {'exp', 'tan', 'base', 'primitive', 'other_linear',
    'other_nonlinear'}.

    """
    if not d.as_expr().has(t):
        if d.is_one:
            return 'base'
        return 'primitive'
    if d.rem(Poly(t, t)).is_zero:
        return 'exp'
    if d.rem(Poly(1 + t**2, t)).is_zero:
        return 'tan'
    if d.degree(t) > 1:
        return 'other_nonlinear'
    return 'other_linear'


def splitfactor(p, DE, coefficientD=False, z=None):
    """
    Splitting factorization.

    Given a derivation D on k[t] and p in k[t], return (p_n, p_s) in
    k[t] x k[t] such that p = p_n*p_s, p_s is special, and each square
    factor of p_n is normal.

    Page. 100

    """
    kinv = [1/x for x in DE.T[:DE.level]]
    if z:
        kinv.append(z)

    One = Poly(1, DE.t, domain=p.domain)
    Dp = derivation(p, DE, coefficientD=coefficientD)
    # XXX: Is this right?
    if p.is_zero:
        return p, One

    if not p.as_expr().has(DE.t):
        s = p.as_poly(*kinv).gcd(Dp.as_poly(*kinv)).as_poly(DE.t)
        n = p.exquo(s)
        return n, s

    if not Dp.is_zero:
        h = p.gcd(Dp).to_field()
        g = p.gcd(p.diff(DE.t)).to_field()
        s = h.exquo(g)

        if s.degree(DE.t) == 0:
            return p, One

        q_split = splitfactor(p.exquo(s), DE, coefficientD=coefficientD)

        return q_split[0], q_split[1]*s
    return p, One


def splitfactor_sqf(p, DE, coefficientD=False, z=None, basic=False):
    """
    Splitting Square-free Factorization

    Given a derivation D on k[t] and p in k[t], returns (N1, ..., Nm)
    and (S1, ..., Sm) in k[t]^m such that p =
    (N1*N2**2*...*Nm**m)*(S1*S2**2*...*Sm**m) is a splitting
    factorization of p and the Ni and Si are square-free and coprime.

    """
    # TODO: This algorithm appears to be faster in every case
    # TODO: Verify this and splitfactor() for multiple extensions
    kkinv = [1/x for x in DE.T[:DE.level]] + DE.T[:DE.level]
    if z:
        kkinv = [z]

    S = []
    N = []
    p_c, p_sqf = p.sqf_list()
    if p_c != 1:
        p_sqf.insert(0, (Poly(p_c, DE.t), 1))
    if p.is_zero:
        return ((p, 1),), ()

    for pi, i in p_sqf:
        Si = pi.as_poly(*kkinv).gcd(derivation(pi, DE,
                                               coefficientD=coefficientD, basic=basic).as_poly(*kkinv)).as_poly(DE.t)
        pi = Poly(pi, DE.t)
        Si = Poly(Si, DE.t)
        Ni = pi.exquo(Si)
        if not Si.is_one:
            S.append((Si, i))
        if not Ni.is_one:
            N.append((Ni, i))

    return tuple(N), tuple(S)


def canonical_representation(a, d, DE):
    """
    Canonical Representation.

    Given a derivation D on k[t] and f = a/d in k(t), return (f_p, f_s,
    f_n) in k[t] x k(t) x k(t) such that f = f_p + f_s + f_n is the
    canonical representation of f (f_p is a polynomial, f_s is reduced
    (has a special denominator), and f_n is simple (has a normal
    denominator).

    """
    # Make d monic
    l = Poly(1/d.LC(), DE.t)
    a *= l
    d *= l

    q, r = a.div(d)
    dn, ds = splitfactor(d, DE)

    b, c = gcdex_diophantine(dn.as_poly(DE.t), ds.as_poly(DE.t), r.as_poly(DE.t))
    b, c = b.as_poly(DE.t), c.as_poly(DE.t)

    return q, (b, ds), (c, dn)


def hermite_reduce(a, d, DE):
    """
    Hermite Reduction - Mack's Linear Version.

    Given a derivation D on k(t) and f = a/d in k(t), returns g, h, r in
    k(t) such that f = Dg + h + r, h is simple, and r is reduced.

    """
    # Make d monic
    l = Poly(1/d.LC(), DE.t)
    a *= l
    d *= l

    fp, fs, fn = canonical_representation(a, d, DE)
    a, d = fn
    l = Poly(1/d.LC(), DE.t)
    a *= l
    d *= l

    ga = Poly(0, DE.t)
    gd = Poly(1, DE.t)

    dd = derivation(d, DE)
    dm = gcd(d, dd).as_poly(DE.t)
    ds, r = d.div(dm)

    while dm.degree(DE.t) > 0:

        ddm = derivation(dm, DE)
        dm2 = gcd(dm, ddm)
        dms, r = dm.div(dm2)
        ds_ddm = ds*ddm
        ds_ddm_dm, r = ds_ddm.div(dm)

        b, c = gcdex_diophantine(-ds_ddm_dm.as_poly(DE.t), dms.as_poly(DE.t), a.as_poly(DE.t))
        b, c = b.as_poly(DE.t), c.as_poly(DE.t)

        db = derivation(b, DE).as_poly(DE.t)
        ds_dms, r = ds.div(dms)
        a = c.as_poly(DE.t) - (db*ds_dms).as_poly(DE.t)

        ga = ga*dm + b*gd
        gd = gd*dm
        ga, gd = ga.cancel(gd, include=True)
        dm = dm2

    d = ds
    q, r = a.div(d)
    ga, gd = ga.cancel(gd, include=True)

    r, d = r.cancel(d, include=True)
    rra = q*fs[1] + fp*fs[1] + fs[0]
    rrd = fs[1]
    rra, rrd = rra.cancel(rrd, include=True)

    return (ga, gd), (r, d), (rra, rrd)


def polynomial_reduce(p, DE):
    """
    Polynomial Reduction.

    Given a derivation D on k(t) and p in k[t] where t is a nonlinear
    monomial over k, return q, r in k[t] such that p = Dq  + r, and
    deg(r) < deg_t(Dt).

    """
    q = Poly(0, DE.t)
    while p.degree(DE.t) >= DE.d.degree(DE.t):
        m = p.degree(DE.t) - DE.d.degree(DE.t) + 1
        q0 = Poly(DE.t**m, DE.t)*Poly(p.as_poly(DE.t).LC() /
                                      (m*DE.d.LC()), DE.t)
        q += q0
        p = p - derivation(q0, DE)

    return q, p


def laurent_series(a, d, F, n, DE):
    """
    Contribution of F to the full partial fraction decomposition of A/D

    Given a field K of characteristic 0 and A,D,F in K[x] with D monic,
    nonzero, coprime with A, and F the factor of multiplicity n in the square-
    free factorization of D, return the principal parts of the Laurent series of
    A/D at all the zeros of F.

    """
    if F.degree() == 0:
        return 0
    Z = _symbols('z', n)
    Z.insert(0, z)
    delta_a = Poly(0, DE.t)
    delta_d = Poly(1, DE.t)

    E = d.quo(F**n)
    ha, hd = (a, E*Poly(z**n, DE.t))
    dF = derivation(F, DE)
    B, _ = gcdex_diophantine(E, F, Poly(1, DE.t))
    C, _ = gcdex_diophantine(dF, F, Poly(1, DE.t))

    # initialization
    F_store = F
    V, DE_D_list, H_list = [], [], []

    for j in range(n):
        # jth derivative of z would be substituted with dfnth/(j+1) where dfnth =(d^n)f/(dx)^n
        F_store = derivation(F_store, DE)
        v = (F_store.as_expr())/(j + 1)
        V.append(v)
        DE_D_list.append(Poly(Z[j + 1], Z[j]))

    DE_new = DifferentialExtension(extension={'D': DE_D_list})  # a differential indeterminate
    for j in range(n):
        zEha = Poly(z**(n + j), DE.t)*E**(j + 1)*ha
        zEhd = hd
        Pa, Pd = cancel((zEha, zEhd))[1], cancel((zEha, zEhd))[2]
        Q = Pa.quo(Pd)
        for i in range(j + 1):
            Q = Q.subs({Z[i]: V[i]})
        Dha = hd*derivation(ha, DE, basic=True) + ha*derivation(hd, DE, basic=True)
        Dha += hd*derivation(ha, DE_new, basic=True) + ha*derivation(hd, DE_new, basic=True)
        Dhd = Poly(j + 1, DE.t)*hd**2
        ha, hd = Dha, Dhd

        Ff, _ = F.div(gcd(F, Q))
        F_stara, F_stard = frac_in(Ff, DE.t)
        if F_stara.degree(DE.t) - F_stard.degree(DE.t) > 0:
            QBC = Poly(Q, DE.t)*B**(1 + j)*C**(n + j)
            H = QBC
            H_list.append(H)
            H = (QBC*F_stard).rem(F_stara)
            alphas = real_roots(F_stara)
            for alpha in list(alphas):
                delta_a = delta_a*Poly((DE.t - alpha)**(n - j), DE.t) + Poly(H.eval(alpha), DE.t)
                delta_d = delta_d*Poly((DE.t - alpha)**(n - j), DE.t)
    return delta_a, delta_d, H_list


def recognize_derivative(a, d, DE, z=None):
    """
    Compute the squarefree factorization of the denominator of f
    and for each Di the polynomial H in K[x] (see Theorem 2.7.1), using the
    LaurentSeries algorithm. Write Di = GiEi where Gj = gcd(Hn, Di) and
    gcd(Ei,Hn) = 1. Since the residues of f at the roots of Gj are all 0, and
    the residue of f at a root alpha of Ei is Hi(a) != 0, f is the derivative of a
    rational function if and only if Ei = 1 for each i, which is equivalent to
    Di | H[-1] for each i.

    """
    flag = True
    a, d = a.cancel(d, include=True)
    _, r = a.div(d)
    _, Sp = splitfactor_sqf(d, DE, coefficientD=True, z=z)

    j = 1
    for (s, _) in Sp:
        *_, H = laurent_series(r, d, s, j, DE)
        g = gcd(d, H[-1]).as_poly()
        if g is not d:
            flag = False
            break
        j = j + 1
    return flag


def recognize_log_derivative(a, d, DE, z=None):
    """
    There exists a v in K(x)* such that f = dv/v
    where f a rational function if and only if f can be written as f = A/D
    where D is squarefree,deg(A) < deg(D), gcd(A, D) = 1,
    and all the roots of the Rothstein-Trager resultant are integers. In that case,
    any of the Rothstein-Trager, Lazard-Rioboo-Trager or Czichowski algorithm
    produces u in K(x) such that du/dx = uf.

    """
    z = z or Dummy('z')
    a, d = a.cancel(d, include=True)
    _, a = a.div(d)

    pz = Poly(z, DE.t)
    Dd = derivation(d, DE)
    q = a - pz*Dd
    r = d.resultant(q)
    r = Poly(r, z)
    _, Sp = splitfactor_sqf(r, DE, coefficientD=True, z=z)

    for s, _ in Sp:
        # TODO also consider the complex roots
        # in case we have complex roots it should turn the flag false
        a = real_roots(s.as_poly(z))

        if any(not j.is_Integer for j in a):
            return False
    return True


def residue_reduce(a, d, DE, z=None, invert=True):
    """
    Lazard-Rioboo-Rothstein-Trager resultant reduction.

    Given a derivation D on k(t) and f in k(t) simple, return g
    elementary over k(t) and a Boolean b in {True, False} such that f -
    Dg in k[t] if b == True or f + h and f + h - Dg do not have an
    elementary integral over k(t) for any h in k<t> (reduced) if b ==
    False.

    Returns (G, b), where G is a tuple of tuples of the form (s_i, S_i),
    such that g = Add(*[RootSum(s_i, lambda z: z*log(S_i(z, t))) for
    S_i, s_i in G]). f - Dg is the remaining integral, which is elementary
    only if b == True, and hence the integral of f is elementary only if
    b == True.

    f - Dg is not calculated in this function because that would require
    explicitly calculating the RootSum.  Use residue_reduce_derivation().

    """
    # TODO: Use log_to_atan() from rationaltools.py
    # If r = residue_reduce(...), then the logarithmic part is given by:
    # sum(RootSum(a[0].as_poly(z),
    #             lambda i: i*log(a[1].as_expr()).subs({z: i})).subs({t: log(x)})
    #     for a in r[0])

    z = z or Dummy('z')
    a, d = a.cancel(d, include=True)
    a, d = a.to_field()*(1/d.LC()), d.to_field()*(1/d.LC())
    kkinv = [1/x for x in DE.T[:DE.level]] + DE.T[:DE.level]

    if a.is_zero:
        return [], True
    _, a = a.div(d)

    pz = Poly(z, DE.t)

    Dd = derivation(d, DE)
    q = a - pz*Dd

    if Dd.degree(DE.t) <= d.degree(DE.t):
        r, R = d.resultant(q, includePRS=True)
    else:
        r, R = q.resultant(d, includePRS=True)

    R_map, H = {}, []
    for i in R:
        R_map[i.degree()] = i

    r = Poly(r, z)
    Np, Sp = splitfactor_sqf(r, DE, coefficientD=True, z=z)

    for s, i in Sp:
        if i == d.degree(DE.t):
            s = Poly(s, z).monic()
            H.append((s, d))
        else:
            h = R_map.get(i)
            if h is None:
                continue
            h_lc = Poly(h.as_poly(DE.t).LC(), DE.t, field=True)

            h_lc_c, h_lc_sqf = h_lc.sqf_list()
            h_lc_sqf.insert(0, (Poly(h_lc_c, DE.t, field=True), 1))

            for a, j in h_lc_sqf:
                h = Poly(h, DE.t, field=True).exquo(Poly(gcd(a, s**j, *kkinv),
                                                         DE.t))

            s = Poly(s, z).monic()

            if invert:
                h_lc = Poly(h.as_poly(DE.t).LC(), DE.t, field=True, expand=False)
                inv, coeffs = h_lc.as_poly(z, field=True).invert(s), [Integer(1)]

                for coeff in h.coeffs()[1:]:
                    L = reduced(inv*coeff, [s])[1]
                    coeffs.append(L.as_expr())

                h = Poly(dict(zip(h.monoms(), coeffs)), DE.t)

            if not s.is_one:
                H.append((s, h))

    b = all(not cancel(i.as_expr()).has(DE.t, z) for i, _ in Np)

    return H, b


def residue_reduce_to_basic(H, DE, z):
    """Converts the tuple returned by residue_reduce() into a Basic expression."""
    # TODO: check what Lambda does with RootOf
    i = Dummy('i')
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    return sum(RootSum(a[0].as_poly(z), Lambda(i, i*log(a[1].as_expr()).subs(
        {z: i}).subs(s))) for a in H)


def residue_reduce_derivation(H, DE, z):
    """
    Computes the derivation of an expression returned by residue_reduce().

    In general, this is a rational function in t, so this returns an
    as_expr() result.

    """
    # TODO: verify that this is correct for multiple extensions
    i = Dummy('i')
    return sympify(sum(RootSum(a[0].as_poly(z),
                       Lambda(i, i*derivation(a[1], DE).as_expr().subs({z: i})/a[1].as_expr().subs({z: i}))) for a in H))


def integrate_primitive_polynomial(p, DE):
    """
    Integration of primitive polynomials.

    Given a primitive monomial t over k, and p in k[t], return q in k[t],
    r in k, and a bool b in {True, False} such that r = p - Dq is in k if b is
    True, or r = p - Dq does not have an elementary integral over k(t) if b is
    False.

    """
    from .prde import limited_integrate

    Zero = Poly(0, DE.t)
    q = Poly(0, DE.t)

    if not p.as_expr().has(DE.t):
        return Zero, p, True

    while True:
        if not p.as_expr().has(DE.t):
            return q, p, True

        Dta, Dtb = frac_in(DE.d, DE.T[DE.level - 1])

        # We had better be integrating the lowest extension (x)
        # with ratint().
        with DecrementLevel(DE):
            a = p.LC()
            aa, ad = frac_in(a, DE.t)

            try:
                (ba, bd), c = limited_integrate(aa, ad, [(Dta, Dtb)], DE)
                if len(c) != 1:
                    raise ValueError('Length of c should  be 1')
            except NonElementaryIntegralExceptionError:
                return q, p, False

        m = p.degree(DE.t)
        q0 = c[0].as_poly(DE.t)*Poly(DE.t**(m + 1)/(m + 1), DE.t) + \
            (ba.as_expr()/bd.as_expr()).as_poly(DE.t)*Poly(DE.t**m, DE.t)

        p = p - derivation(q0, DE)
        q = q + q0


def integrate_primitive(a, d, DE, z=None):
    """
    Integration of primitive functions.

    Given a primitive monomial t over k and f in k(t), return g elementary over
    k(t), i in k(t), and b in {True, False} such that i = f - Dg is in k if b
    is True or i = f - Dg does not have an elementary integral over k(t) if b
    is False.

    This function returns a Basic expression for the first argument.  If b is
    True, the second argument is Basic expression in k to recursively integrate.
    If b is False, the second argument is an unevaluated Integral, which has
    been proven to be nonelementary.

    """
    # XXX: a and d must be canceled, or this might return incorrect results
    z = z or Dummy('z')
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    g1, h, r = hermite_reduce(a, d, DE)
    g2, b = residue_reduce(h[0], h[1], DE, z=z)
    if not b:
        i = cancel(a.as_expr()/d.as_expr() - (g1[1]*derivation(g1[0], DE) -
                                              g1[0]*derivation(g1[1], DE)).as_expr()/(g1[1]**2).as_expr() -
                   residue_reduce_derivation(g2, DE, z))
        i = NonElementaryIntegral(cancel(i).subs(s), DE.x)
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
                residue_reduce_to_basic(g2, DE, z), i, b)

    # h - Dg2 + r
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
                                                                         DE, z) + r[0].as_expr()/r[1].as_expr())
    p = p.as_poly(DE.t)

    q, i, b = integrate_primitive_polynomial(p, DE)

    ret = ((g1[0].as_expr()/g1[1].as_expr() + q.as_expr()).subs(s) +
           residue_reduce_to_basic(g2, DE, z))
    if not b:
        # TODO: This does not do the right thing when b is False
        i = NonElementaryIntegral(cancel(i.as_expr()).subs(s), DE.x)
    else:
        i = cancel(i.as_expr())

    return ret, i, b


def integrate_hyperexponential_polynomial(p, DE, z):
    """
    Integration of hyperexponential polynomials.

    Given a hyperexponential monomial t over k and p in k[t, 1/t], return q in
    k[t, 1/t] and a bool b in {True, False} such that p - Dq in k if b is True,
    or p - Dq does not have an elementary integral over k(t) if b is False.

    """
    from .rde import rischDE

    t1 = DE.t
    dtt = DE.d.exquo(Poly(DE.t, DE.t))
    qa = Poly(0, DE.t)
    qd = Poly(1, DE.t)
    b = True

    if p.is_zero:
        return qa, qd, b

    with DecrementLevel(DE):
        for i in range(-p.degree(z), p.degree(t1) + 1):
            if not i:
                continue
            if i < 0:
                # If you get AttributeError: 'NoneType' object has no attribute 'nth'
                # then this should really not have expand=False
                # But it shouldn't happen because p is already a Poly in t and z
                a = p.as_poly(z, expand=False).coeff_monomial((-i,))
            else:
                # If you get AttributeError: 'NoneType' object has no attribute 'nth'
                # then this should really not have expand=False
                a = p.as_poly(t1, expand=False).coeff_monomial((i,))

            aa, ad = frac_in(a, DE.t, field=True)
            aa, ad = aa.cancel(ad, include=True)
            iDt = Poly(i, t1)*dtt
            iDta, iDtd = frac_in(iDt, DE.t, field=True)
            try:
                va, vd = rischDE(iDta, iDtd, Poly(aa, DE.t), Poly(ad, DE.t), DE)
                va, vd = frac_in((va, vd), t1)
            except NonElementaryIntegralExceptionError:
                b = False
            else:
                qa = qa*vd + va*Poly(t1**i)*qd
                qd *= vd

    return qa, qd, b


def integrate_hyperexponential(a, d, DE, z=None, conds='piecewise'):
    """
    Integration of hyperexponential functions.

    Given a hyperexponential monomial t over k and f in k(t), return g
    elementary over k(t), i in k(t), and a bool b in {True, False} such that
    i = f - Dg is in k if b is True or i = f - Dg does not have an elementary
    integral over k(t) if b is False.

    This function returns a Basic expression for the first argument.  If b is
    True, the second argument is Basic expression in k to recursively integrate.
    If b is False, the second argument is an unevaluated Integral, which has
    been proven to be nonelementary.

    """
    # XXX: a and d must be canceled, or this might return incorrect results
    z = z or Dummy('z')
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    g1, h, r = hermite_reduce(a, d, DE)
    g2, b = residue_reduce(h[0], h[1], DE, z=z)
    if not b:
        i = cancel(a.as_expr()/d.as_expr() - (g1[1]*derivation(g1[0], DE) -
                                              g1[0]*derivation(g1[1], DE)).as_expr()/(g1[1]**2).as_expr() -
                   residue_reduce_derivation(g2, DE, z))
        i = NonElementaryIntegral(cancel(i.subs(s)), DE.x)
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
                residue_reduce_to_basic(g2, DE, z), i, b)

    # p should be a polynomial in t and 1/t, because Sirr == k[t, 1/t]
    # h - Dg2 + r
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
                                                                         DE, z) + r[0].as_expr()/r[1].as_expr())
    pp = as_poly_1t(p, DE.t, z)

    qa, qd, b = integrate_hyperexponential_polynomial(pp, DE, z)

    i = pp.coeff_monomial(1)

    ret = ((g1[0].as_expr()/g1[1].as_expr()).subs(s)
           + residue_reduce_to_basic(g2, DE, z))

    qas = qa.as_expr().subs(s)
    qds = qd.as_expr().subs(s)
    if conds == 'piecewise' and DE.x not in qds.free_symbols:
        # We have to be careful if the exponent is Integer(0)!

        # XXX: Does qd = 0 always necessarily correspond to the exponential
        # equaling 1?
        ret += Piecewise(
            (integrate((p - i).subs({DE.t: 1}).subs(s), DE.x), Eq(qds, 0)),
            (qas/qds, True))
    else:
        ret += qas/qds

    if not b:
        i = p - (qd*derivation(qa, DE) - qa*derivation(qd, DE)).as_expr() /\
            (qd**2).as_expr()
        i = NonElementaryIntegral(cancel(i).subs(s), DE.x)
    return ret, i, b


def integrate_hypertangent_polynomial(p, DE):
    """
    Integration of hypertangent polynomials.

    Given a differential field k such that sqrt(-1) is not in k, a
    hypertangent monomial t over k, and p in k[t], return q in k[t] and
    c in k such that p - Dq - c*D(t**2 + 1)/(t**1 + 1) is in k and p -
    Dq does not have an elementary integral over k(t) if Dc != 0.

    """
    # XXX: Make sure that sqrt(-1) is not in k.
    q, r = polynomial_reduce(p, DE)
    a = DE.d.exquo(Poly(DE.t**2 + 1, DE.t))
    c = Poly(r.coeff_monomial((1,))/(2*a.as_expr()), DE.t)
    return q, c


def integrate_nonlinear_no_specials(a, d, DE, z=None):
    """
    Integration of nonlinear monomials with no specials.

    Given a nonlinear monomial t over k such that Sirr ({p in k[t] | p is
    special, monic, and irreducible}) is empty, and f in k(t), returns g
    elementary over k(t) and a Boolean b in {True, False} such that f - Dg is
    in k if b == True, or f - Dg does not have an elementary integral over k(t)
    if b == False.

    This function is applicable to all nonlinear extensions, but in the case
    where it returns b == False, it will only have proven that the integral of
    f - Dg is nonelementary if Sirr is empty.

    This function returns a Basic expression.

    """
    # TODO: Integral from k?
    # TODO: split out nonelementary integral
    # XXX: a and d must be canceled, or this might not return correct results
    z = z or Dummy('z')
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    g1, h, r = hermite_reduce(a, d, DE)
    g2, b = residue_reduce(h[0], h[1], DE, z=z)
    if not b:
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
                residue_reduce_to_basic(g2, DE, z), b)

    # Because f has no specials, this should be a polynomial in t, or else
    # there is a bug.
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
                                                                         DE, z).as_expr() + r[0].as_expr()/r[1].as_expr()).as_poly(DE.t)
    q1, q2 = polynomial_reduce(p, DE)

    if q2.as_expr().has(DE.t):
        b = False
    else:
        b = True

    ret = (cancel(g1[0].as_expr()/g1[1].as_expr() + q1.as_expr()).subs(s) +
           residue_reduce_to_basic(g2, DE, z))
    return ret, b


class NonElementaryIntegral(Integral):
    """
    Represents a nonelementary Integral.

    If the result of integrate() is an instance of this class, it is
    guaranteed to be nonelementary.  Note that integrate() by default will try
    to find any closed-form solution, even in terms of special functions which
    may themselves not be elementary.  To make integrate() only give
    elementary solutions, or, in the cases where it can prove the integral to
    be nonelementary, instances of this class, use integrate(risch=True).
    In this case, integrate() may raise NotImplementedError if it cannot make
    such a determination.

    integrate() uses the deterministic Risch algorithm to integrate elementary
    functions or prove that they have no elementary integral.  In some cases,
    this algorithm can split an integral into an elementary and nonelementary
    part, so that the result of integrate will be the sum of an elementary
    expression and a NonElementaryIntegral.

    Examples
    ========

    >>> a = integrate(exp(-x**2), x, risch=True)
    >>> a
    Integral(E**(-x**2), x)
    >>> type(a)
    <class 'diofant.integrals.risch.NonElementaryIntegral'>

    >>> expr = (2*log(x)**2 - log(x) - x**2)/(log(x)**3 - x**2*log(x))
    >>> b = integrate(expr, x, risch=True)
    >>> b
    -log(-x + log(x))/2 + log(x + log(x))/2 + Integral(1/log(x), x)
    >>> type(b.atoms(Integral).pop())
    <class 'diofant.integrals.risch.NonElementaryIntegral'>

    """

    # TODO: This is useful in and of itself, because isinstance(result,
    # NonElementaryIntegral) will tell if the integral has been proven to be
    # elementary. But should we do more?  Perhaps a no-op .doit() if
    # elementary=True?  Or maybe some information on why the integral is
    # nonelementary.


def risch_integrate(f, x, extension=None, handle_first='log',
                    separate_integral=False, rewrite_complex=False,
                    conds='piecewise'):
    r"""
    The Risch Integration Algorithm.

    Only transcendental functions are supported.  Currently, only exponentials
    and logarithms are supported, but support for trigonometric functions is
    forthcoming.

    If this function returns an unevaluated Integral in the result, it means
    that it has proven that integral to be nonelementary.  Any errors will
    result in raising NotImplementedError.  The unevaluated Integral will be
    an instance of NonElementaryIntegral, a subclass of Integral.

    handle_first may be either 'exp' or 'log'.  This changes the order in
    which the extension is built, and may result in a different (but
    equivalent) solution (for an example of this, see issue sympy/sympy#5109).  It is also
    possible that the integral may be computed with one but not the other,
    because not all cases have been implemented yet.  It defaults to 'log' so
    that the outer extension is exponential when possible, because more of the
    exponential case has been implemented.

    If separate_integral is True, the result is returned as a tuple (ans, i),
    where the integral is ans + i, ans is elementary, and i is either a
    NonElementaryIntegral or 0.  This useful if you want to try further
    integrating the NonElementaryIntegral part using other algorithms to
    possibly get a solution in terms of special functions.  It is False by
    default.

    Examples
    ========

    First, we try integrating exp(-x**2). Except for a constant factor of
    2/sqrt(pi), this is the famous error function.

    >>> pprint(risch_integrate(exp(-x**2), x), use_unicode=False)
      /
     |
     |    2
     |  -x
     | E    dx
     |
    /

    The unevaluated Integral in the result means that risch_integrate() has
    proven that exp(-x**2) does not have an elementary anti-derivative.

    In many cases, risch_integrate() can split out the elementary
    anti-derivative part from the nonelementary anti-derivative part.
    For example,

    >>> pprint(risch_integrate((2*log(x)**2 - log(x) - x**2)/(log(x)**3 -
    ...                        x**2*log(x)), x), use_unicode=False)
                                             /
                                            |
      log(-x + log(x))   log(x + log(x))    |   1
    - ---------------- + --------------- +  | ------ dx
             2                  2           | log(x)
                                            |
                                           /

    This means that it has proven that the integral of 1/log(x) is
    nonelementary.  This function is also known as the logarithmic integral,
    and is often denoted as Li(x).

    risch_integrate() currently only accepts purely transcendental functions
    with exponentials and logarithms, though note that this can include
    nested exponentials and logarithms, as well as exponentials with bases
    other than E.

    >>> pprint(risch_integrate(exp(x)*exp(exp(x)), x), use_unicode=False)
     / x\
     \E /
    E
    >>> pprint(risch_integrate(exp(exp(x)), x), use_unicode=False)
      /
     |
     |  / x\
     |  \E /
     | E     dx
     |
    /

    >>> pprint(risch_integrate(x*x**x*log(x) + x**x + x*x**x, x), use_unicode=False)
       x
    x*x
    >>> pprint(risch_integrate(x**x, x), use_unicode=False)
      /
     |
     |  x
     | x  dx
     |
    /

    >>> pprint(risch_integrate(-1/(x*log(x)*log(log(x))**2), x), use_unicode=False)
         1
    -----------
    log(log(x))

    """
    f = sympify(f)

    DE = extension or DifferentialExtension(f, x, handle_first=handle_first, rewrite_complex=rewrite_complex)
    fa, fd = DE.fa, DE.fd

    result = Integer(0)
    for case in reversed(DE.cases):
        if not DE.fa.as_expr().has(DE.t) and not fd.as_expr().has(DE.t) and not case == 'base':
            DE.decrement_level()
            fa, fd = frac_in((fa, fd), DE.t)
            continue

        fa, fd = fa.cancel(fd, include=True)
        if case == 'exp':
            ans, i, b = integrate_hyperexponential(fa, fd, DE, conds=conds)
        elif case == 'primitive':
            ans, i, b = integrate_primitive(fa, fd, DE)
        elif case == 'base':
            # XXX: We can't call ratint() directly here because it doesn't
            # handle polynomials correctly.
            ans = integrate(fa.as_expr()/fd.as_expr(), DE.x, risch=False)
            b = False
            i = Integer(0)
        else:
            raise NotImplementedError('Only exponential and logarithmic '
                                      'extensions are currently supported.')

        result += ans
        if b:
            DE.decrement_level()
            fa, fd = frac_in(i, DE.t)
        else:
            result = result.subs(DE.backsubs)
            if not i.is_zero:
                i = NonElementaryIntegral(i.function.subs(DE.backsubs), i.limits)
            if not separate_integral:
                result += i
                return result
            if isinstance(i, NonElementaryIntegral):
                return result, i
            return result, 0
