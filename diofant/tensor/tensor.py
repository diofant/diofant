"""
This module defines tensors with abstract index notation.

The abstract index notation has been first formalized by Penrose.

Tensor indices are formal objects, with a tensor type; there is no
notion of index range, it is only possible to assign the dimension,
used to trace the Kronecker delta; the dimension can be a Symbol.

The Einstein summation convention is used.
The covariant indices are indicated with a minus sign in front of the index.

For instance the tensor ``t = p(a)*A(b,c)*q(-c)`` has the index ``c``
contracted.

A tensor expression ``t`` can be called; called with its
indices in sorted order it is equal to itself:
in the above example ``t(a, b) == t``;
one can call ``t`` with different indices; ``t(c, d) == p(c)*A(d,a)*q(-a)``.

The contracted indices are dummy indices, internally they have no name,
the indices being represented by a graph-like structure.

Tensors are put in canonical form using ``canon_bp``, which uses
the Butler-Portugal algorithm for canonicalization using the monoterm
symmetries of the tensors.

If there is a (anti)symmetric metric, the indices can be raised and
lowered when the tensor is put in canonical form.
"""

import collections
import functools

from ..combinatorics.tensor_can import (bsgs_direct_product, canonicalize,
                                        get_symmetric_group_sgs, riemann_bsgs)
from ..core import (Add, Basic, Integer, Rational, Symbol, Tuple, symbols,
                    sympify)
from ..core.sympify import CantSympify
from ..external import import_module
from ..matrices import Matrix, eye
from ..utilities.decorator import doctest_depends_on


class TIDS(CantSympify):
    """
    Tensor-index data structure. This contains internal data structures about
    components of a tensor expression, its free and dummy indices.

    To create a ``TIDS`` object via the standard constructor, the required
    arguments are

    WARNING: this class is meant as an internal representation of tensor data
    structures and should not be directly accessed by end users.

    Parameters
    ==========

    components : ``TensorHead`` objects representing the components of the tensor expression.

    free : Free indices in their internal representation.

    dum : Dummy indices in their internal representation.

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> m0, m1, m2, m3 = tensor_indices('m0 m1 m2 m3', Lorentz)
    >>> T = tensorhead('T', [Lorentz]*4, [[1]*4])
    >>> TIDS([T], [(m0, 0, 0), (m3, 3, 0)], [(1, 2, 0, 0)])
    TIDS([T(Lorentz,Lorentz,Lorentz,Lorentz)], [(m0, 0, 0), (m3, 3, 0)], [(1, 2, 0, 0)])

    Notes
    =====

    In short, this has created the components, free and dummy indices for
    the internal representation of a tensor T(m0, m1, -m1, m3).

    Free indices are represented as a list of triplets. The elements of
    each triplet identify a single free index and are

    1. TensorIndex object
    2. position inside the component
    3. component number

    Dummy indices are represented as a list of 4-plets. Each 4-plet stands
    for couple for contracted indices, their original TensorIndex is not
    stored as it is no longer required. The four elements of the 4-plet
    are

    1. position inside the component of the first index.
    2. position inside the component of the second index.
    3. component number of the first index.
    4. component number of the second index.

    """

    def __init__(self, components, free, dum):
        """Initialize self."""
        self.components = components
        self.free = free
        self.dum = dum
        self._ext_rank = len(self.free) + 2*len(self.dum)
        self.dum.sort(key=lambda x: (x[2], x[0]))

    def get_tensors(self):
        """
        Get a list of ``Tensor`` objects having the same ``TIDS`` if multiplied
        by one another.

        """
        indices = self.get_indices()
        components = self.components
        tensors = [None for i in components]  # pre-allocate list
        ind_pos = 0
        for i, component in enumerate(components):
            prev_pos = ind_pos
            ind_pos += component.rank
            tensors[i] = Tensor(component, indices[prev_pos:ind_pos])
        return tensors

    def get_components_with_free_indices(self):
        """
        Get a list of components with their associated indices.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0 m1 m2 m3', Lorentz)
        >>> T = tensorhead('T', [Lorentz]*4, [[1]*4])
        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> t = TIDS.from_components_and_indices([T], [m0, m1, -m1, m3])
        >>> t.get_components_with_free_indices()
        [(T(Lorentz,Lorentz,Lorentz,Lorentz), [(m0, 0, 0), (m3, 3, 0)])]
        >>> t2 = (A(m0)*A(-m0))._tids
        >>> t2.get_components_with_free_indices()
        [(A(Lorentz), []), (A(Lorentz), [])]
        >>> t3 = (A(m0)*A(-m1)*A(-m0)*A(m1))._tids
        >>> t3.get_components_with_free_indices()
        [(A(Lorentz), []), (A(Lorentz), []), (A(Lorentz), []), (A(Lorentz), [])]
        >>> t4 = (A(m0)*A(m1)*A(-m0))._tids
        >>> t4.get_components_with_free_indices()
        [(A(Lorentz), []), (A(Lorentz), [(m1, 0, 1)]), (A(Lorentz), [])]
        >>> t5 = (A(m0)*A(m1)*A(m2))._tids
        >>> t5.get_components_with_free_indices()
        [(A(Lorentz), [(m0, 0, 0)]), (A(Lorentz), [(m1, 0, 1)]), (A(Lorentz), [(m2, 0, 2)])]

        """
        components = self.components
        ret_comp = []

        free_counter = 0
        if len(self.free) == 0:
            return [(comp, []) for comp in components]

        for i, comp in enumerate(components):
            c_free = []
            while free_counter < len(self.free):
                if not self.free[free_counter][2] == i:
                    break

                c_free.append(self.free[free_counter])
                free_counter += 1

                if free_counter >= len(self.free):
                    break
            ret_comp.append((comp, c_free))

        return ret_comp

    @staticmethod
    def from_components_and_indices(components, indices):
        """
        Create a new ``TIDS`` object from ``components`` and ``indices``

        ``components``  ``TensorHead`` objects representing the components
                        of the tensor expression.

        ``indices``     ``TensorIndex`` objects, the indices. Contractions are
                        detected upon construction.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0 m1 m2 m3', Lorentz)
        >>> T = tensorhead('T', [Lorentz]*4, [[1]*4])
        >>> TIDS.from_components_and_indices([T], [m0, m1, -m1, m3])
        TIDS([T(Lorentz,Lorentz,Lorentz,Lorentz)], [(m0, 0, 0), (m3, 3, 0)], [(1, 2, 0, 0)])

        In case of many components the same indices have slightly different
        indexes:

        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> TIDS.from_components_and_indices([A]*4, [m0, m1, -m1, m3])
        TIDS([A(Lorentz), A(Lorentz), A(Lorentz), A(Lorentz)], [(m0, 0, 0), (m3, 0, 3)], [(0, 0, 1, 2)])

        """
        tids = None
        cur_pos = 0
        for i in components:
            tids_sing = TIDS([i], *TIDS.free_dum_from_indices(*indices[cur_pos:cur_pos+i.rank]))
            if tids is None:
                tids = tids_sing
            else:
                tids *= tids_sing
            cur_pos += i.rank

        if tids is None:
            tids = TIDS([], [], [])

        tids.free.sort(key=lambda x: x[0].name)
        tids.dum.sort()

        return tids

    @staticmethod
    def free_dum_from_indices(*indices):
        """
        Convert ``indices`` into ``free``, ``dum`` for single component tensor

        ``free``     list of tuples ``(index, pos, 0)``,
                     where ``pos`` is the position of index in
                     the list of indices formed by the component tensors

        ``dum``      list of tuples ``(pos_contr, pos_cov, 0, 0)``

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0 m1 m2 m3', Lorentz)
        >>> TIDS.free_dum_from_indices(m0, m1, -m1, m3)
        ([(m0, 0, 0), (m3, 3, 0)], [(1, 2, 0, 0)])

        """
        n = len(indices)
        if n == 1:
            return [(indices[0], 0, 0)], []

        # find the positions of the free indices and of the dummy indices
        free = [True]*len(indices)
        index_dict = {}
        dum = []
        for i, index in enumerate(indices):
            name = index._name
            typ = index._tensortype
            contr = index._is_up
            if (name, typ) in index_dict:
                # found a pair of dummy indices
                is_contr, pos = index_dict[(name, typ)]
                # check consistency and update free
                if is_contr:
                    if contr:
                        raise ValueError(f'two equal contravariant indices in slots {pos:d} and {i:d}')
                    free[pos] = False
                    free[i] = False
                else:
                    if contr:
                        free[pos] = False
                        free[i] = False
                    else:
                        raise ValueError(f'two equal covariant indices in slots {pos:d} and {i:d}')
                if contr:
                    dum.append((i, pos, 0, 0))
                else:
                    dum.append((pos, i, 0, 0))
            else:
                index_dict[(name, typ)] = index._is_up, i

        free = [(index, i, 0) for i, index in enumerate(indices) if free[i]]
        free.sort()
        return free, dum

    @staticmethod
    def _check_matrix_indices(f_free, g_free, nc1):
        # This "private" method checks matrix indices.
        # Matrix indices are special as there are only two, and observe
        # anomalous substitution rules to determine contractions.

        dum = []
        # make sure that free indices appear in the same order as in their component:
        f_free.sort(key=lambda x: (x[2], x[1]))
        g_free.sort(key=lambda x: (x[2], x[1]))
        matrix_indices_storage = {}
        transform_right_to_left = {}
        f_pop_pos = []
        g_pop_pos = []
        for free_pos, (ind, i, c) in enumerate(f_free):
            index_type = ind._tensortype
            if ind not in (index_type.auto_left, -index_type.auto_right):
                continue
            matrix_indices_storage[ind] = (free_pos, i, c)

        for free_pos, (ind, i, c) in enumerate(g_free):
            index_type = ind._tensortype
            if ind not in (index_type.auto_left, -index_type.auto_right):
                continue

            if ind == index_type.auto_left:
                if -index_type.auto_right in matrix_indices_storage:
                    other_pos, other_i, other_c = matrix_indices_storage.pop(-index_type.auto_right)
                    dum.append((other_i, i, other_c, c + nc1))
                    # mark to remove other_pos and free_pos from free:
                    g_pop_pos.append(free_pos)
                    f_pop_pos.append(other_pos)
                    continue
                if ind in matrix_indices_storage:
                    other_pos, other_i, other_c = matrix_indices_storage.pop(ind)
                    dum.append((other_i, i, other_c, c + nc1))
                    # mark to remove other_pos and free_pos from free:
                    g_pop_pos.append(free_pos)
                    f_pop_pos.append(other_pos)
                    transform_right_to_left[-index_type.auto_right] = c
                    continue

            if ind in transform_right_to_left:
                other_c = transform_right_to_left.pop(ind)
                if c == other_c:
                    g_free[free_pos] = (index_type.auto_left, i, c)

        for i in sorted(f_pop_pos, reverse=True):
            f_free.pop(i)
        for i in sorted(g_pop_pos, reverse=True):
            g_free.pop(i)
        return dum

    @staticmethod
    def mul(f, g):
        """
        The algorithms performing the multiplication of two ``TIDS`` instances.

        In short, it forms a new ``TIDS`` object, joining components and indices,
        checking that abstract indices are compatible, and possibly contracting
        them.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0 m1 m2 m3', Lorentz)
        >>> T = tensorhead('T', [Lorentz]*4, [[1]*4])
        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> tids_1 = TIDS.from_components_and_indices([T], [m0, m1, -m1, m3])
        >>> tids_2 = TIDS.from_components_and_indices([A], [m2])
        >>> tids_1 * tids_2
        TIDS([T(Lorentz,Lorentz,Lorentz,Lorentz), A(Lorentz)],
             [(m0, 0, 0), (m3, 3, 0), (m2, 0, 1)], [(1, 2, 0, 0)])

        In this case no contraction has been performed.

        >>> tids_3 = TIDS.from_components_and_indices([A], [-m3])
        >>> tids_1 * tids_3
        TIDS([T(Lorentz,Lorentz,Lorentz,Lorentz), A(Lorentz)],
             [(m0, 0, 0)], [(1, 2, 0, 0), (3, 0, 0, 1)])

        Free indices ``m3`` and ``-m3`` are identified as a contracted couple, and are
        therefore transformed into dummy indices.

        A wrong index construction (for example, trying to contract two
        contravariant indices or using indices multiple times) would result in
        an exception:

        >>> tids_4 = TIDS.from_components_and_indices([A], [m3])
        >>> # This raises an exception:
        >>> # tids_1 * tids_4

        """
        def index_up(u):
            return u if u.is_up else -u

        f_free = f.free[:]
        g_free = g.free[:]
        nc1 = len(f.components)
        dum = TIDS._check_matrix_indices(f_free, g_free, nc1)

        # find out which free indices of f and g are contracted
        free_dict1 = {i if i.is_up else -i: (pos, cpos, i) for i, pos, cpos in f_free}
        free_dict2 = {i if i.is_up else -i: (pos, cpos, i) for i, pos, cpos in g_free}
        free_names = set(free_dict1) & set(free_dict2)
        # find the new `free` and `dum`

        dum2 = [(i1, i2, c1 + nc1, c2 + nc1) for i1, i2, c1, c2 in g.dum]
        free1 = [(ind, i, c) for ind, i, c in f_free if index_up(ind) not in free_names]
        free2 = [(ind, i, c + nc1) for ind, i, c in g_free if index_up(ind) not in free_names]
        free = free1 + free2
        dum.extend(f.dum + dum2)
        for name in free_names:
            ipos1, cpos1, ind1 = free_dict1[name]
            ipos2, cpos2, ind2 = free_dict2[name]
            cpos2 += nc1
            if ind1._is_up == ind2._is_up:
                raise ValueError(f'wrong index construction {ind1}')
            if ind1._is_up:
                new_dummy = (ipos1, ipos2, cpos1, cpos2)
            else:
                new_dummy = (ipos2, ipos1, cpos2, cpos1)
            dum.append(new_dummy)
        return f.components + g.components, free, dum

    def __mul__(self, other):
        return TIDS(*self.mul(self, other))

    def __str__(self):
        from ..printing import sstr
        return f'TIDS({sstr(self.components)}, {sstr(self.free)}, {sstr(self.dum)})'

    def sorted_components(self):
        """
        Returns a ``TIDS`` with sorted components

        The sorting is done taking into account the commutation group
        of the component tensors.

        """
        from ..combinatorics.permutations import _af_invert
        cv = list(zip(self.components, range(len(self.components))))
        sign = 1
        n = len(cv) - 1
        for i in range(n):
            for j in range(n, i, -1):
                c = cv[j-1][0].commutes_with(cv[j][0])
                if c not in [0, 1]:
                    continue
                if (cv[j-1][0]._types, cv[j-1][0]._name) > \
                        (cv[j][0]._types, cv[j][0]._name):
                    cv[j-1], cv[j] = cv[j], cv[j-1]
                    if c:
                        sign = -sign

        # perm_inv[new_pos] = old_pos
        components = [x[0] for x in cv]
        perm_inv = [x[1] for x in cv]
        perm = _af_invert(perm_inv)
        free = [(ind, i, perm[c]) for ind, i, c in self.free]
        free.sort()
        dum = [(i1, i2, perm[c1], perm[c2]) for i1, i2, c1, c2 in self.dum]
        dum.sort(key=lambda x: components[x[2]].index_types[x[0]])

        return TIDS(components, free, dum), sign

    def _get_sorted_free_indices_for_canon(self):
        sorted_free = self.free[:]
        sorted_free.sort(key=lambda x: x[0])
        return sorted_free

    def _get_sorted_dum_indices_for_canon(self):
        return sorted(self.dum, key=lambda x: (x[2], x[0]))

    def canon_args(self):
        """
        Returns ``(g, dummies, msym, v)``, the entries of ``canonicalize``

        see ``canonicalize`` in ``tensor_can.py``

        """
        # to be called after sorted_components
        from ..combinatorics.permutations import _af_new
#       types = list(set(self._types))
#       types.sort(key = lambda x: x._name)
        n = self._ext_rank
        g = [None]*n + [n, n+1]
        pos = 0
        vpos = []
        components = self.components
        for t in components:
            vpos.append(pos)
            pos += t._rank
        # ordered indices: first the free indices, ordered by types
        # then the dummy indices, ordered by types and contravariant before
        # covariant
        # g[position in tensor] = position in ordered indices
        for i, (_, ipos, cpos) in enumerate(self._get_sorted_free_indices_for_canon()):
            pos = vpos[cpos] + ipos
            g[pos] = i
        pos = len(self.free)
        j = len(self.free)
        dummies = []
        prev = None
        a = []
        msym = []
        for ipos1, ipos2, cpos1, cpos2 in self._get_sorted_dum_indices_for_canon():
            pos1 = vpos[cpos1] + ipos1
            pos2 = vpos[cpos2] + ipos2
            g[pos1] = j
            g[pos2] = j + 1
            j += 2
            typ = components[cpos1].index_types[ipos1]
            if typ != prev:
                if a:
                    dummies.append(a)
                a = [pos, pos + 1]
                prev = typ
                msym.append(typ.metric_antisym)
            else:
                a.extend([pos, pos + 1])
            pos += 2
        if a:
            dummies.append(a)
        numtyp = []
        prev = None
        for t in components:
            if t == prev:
                numtyp[-1][1] += 1
            else:
                prev = t
                numtyp.append([prev, 1])
        v = []
        for h, n in numtyp:
            if h._comm in (0, 1):
                comm = h._comm
            else:
                comm = TensorManager.get_comm(h._comm, h._comm)
            v.append((h._symmetry.base, h._symmetry.generators, n, comm))
        return _af_new(g), dummies, msym, v

    def perm2tensor(self, g, canon_bp=False):
        """
        Returns a ``TIDS`` instance corresponding to the permutation ``g``

        ``g``  permutation corresponding to the tensor in the representation
        used in canonicalization

        ``canon_bp``   if True, then ``g`` is the permutation
        corresponding to the canonical form of the tensor

        """
        vpos = []
        components = self.components
        pos = 0
        for t in components:
            vpos.append(pos)
            pos += t._rank
        sorted_free = [i[0] for i in self._get_sorted_free_indices_for_canon()]
        nfree = len(sorted_free)
        rank = self._ext_rank
        dum = [[None]*4 for i in range((rank - nfree)//2)]
        free = []
        icomp = -1
        for i in range(rank):
            if i in vpos:
                icomp += vpos.count(i)
                pos0 = i
            ipos = i - pos0
            gi = g[i]
            if gi < nfree:
                ind = sorted_free[gi]
                free.append((ind, ipos, icomp))
            else:
                j = gi - nfree
                idum, cov = divmod(j, 2)
                if cov:
                    dum[idum][1] = ipos
                    dum[idum][3] = icomp
                else:
                    dum[idum][0] = ipos
                    dum[idum][2] = icomp
        dum = [tuple(x) for x in dum]

        return TIDS(components, free, dum)

    def get_indices(self):
        """Get a list of indices, creating new tensor indices to complete dummy indices."""
        components = self.components
        free = self.free
        dum = self.dum
        indices = [None]*self._ext_rank
        start = 0
        pos = 0
        vpos = []
        for t in components:
            vpos.append(pos)
            pos += t.rank
        cdt = collections.defaultdict(int)
        # if the free indices have names with dummy_fmt, start with an
        # index higher than those for the dummy indices
        # to avoid name collisions
        for indx, ipos, cpos in free:
            if indx._name.split('_')[0] == indx._tensortype._dummy_fmt[:-3]:
                cdt[indx._tensortype] = max(cdt[indx._tensortype], int(indx._name.split('_')[1]) + 1)
            start = vpos[cpos]
            indices[start + ipos] = indx
        for ipos1, ipos2, cpos1, cpos2 in dum:
            start1 = vpos[cpos1]
            start2 = vpos[cpos2]
            typ1 = components[cpos1].index_types[ipos1]
            assert typ1 == components[cpos2].index_types[ipos2]
            fmt = typ1._dummy_fmt
            nd = cdt[typ1]
            indices[start1 + ipos1] = TensorIndex(fmt % nd, typ1)
            indices[start2 + ipos2] = TensorIndex(fmt % nd, typ1, False)
            cdt[typ1] += 1
        return indices

    def contract_metric(self, g):
        """
        Returns new TIDS and sign.

        Sign is either 1 or -1, to correct the sign after metric contraction
        (for spinor indices).

        """
        components = self.components
        antisym = g.index_types[0].metric_antisym
        # if not any(x == g for x in components):
        #    return self
        # list of positions of the metric ``g``
        gpos = [i for i, x in enumerate(components) if x == g]
        if not gpos:
            return self, 1
        sign = 1
        dum = self.dum[:]
        free = self.free[:]
        elim = set()
        for gposx in gpos:
            if gposx in elim:
                continue
            free1 = [x for x in free if x[-1] == gposx]
            dum1 = [x for x in dum if gposx in (x[-2], x[-1])]
            if not dum1:
                continue
            elim.add(gposx)
            if len(dum1) == 2:
                if not antisym:
                    dum10, dum11 = dum1
                    if dum10[3] == gposx:
                        # the index with pos p0 and component c0 is contravariant
                        c0 = dum10[2]
                        p0 = dum10[0]
                    else:
                        # the index with pos p0 and component c0 is covariant
                        c0 = dum10[3]
                        p0 = dum10[1]
                    if dum11[3] == gposx:
                        # the index with pos p1 and component c1 is contravariant
                        c1 = dum11[2]
                        p1 = dum11[0]
                    else:
                        # the index with pos p1 and component c1 is covariant
                        c1 = dum11[3]
                        p1 = dum11[1]
                    dum.append((p0, p1, c0, c1))
                else:
                    dum10, dum11 = dum1
                    # change the sign to bring the indices of the metric to contravariant
                    # form; change the sign if dum10 has the metric index in position 0
                    if dum10[3] == gposx:
                        # the index with pos p0 and component c0 is contravariant
                        c0 = dum10[2]
                        p0 = dum10[0]
                        if dum10[1] == 1:
                            sign = -sign
                    else:
                        # the index with pos p0 and component c0 is covariant
                        c0 = dum10[3]
                        p0 = dum10[1]
                        if dum10[0] == 0:
                            sign = -sign
                    if dum11[3] == gposx:
                        # the index with pos p1 and component c1 is contravariant
                        c1 = dum11[2]
                        p1 = dum11[0]
                        sign = -sign
                    else:
                        # the index with pos p1 and component c1 is covariant
                        c1 = dum11[3]
                        p1 = dum11[1]
                    dum.append((p0, p1, c0, c1))

            elif len(dum1) == 1:
                if not antisym:
                    dp0, dp1, dc0, dc1 = dum1[0]
                    if dc0 == dc1:
                        # g(i, -i)
                        typ = g.index_types[0]
                        if typ._dim is None:
                            raise ValueError('dimension not assigned')
                        sign = sign*typ._dim

                    else:
                        # g(i0, i1)*p(-i1)
                        if dc0 == gposx:
                            p1 = dp1
                            c1 = dc1
                        else:
                            p1 = dp0
                            c1 = dc0
                        ind, _, c = free1[0]
                        free.append((ind, p1, c1))
                else:
                    dp0, dp1, dc0, dc1 = dum1[0]
                    if dc0 == dc1:
                        # g(i, -i)
                        typ = g.index_types[0]
                        if typ._dim is None:
                            raise ValueError('dimension not assigned')
                        sign = sign*typ._dim

                        if dp0 < dp1:
                            # g(i, -i) = -D with antisymmetric metric
                            sign = -sign
                    else:
                        # g(i0, i1)*p(-i1)
                        if dc0 == gposx:
                            p1 = dp1
                            c1 = dc1
                            if dp0 == 0:
                                sign = -sign
                        else:
                            p1 = dp0
                            c1 = dc0
                        ind, _, c = free1[0]
                        free.append((ind, p1, c1))
            dum = [x for x in dum if x not in dum1]
            free = [x for x in free if x not in free1]

        shift = 0
        shifts = [0]*len(components)
        for i in range(len(components)):
            if i in elim:
                shift += 1
                continue
            shifts[i] = shift
        free = [(ind, p, c - shifts[c]) for (ind, p, c) in free if c not in elim]
        dum = [(p0, p1, c0 - shifts[c0], c1 - shifts[c1]) for i, (p0, p1, c0, c1) in enumerate(dum) if c0 not in elim and c1 not in elim]
        components = [c for i, c in enumerate(components) if i not in elim]
        tids = TIDS(components, free, dum)
        return tids, sign


class _TensorDataLazyEvaluator(CantSympify):
    """
    EXPERIMENTAL: do not rely on this class, it may change without deprecation
    warnings in future versions of Diofant.

    This object contains the logic to associate components data to a tensor
    expression. Components data are set via the ``.data`` property of tensor
    expressions, is stored inside this class as a mapping between the tensor
    expression and the ``ndarray``.

    Computations are executed lazily: whereas the tensor expressions can have
    contractions, tensor products, and additions, components data are not
    computed until they are accessed by reading the ``.data`` property
    associated to the tensor expression.

    """

    _substitutions_dict: dict[object, object] = {}
    _substitutions_dict_tensmul: dict[object, object] = {}

    def __getitem__(self, key):
        dat = self._get(key)
        if dat is None:
            return

        numpy = import_module('numpy')
        if not isinstance(dat, numpy.ndarray):
            return dat

        if dat.ndim == 0:
            return dat[()]
        if dat.ndim == 1 and dat.size == 1:
            return dat[0]
        return dat

    def _get(self, key):
        """
        Retrieve ``data`` associated with ``key``.

        This algorithm looks into ``self._substitutions_dict`` for all
        ``TensorHead`` in the ``TensExpr`` (or just ``TensorHead`` if key is a
        TensorHead instance). It reconstructs the components data that the
        tensor expression should have by performing on components data the
        operations that correspond to the abstract tensor operations applied.

        Metric tensor is handled in a different manner: it is pre-computed in
        ``self._substitutions_dict_tensmul``.

        """
        if key in self._substitutions_dict:
            return self._substitutions_dict[key]

        if isinstance(key, TensorHead):
            return

        if isinstance(key, Tensor):
            # special case to handle metrics. Metric tensors cannot be
            # constructed through contraction by the metric, their
            # components show if they are a matrix or its inverse.
            signature = tuple(i.is_up for i in key.get_indices())
            srch = (key.component,) + signature
            if srch in self._substitutions_dict_tensmul:
                return self._substitutions_dict_tensmul[srch]
            return self.data_tensmul_from_tensorhead(key, key.component)

        if isinstance(key, TensMul):
            tensmul_list = key.split()
            if len(tensmul_list) == 1 and len(tensmul_list[0].components) == 1:
                # special case to handle metrics. Metric tensors cannot be
                # constructed through contraction by the metric, their
                # components show if they are a matrix or its inverse.
                signature = tuple(i.is_up for i in tensmul_list[0].get_indices())
                srch = (tensmul_list[0].components[0],) + signature
                if srch in self._substitutions_dict_tensmul:
                    return self._substitutions_dict_tensmul[srch]
            data_list = [self.data_tensmul_from_tensorhead(i, i.components[0]) for i in tensmul_list]
            if all(i is None for i in data_list):
                return
            if any(i is None for i in data_list):
                raise ValueError('Mixing tensors with associated components '
                                 'data with tensors without components data')
            data_result, _ = self.data_product_tensors(data_list, tensmul_list)
            return data_result

        if isinstance(key, TensAdd):
            sumvar = Integer(0)
            data_list = []
            free_args_list = []
            for arg in key.args:
                if isinstance(arg, TensExpr):
                    data_list.append(arg.data)
                    free_args_list.append([x[0] for x in arg.free])
                else:
                    data_list.append(arg)
                    free_args_list.append([])
            if all(i is None for i in data_list):
                return
            if any(i is None for i in data_list):
                raise ValueError('Mixing tensors with associated components '
                                 'data with tensors without components data')

            numpy = import_module('numpy')
            for data, free_args in zip(data_list, free_args_list):
                if len(free_args) < 2:
                    sumvar += data
                else:
                    free_args_pos = {y: x for x, y in enumerate(free_args)}
                    axes = [free_args_pos[arg] for arg in key.free_args]
                    sumvar += numpy.transpose(data, axes)
            return sumvar

    def data_tensorhead_from_tensmul(self, data, tensmul, tensorhead):
        """
        This method is used when assigning components data to a ``TensMul``
        object, it converts components data to a fully contravariant ndarray,
        which is then stored according to the ``TensorHead`` key.

        """
        if data is not None:
            return self._correct_signature_from_indices(
                data,
                tensmul.get_indices(),
                tensmul.free,
                tensmul.dum,
                True)

    def data_tensmul_from_tensorhead(self, tensmul, tensorhead):
        """
        This method corrects the components data to the right signature
        (covariant/contravariant) using the metric associated with each
        ``TensorIndexType``.

        """
        if tensorhead.data is not None:
            return self._correct_signature_from_indices(
                tensorhead.data,
                tensmul.get_indices(),
                tensmul.free,
                tensmul.dum)

    def data_product_tensors(self, data_list, tensmul_list):
        """
        Given a ``data_list``, list of ``ndarray``'s and a ``tensmul_list``,
        list of ``TensMul`` instances, compute the resulting ``ndarray``,
        after tensor products and contractions.

        """
        def data_mul(f, g):
            """
            Multiplies two ``ndarray`` objects, it first calls ``TIDS.mul``,
            then checks which indices have been contracted, and finally
            contraction operation on data, according to the contracted indices.

            """
            data1, tensmul1 = f
            data2, tensmul2 = g
            components, free, dum = TIDS.mul(tensmul1, tensmul2)
            data = _TensorDataLazyEvaluator._contract_ndarray(tensmul1.free, tensmul2.free, data1, data2)
            # TODO: do this more efficiently... maybe by just passing an index list
            # to .data_product_tensor(...)
            return data, TensMul.from_TIDS(Integer(1), TIDS(components, free, dum))

        return functools.reduce(data_mul, zip(data_list, tensmul_list))

    def _assign_data_to_tensor_expr(self, key, data):
        if isinstance(key, TensAdd):
            raise ValueError('cannot assign data to TensAdd')
        # here it is assumed that `key` is a `TensMul` instance.
        if len(key.components) != 1:
            raise ValueError('cannot assign data to TensMul with multiple components')
        tensorhead = key.components[0]
        newdata = self.data_tensorhead_from_tensmul(data, key, tensorhead)
        return tensorhead, newdata

    def _check_permutations_on_data(self, tens, data):
        import numpy

        if isinstance(tens, TensorHead):
            rank = tens.rank
            generators = tens.symmetry.generators
        elif isinstance(tens, Tensor):
            rank = tens.rank
            generators = tens.components[0].symmetry.generators
        elif isinstance(tens, TensorIndexType):
            rank = tens.metric.rank
            generators = tens.metric.symmetry.generators

        # Every generator is a permutation, check that by permuting the array
        # by that permutation, the array will be the same, except for a
        # possible sign change if the permutation admits it.
        for gener in generators:
            sign_change = +1 if (gener(rank) == rank) else -1
            data_swapped = data
            last_data = data
            permute_axes = list(map(gener, range(rank)))
            # the order of a permutation is the number of times to get the
            # identity by applying that permutation.
            for _ in range(gener.order()-1):
                data_swapped = numpy.transpose(data_swapped, permute_axes)
                # if any value in the difference array is non-zero, raise an error:
                if (last_data - sign_change*data_swapped).any():
                    raise ValueError('Component data symmetry structure error')
                last_data = data_swapped

    def __setitem__(self, key, value):
        """
        Set the components data of a tensor object/expression.

        Components data are transformed to the all-contravariant form and stored
        with the corresponding ``TensorHead`` object. If a ``TensorHead`` object
        cannot be uniquely identified, it will raise an error.

        """
        data = _TensorDataLazyEvaluator.parse_data(value)
        self._check_permutations_on_data(key, data)

        # TensorHead and TensorIndexType can be assigned data directly, while
        # TensMul must first convert data to a fully contravariant form, and
        # assign it to its corresponding TensorHead single component.
        if not isinstance(key, (TensorHead, TensorIndexType)):
            key, data = self._assign_data_to_tensor_expr(key, data)

        if isinstance(key, TensorHead):
            for dim, indextype in zip(data.shape, key.index_types):
                if indextype.data is None:
                    raise ValueError(f'index type {indextype} has no components data'
                                     ' associated (needed to raise/lower index')
                if indextype.dim is None:
                    continue
                if dim != indextype.dim:
                    raise ValueError('wrong dimension of ndarray')
        self._substitutions_dict[key] = data

    def __delitem__(self, key):
        del self._substitutions_dict[key]

    def __contains__(self, key):
        return key in self._substitutions_dict

    @staticmethod
    def _contract_ndarray(free1, free2, ndarray1, ndarray2):
        numpy = import_module('numpy')

        def ikey(x):
            return x[2], x[1]

        free1 = free1[:]
        free2 = free2[:]
        free1.sort(key=ikey)
        free2.sort(key=ikey)
        self_free = [_[0] for _ in free1]
        axes1 = []
        axes2 = []
        for jpos, jindex in enumerate(free2):
            if -jindex[0] in self_free:
                nidx = self_free.index(-jindex[0])
            else:
                continue
            axes1.append(nidx)
            axes2.append(jpos)

        contracted_ndarray = numpy.tensordot(
            ndarray1,
            ndarray2,
            (axes1, axes2)
        )
        return contracted_ndarray

    def add_metric_data(self, metric, data):
        """
        Assign data to the ``metric`` tensor. The metric tensor behaves in an
        anomalous way when raising and lowering indices.

        A fully covariant metric is the inverse transpose of the fully
        contravariant metric (it is meant matrix inverse). If the metric is
        symmetric, the transpose is not necessary and mixed
        covariant/contravariant metrics are Kronecker deltas.

        """
        # hard assignment, data should not be added to `TensorHead` for metric:
        # the problem with `TensorHead` is that the metric is anomalous, i.e.
        # raising and lowering the index means considering the metric or its
        # inverse, this is not the case for other tensors.
        self._substitutions_dict_tensmul[metric, True, True] = data
        inverse_transpose = self.inverse_transpose_matrix(data)
        # in symmetric spaces, the traspose is the same as the original matrix,
        # the full covariant metric tensor is the inverse transpose, so this
        # code will be able to handle non-symmetric metrics.
        self._substitutions_dict_tensmul[metric, False, False] = inverse_transpose
        # now mixed cases, these are identical to the unit matrix if the metric
        # is symmetric.
        m = Matrix(data)
        invt = Matrix(inverse_transpose)
        self._substitutions_dict_tensmul[metric, True, False] = m * invt
        self._substitutions_dict_tensmul[metric, False, True] = invt * m

    @staticmethod
    def _flip_index_by_metric(data, metric, pos):
        numpy = import_module('numpy')

        data = numpy.tensordot(
            metric,
            data,
            (1, pos))
        return numpy.rollaxis(data, 0, pos+1)

    @staticmethod
    def inverse_matrix(ndarray):
        m = Matrix(ndarray).inv()
        return _TensorDataLazyEvaluator.parse_data(m)

    @staticmethod
    def inverse_transpose_matrix(ndarray):
        m = Matrix(ndarray).inv().T
        return _TensorDataLazyEvaluator.parse_data(m)

    @staticmethod
    def _correct_signature_from_indices(data, indices, free, dum, inverse=False):
        """
        Utility function to correct the values inside the components data
        ndarray according to whether indices are covariant or contravariant.

        It uses the metric matrix to lower values of covariant indices.

        """
        numpy = import_module('numpy')
        # change the ndarray values according covariantness/contravariantness of the indices
        # use the metric
        for i, indx in enumerate(indices):
            if not indx.is_up and not inverse:
                data = _TensorDataLazyEvaluator._flip_index_by_metric(data, indx._tensortype.data, i)
            elif not indx.is_up and inverse:
                data = _TensorDataLazyEvaluator._flip_index_by_metric(
                    data,
                    _TensorDataLazyEvaluator.inverse_matrix(indx._tensortype.data),
                    i
                )

        if len(dum) > 0:
            # perform contractions
            axes1 = []
            axes2 = []
            for i, indx1 in enumerate(indices):
                try:
                    nd = indices[:i].index(-indx1)
                except ValueError:
                    continue
                axes1.append(nd)
                axes2.append(i)

            for ax1, ax2 in zip(axes1, axes2):
                data = numpy.trace(data, axis1=ax1, axis2=ax2)
        return data

    @staticmethod
    @doctest_depends_on(modules=('numpy',))
    def parse_data(data):
        """
        Transform ``data`` to a numpy ndarray. The parameter ``data`` may
        contain data in various formats, e.g. nested lists, diofant ``Matrix``,
        and so on.

        Examples
        ========

        >>> print(str(_TensorDataLazyEvaluator.parse_data([1, 3, -6, 12])))
        [1 3 -6 12]

        >>> print(str(_TensorDataLazyEvaluator.parse_data([[1, 2], [4, 7]])))
        [[1 2]
         [4 7]]

        """
        numpy = import_module('numpy')

        if (numpy is not None) and (not isinstance(data, numpy.ndarray)):
            vsympify = numpy.vectorize(sympify)
            data = vsympify(numpy.array(data))
        return data


_tensor_data_substitution_dict = _TensorDataLazyEvaluator()


class _TensorManager:
    """
    Class to manage tensor properties.

    Notes
    =====

    Tensors belong to tensor commutation groups; each group has a label
    ``comm``; there are predefined labels:

    ``0``   tensors commuting with any other tensor

    ``1``   tensors anticommuting among themselves

    ``2``   tensors not commuting, apart with those with ``comm=0``

    Other groups can be defined using ``set_comm``; tensors in those
    groups commute with those with ``comm=0``; by default they
    do not commute with any other group.

    """

    def __init__(self):
        self._comm_init()

    def _comm_init(self):
        self._comm = [{} for i in range(3)]
        for i in range(3):
            self._comm[0][i] = 0
            self._comm[i][0] = 0
        self._comm[1][1] = 1
        self._comm[2][1] = None
        self._comm[1][2] = None
        self._comm_symbols2i = {0: 0, 1: 1, 2: 2}
        self._comm_i2symbol = {0: 0, 1: 1, 2: 2}

    @property
    def comm(self):
        return self._comm

    def comm_symbols2i(self, i):
        """
        Get the commutation group number corresponding to ``i``

        ``i`` can be a symbol or a number or a string

        If ``i`` is not already defined its commutation group number
        is set.

        """
        if i not in self._comm_symbols2i:
            n = len(self._comm)
            self._comm.append({})
            self._comm[n][0] = 0
            self._comm[0][n] = 0
            self._comm_symbols2i[i] = n
            self._comm_i2symbol[n] = i
            return n
        return self._comm_symbols2i[i]

    def comm_i2symbol(self, i):
        """Returns the symbol corresponding to the commutation group number."""
        return self._comm_i2symbol[i]

    def set_comm(self, i, j, c):
        """
        Set the commutation parameter ``c`` for commutation groups ``i, j``

        Parameters
        ==========

        i, j : symbols representing commutation groups

        c  :  group commutation number

        Notes
        =====

        ``i, j`` can be symbols, strings or numbers,
        apart from ``0, 1`` and ``2`` which are reserved respectively
        for commuting, anticommuting tensors and tensors not commuting
        with any other group apart with the commuting tensors.
        For the remaining cases, use this method to set the commutation rules;
        by default ``c=None``.

        The group commutation number ``c`` is assigned in correspondence
        to the group commutation symbols; it can be

        0        commuting

        1        anticommuting

        None     no commutation property

        Examples
        ========

        ``G`` and ``GH`` do not commute with themselves and commute with
        each other; A is commuting.

        >>> Lorentz = TensorIndexType('Lorentz')
        >>> i0, i1, i2, i3, i4 = tensor_indices('i0:5', Lorentz)
        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> G = tensorhead('G', [Lorentz], [[1]], 'Gcomm')
        >>> GH = tensorhead('GH', [Lorentz], [[1]], 'GHcomm')
        >>> TensorManager.set_comm('Gcomm', 'GHcomm', 0)
        >>> (GH(i1)*G(i0)).canon_bp()
        G(i0)*GH(i1)
        >>> (G(i1)*G(i0)).canon_bp()
        G(i1)*G(i0)
        >>> (G(i1)*A(i0)).canon_bp()
        A(i0)*G(i1)

        """
        if c not in (0, 1, None):
            raise ValueError('`c` can assume only the values 0, 1 or None')

        if i not in self._comm_symbols2i:
            n = len(self._comm)
            self._comm.append({})
            self._comm[n][0] = 0
            self._comm[0][n] = 0
            self._comm_symbols2i[i] = n
            self._comm_i2symbol[n] = i
        if j not in self._comm_symbols2i:
            n = len(self._comm)
            self._comm.append({})
            self._comm[0][n] = 0
            self._comm[n][0] = 0
            self._comm_symbols2i[j] = n
            self._comm_i2symbol[n] = j
        ni = self._comm_symbols2i[i]
        nj = self._comm_symbols2i[j]
        self._comm[ni][nj] = c
        self._comm[nj][ni] = c

    def set_comms(self, *args):
        """
        Set the commutation group numbers ``c`` for symbols ``i, j``

        Parameters
        ==========

        args : sequence of ``(i, j, c)``

        """
        for i, j, c in args:
            self.set_comm(i, j, c)

    def get_comm(self, i, j):
        """
        Return the commutation parameter for commutation group numbers ``i, j``

        see ``_TensorManager.set_comm``

        """
        return self._comm[i].get(j, 0 if i == 0 or j == 0 else None)

    def clear(self):
        """Clear the TensorManager."""
        self._comm_init()


TensorManager = _TensorManager()


@doctest_depends_on(modules=('numpy',))
class TensorIndexType(Basic):
    """
    A TensorIndexType is characterized by its name and its metric.

    Parameters
    ==========

    name : name of the tensor type

    metric : metric symmetry or metric object or ``None``


    dim : dimension, it can be a symbol or an integer or ``None``

    eps_dim : dimension of the epsilon tensor

    dummy_fmt : name of the head of dummy indices

    Attributes
    ==========

    ``name``
    ``metric_name`` : str
        it is 'metric' or metric.name
    ``metric_antisym``
    ``metric`` : TensorType
        the metric tensor
    ``delta`` : ``Kronecker delta``
    ``epsilon`` : the ``Levi-Civita epsilon`` tensor
    ``dim``
    ``dim_eps``
    ``dummy_fmt``
    ``data`` : a property to add ``ndarray`` values, to work in a specified basis.

    Notes
    =====

    The ``metric`` parameter can be:
    ``metric = False`` symmetric metric (in Riemannian geometry)

    ``metric = True`` antisymmetric metric (for spinor calculus)

    ``metric = None``  there is no metric

    ``metric`` can be an object having ``name`` and ``antisym`` attributes.


    If there is a metric the metric is used to raise and lower indices.

    In the case of antisymmetric metric, the following raising and
    lowering conventions will be adopted:

    ``psi(a) = g(a, b)*psi(-b); chi(-a) = chi(b)*g(-b, -a)``

    ``g(-a, b) = delta(-a, b); g(b, -a) = -delta(a, -b)``

    where ``delta(-a, b) = delta(b, -a)`` is the ``Kronecker delta``
    (see ``TensorIndex`` for the conventions on indices).

    If there is no metric it is not possible to raise or lower indices;
    e.g. the index of the defining representation of ``SU(N)``
    is 'covariant' and the conjugate representation is
    'contravariant'; for ``N > 2`` they are linearly independent.

    ``eps_dim`` is by default equal to ``dim``, if the latter is an integer;
    else it can be assigned (for use in naive dimensional regularization);
    if ``eps_dim`` is not an integer ``epsilon`` is ``None``.

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> Lorentz.metric
    metric(Lorentz,Lorentz)

    Examples with metric components data added, this means it is working on a
    fixed basis:

    >>> Lorentz.data = [1, -1, -1, -1]
    >>> print(sstr(Lorentz))
    TensorIndexType(Lorentz, 0)
    >>> print(str(Lorentz.data))
    [[1 0 0 0]
    [0 -1 0 0]
    [0 0 -1 0]
    [0 0 0 -1]]

    """

    def __new__(cls, name, metric=False, dim=None, eps_dim=None,
                dummy_fmt=None):

        if isinstance(name, str):
            name = Symbol(name)
        obj = Basic.__new__(cls, name, Integer(1) if metric else Integer(0))
        obj._name = str(name)
        if not dummy_fmt:
            obj._dummy_fmt = f'{obj.name}_%d'
        else:
            obj._dummy_fmt = f'{dummy_fmt}_%d'
        if metric is None:
            obj.metric_antisym = None
            obj.metric = None
        else:
            if metric in (True, False, 0, 1):
                metric_name = 'metric'
                obj.metric_antisym = metric
            else:
                metric_name = metric.name
                obj.metric_antisym = metric.antisym
            sym2 = TensorSymmetry(get_symmetric_group_sgs(2, obj.metric_antisym))
            S2 = TensorType([obj]*2, sym2)
            obj.metric = S2(metric_name)
            obj.metric._matrix_behavior = True

        obj._dim = dim
        obj._delta = obj.get_kronecker_delta()
        obj._eps_dim = eps_dim if eps_dim else dim
        obj._epsilon = obj.get_epsilon()
        obj._autogenerated = []
        obj._auto_left = None
        obj._auto_right = None
        return obj

    @property
    def auto_right(self):
        if self._auto_right is None:
            self._auto_right = TensorIndex('auto_right', self)
        return self._auto_right

    @property
    def auto_left(self):
        if self._auto_left is None:
            self._auto_left = TensorIndex('auto_left', self)
        return self._auto_left

    @property
    def data(self):
        return _tensor_data_substitution_dict[self]

    @data.setter
    def data(self, data):
        # This assignment is a bit controversial, should metric components be assigned
        # to the metric only or also to the TensorIndexType object? The advantage here
        # is the ability to assign a 1D array and transform it to a 2D diagonal array.
        numpy = import_module('numpy')
        data = _TensorDataLazyEvaluator.parse_data(data)
        if data.ndim > 2:
            raise ValueError('data have to be of rank 1 (diagonal metric) or 2.')
        if data.ndim == 1:
            if self.dim is not None:
                nda_dim = data.shape[0]
                if nda_dim != self.dim:
                    raise ValueError('Dimension mismatch')

            dim = data.shape[0]
            newndarray = numpy.zeros((dim, dim), dtype=object)
            for i, val in enumerate(data):
                newndarray[i, i] = val
            data = newndarray
        dim1, dim2 = data.shape
        if dim1 != dim2:
            raise ValueError('Non-square matrix tensor.')
        if self.dim is not None:
            if self.dim != dim1:
                raise ValueError('Dimension mismatch')
        _tensor_data_substitution_dict[self] = data
        _tensor_data_substitution_dict.add_metric_data(self.metric, data)
        delta = self.get_kronecker_delta()
        i1 = TensorIndex('i1', self)
        i2 = TensorIndex('i2', self)
        delta(i1, -i2).data = _TensorDataLazyEvaluator.parse_data(eye(dim1))

    @data.deleter
    def data(self):
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]

    @property
    def name(self):
        return self._name

    @property
    def dim(self):
        return self._dim

    @property
    def delta(self):
        return self._delta

    @property
    def eps_dim(self):
        return self._eps_dim

    @property
    def epsilon(self):
        return self._epsilon

    def get_kronecker_delta(self):
        sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
        S2 = TensorType([self]*2, sym2)
        delta = S2('KD')
        delta._matrix_behavior = True
        return delta

    def get_epsilon(self):
        if not isinstance(self._eps_dim, int):
            return
        sym = TensorSymmetry(get_symmetric_group_sgs(self._eps_dim, 1))
        Sdim = TensorType([self]*self._eps_dim, sym)
        epsilon = Sdim('Eps')
        return epsilon

    def __lt__(self, other):
        return self.name < other.name

    def __str__(self):
        return self.name

    __repr__ = __str__


@doctest_depends_on(modules=('numpy',))
class TensorIndex(Basic):
    """
    Represents an abstract tensor index.

    Parameters
    ==========

    name : name of the index, or ``True`` if you want it to be automatically assigned
    tensortype : ``TensorIndexType`` of the index
    is_up :  flag for contravariant index

    Attributes
    ==========

    ``name``
    ``tensortype``
    ``is_up``

    Notes
    =====

    Tensor indices are contracted with the Einstein summation convention.

    An index can be in contravariant or in covariant form; in the latter
    case it is represented prepending a ``-`` to the index name.

    Dummy indices have a name with head given by ``tensortype._dummy_fmt``

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> i = TensorIndex('i', Lorentz)
    >>> i
    i
    >>> sym1 = TensorSymmetry(*get_symmetric_group_sgs(1))
    >>> S1 = TensorType([Lorentz], sym1)
    >>> A, B = S1('A B')
    >>> A(i)*B(-i)
    A(L_0)*B(-L_0)

    If you want the index name to be automatically assigned, just put ``True``
    in the ``name`` field, it will be generated using the reserved character
    ``_`` in front of its name, in order to avoid conflicts with possible
    existing indices:

    >>> i0 = TensorIndex(True, Lorentz)
    >>> i0
    _i0
    >>> i1 = TensorIndex(True, Lorentz)
    >>> i1
    _i1
    >>> A(i0)*B(-i1)
    A(_i0)*B(-_i1)
    >>> A(i0)*B(-i0)
    A(L_0)*B(-L_0)

    """

    def __new__(cls, name, tensortype, is_up=True):
        if isinstance(name, str):
            name_symbol = Symbol(name)
        elif isinstance(name, Symbol):
            name_symbol = name
        elif name is True:
            name = f'_i{len(tensortype._autogenerated)}'
            name_symbol = Symbol(name)
            tensortype._autogenerated.append(name_symbol)
        else:
            raise ValueError('invalid name')

        obj = Basic.__new__(cls, name_symbol, tensortype, Integer(1) if is_up else Integer(0))
        obj._name = str(name)
        obj._tensortype = tensortype
        obj._is_up = is_up
        return obj

    @property
    def name(self):
        return self._name

    @property
    def tensortype(self):
        return self._tensortype

    @property
    def is_up(self):
        return self._is_up

    def _print(self):
        s = self._name
        if not self._is_up:
            s = f'-{s}'
        return s

    def __lt__(self, other):
        return (self._tensortype, self._name) < (other._tensortype, other._name)

    def __neg__(self):
        t1 = TensorIndex(self._name, self._tensortype,
                         (not self._is_up))
        return t1


def tensor_indices(s, typ):
    """
    Returns list of tensor indices given their names and their types

    Parameters
    ==========

    s : string of comma separated names of indices

    typ : list of ``TensorIndexType`` of the indices

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b, c, d = tensor_indices('a b c d', Lorentz)

    """
    if isinstance(s, str):
        a = [x.name for x in symbols(s, seq=True)]
    else:
        raise ValueError('expecting a string')

    tilist = [TensorIndex(i, typ) for i in a]
    if len(tilist) == 1:
        return tilist[0]
    return tilist


@doctest_depends_on(modules=('numpy',))
class TensorSymmetry(Basic):
    """
    Monoterm symmetry of a tensor

    Parameters
    ==========

    bsgs : tuple ``(base, sgs)`` BSGS of the symmetry of the tensor

    Attributes
    ==========

    ``base`` : Tuple
        base of the BSGS
    ``generators`` : Tuple
        generators of the BSGS
    ``rank`` : Tuple
        rank of the tensor

    Notes
    =====

    A tensor can have an arbitrary monoterm symmetry provided by its BSGS.
    Multiterm symmetries, like the cyclic symmetry of the Riemann tensor,
    are not covered.

    See Also
    ========

    diofant.combinatorics.tensor_can.get_symmetric_group_sgs

    Examples
    ========

    Define a symmetric tensor

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')

    """

    def __new__(cls, *args, **kw_args):
        if len(args) == 1:
            base, generators = args[0]
        elif len(args) == 2:
            base, generators = args
        else:
            raise TypeError('bsgs required, either two separate parameters or one tuple')

        if not isinstance(base, Tuple):
            base = Tuple(*base)
        if not isinstance(generators, Tuple):
            generators = Tuple(*generators)
        obj = Basic.__new__(cls, base, generators, **kw_args)
        return obj

    @property
    def base(self):
        return self.args[0]

    @property
    def generators(self):
        return self.args[1]

    @property
    def rank(self):
        return self.args[1][0].size - 2


def tensorsymmetry(*args):
    """
    Return a ``TensorSymmetry`` object.

    One can represent a tensor with any monoterm slot symmetry group
    using a BSGS.

    ``args`` can be a BSGS
    ``args[0]``    base
    ``args[1]``    sgs

    Usually tensors are in (direct products of) representations
    of the symmetric group;
    ``args`` can be a list of lists representing the shapes of Young tableaux

    Notes
    =====

    For instance:
    ``[[1]]``       vector
    ``[[1]*n]``     symmetric tensor of rank ``n``
    ``[[n]]``       antisymmetric tensor of rank ``n``
    ``[[2, 2]]``    monoterm slot symmetry of the Riemann tensor
    ``[[1],[1]]``   vector*vector
    ``[[2],[1],[1]`` (antisymmetric tensor)*vector*vector

    Notice that with the shape ``[2, 2]`` we associate only the monoterm
    symmetries of the Riemann tensor; this is an abuse of notation,
    since the shape ``[2, 2]`` corresponds usually to the irreducible
    representation characterized by the monoterm symmetries and by the
    cyclic symmetry.

    Examples
    ========

    Symmetric tensor using a Young tableau

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = tensorsymmetry([1, 1])
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')

    Symmetric tensor using a ``BSGS`` (base, strong generator set)

    >>> sym2 = tensorsymmetry(*get_symmetric_group_sgs(2))
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')

    """
    from ..combinatorics import Permutation

    def tableau2bsgs(a):
        if len(a) == 1:
            # antisymmetric vector
            n = a[0]
            bsgs = get_symmetric_group_sgs(n, 1)
        else:
            if all(x == 1 for x in a):
                # symmetric vector
                n = len(a)
                bsgs = get_symmetric_group_sgs(n)
            elif a == [2, 2]:
                bsgs = riemann_bsgs
            else:
                raise NotImplementedError
        return bsgs

    if not args:
        return TensorSymmetry(Tuple(), Tuple(Permutation(1)))

    if len(args) == 2 and isinstance(args[1][0], Permutation):
        return TensorSymmetry(args)
    base, sgs = tableau2bsgs(args[0])
    for a in args[1:]:
        basex, sgsx = tableau2bsgs(a)
        base, sgs = bsgs_direct_product(base, sgs, basex, sgsx)
    return TensorSymmetry(Tuple(base, sgs))


@doctest_depends_on(modules=('numpy',))
class TensorType(Basic):
    """
    Class of tensor types.

    Parameters
    ==========

    index_types : list of ``TensorIndexType`` of the tensor indices
    symmetry : ``TensorSymmetry`` of the tensor

    Attributes
    ==========

    ``index_types``
    ``symmetry``
    ``types`` : list of ``TensorIndexType`` without repetitions

    Examples
    ========

    Define a symmetric tensor

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = tensorsymmetry([1, 1])
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')

    """

    is_commutative = False

    def __new__(cls, index_types, symmetry, **kw_args):
        assert symmetry.rank == len(index_types)
        obj = Basic.__new__(cls, Tuple(*index_types), symmetry, **kw_args)
        return obj

    @property
    def index_types(self):
        return self.args[0]

    @property
    def symmetry(self):
        return self.args[1]

    @property
    def types(self):
        return sorted(set(self.index_types), key=lambda x: x.name)

    def __str__(self):
        return f'TensorType({[str(x) for x in self.index_types]})'

    def __call__(self, s, comm=0, matrix_behavior=0):
        """
        Return a TensorHead object or a list of TensorHead objects.

        ``s``  name or string of names

        ``comm``: commutation group number
        see ``_TensorManager.set_comm``

        Examples
        ========

        Define symmetric tensors ``V``, ``W`` and ``G``, respectively
        commuting, anticommuting and with no commutation symmetry

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b = tensor_indices('a b', Lorentz)
        >>> sym2 = tensorsymmetry([1]*2)
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> V = S2('V')
        >>> W = S2('W', 1)
        >>> G = S2('G', 2)
        >>> canon_bp(V(a, b)*V(-b, -a))
        V(L_0, L_1)*V(-L_0, -L_1)
        >>> canon_bp(W(a, b)*W(-b, -a))
        0

        """
        if isinstance(s, str):
            names = [x.name for x in symbols(s, seq=True)]
        else:
            raise ValueError('expecting a string')
        if len(names) == 1:
            return TensorHead(names[0], self, comm, matrix_behavior=matrix_behavior)
        return [TensorHead(name, self, comm, matrix_behavior=matrix_behavior) for name in names]


def tensorhead(name, typ, sym, comm=0, matrix_behavior=0):
    """
    Function generating tensorhead(s).

    Parameters
    ==========

    name : name or sequence of names (as in ``symbol``)

    typ :  index types

    sym :  same as ``*args`` in ``tensorsymmetry``

    comm : commutation group number
    see ``_TensorManager.set_comm``


    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b = tensor_indices('a b', Lorentz)
    >>> A = tensorhead('A', [Lorentz]*2, [[1]*2])
    >>> A(a, -b)
    A(a, -b)

    """
    sym = tensorsymmetry(*sym)
    S = TensorType(typ, sym)
    th = S(name, comm, matrix_behavior=matrix_behavior)
    return th


@doctest_depends_on(modules=('numpy',))
class TensorHead(Basic):
    r"""
    Tensor head of the tensor

    Parameters
    ==========

    name : name of the tensor

    typ : list of TensorIndexType

    comm : commutation group number

    Attributes
    ==========

    ``name``
    ``index_types``
    ``rank``
    ``types``  :  equal to ``typ.types``
    ``symmetry`` : equal to ``typ.symmetry``
    ``comm`` : int
        commutation group

    Notes
    =====

    A ``TensorHead`` belongs to a commutation group, defined by a
    symbol on number ``comm`` (see ``_TensorManager.set_comm``);
    tensors in a commutation group have the same commutation properties;
    by default ``comm`` is ``0``, the group of the commuting tensors.

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> A = tensorhead('A', [Lorentz, Lorentz], [[1], [1]])

    Examples with ndarray values, the components data assigned to the
    ``TensorHead`` object are assumed to be in a fully-contravariant
    representation. In case it is necessary to assign components data which
    represents the values of a non-fully covariant tensor, see the other
    examples.

    >>> Lorentz.data = [1, -1, -1, -1]
    >>> i0, i1 = tensor_indices('i0:2', Lorentz)
    >>> A.data = [[j+2*i for j in range(4)] for i in range(4)]

    in order to retrieve data, it is also necessary to specify abstract indices
    enclosed by round brackets, then numerical indices inside square brackets.

    >>> A(i0, i1)[0, 0]
    0
    >>> A(i0, i1)[2, 3] == 3+2*2
    True

    Notice that square brackets create a valued tensor expression instance:

    >>> A(i0, i1)
    A(i0, i1)

    To view the data, just type:

    >>> print(str(A.data))
    [[0 1 2 3]
     [2 3 4 5]
     [4 5 6 7]
     [6 7 8 9]]

    Turning to a tensor expression, covariant indices get the corresponding
    components data corrected by the metric:

    >>> print(str(A(i0, -i1).data))
    [[0 -1 -2 -3]
     [2 -3 -4 -5]
     [4 -5 -6 -7]
     [6 -7 -8 -9]]

    >>> print(str(A(-i0, -i1).data))
    [[0 -1 -2 -3]
     [-2 3 4 5]
     [-4 5 6 7]
     [-6 7 8 9]]

    while if all indices are contravariant, the ``ndarray`` remains the same

    >>> print(str(A(i0, i1).data))
     [[0 1 2 3]
     [2 3 4 5]
     [4 5 6 7]
     [6 7 8 9]]

    When all indices are contracted and components data are added to the tensor,
    accessing the data will return a scalar, no numpy object. In fact, numpy
    ndarrays are dropped to scalars if they contain only one element.

    >>> A(i0, -i0)
    A(L_0, -L_0)
    >>> A(i0, -i0).data
    -18

    It is also possible to assign components data to an indexed tensor, i.e. a
    tensor with specified covariant and contravariant components. In this
    example, the covariant components data of the Electromagnetic tensor are
    injected into `A`:

    >>> Ex, Ey, Ez, Bx, By, Bz = symbols('E_x E_y E_z B_x B_y B_z')
    >>> c = symbols('c', positive=True)

    Let's define `F`, an antisymmetric tensor, we have to assign an
    antisymmetric matrix to it, because `[[2]]` stands for the Young tableau
    representation of an antisymmetric set of two elements:

    >>> F = tensorhead('A', [Lorentz, Lorentz], [[2]])
    >>> F(-i0, -i1).data = [[0, Ex/c, Ey/c, Ez/c],
    ...                     [-Ex/c, 0, -Bz, By],
    ...                     [-Ey/c, Bz, 0, -Bx],
    ...                     [-Ez/c, -By, Bx, 0]]

    Now it is possible to retrieve the contravariant form of the Electromagnetic
    tensor:

    >>> print(str(F(i0, i1).data))
    [[0 -E_x/c -E_y/c -E_z/c]
     [E_x/c 0 -B_z B_y]
     [E_y/c B_z 0 -B_x]
     [E_z/c -B_y B_x 0]]

    and the mixed contravariant-covariant form:

    >>> print(str(F(i0, -i1).data))
    [[0 E_x/c E_y/c E_z/c]
     [E_x/c 0 B_z -B_y]
     [E_y/c -B_z 0 B_x]
     [E_z/c B_y -B_x 0]]

    To convert the numpy's ndarray to a diofant matrix, just cast:

    >>> Matrix(F.data)
    Matrix([
    [    0, -E_x/c, -E_y/c, -E_z/c],
    [E_x/c,      0,   -B_z,    B_y],
    [E_y/c,    B_z,      0,   -B_x],
    [E_z/c,   -B_y,    B_x,      0]])

    Still notice, in this last example, that accessing components data from a
    tensor without specifying the indices is equivalent to assume that all
    indices are contravariant.

    It is also possible to store symbolic components data inside a tensor, for
    example, define a four-momentum-like tensor:

    >>> P = tensorhead('P', [Lorentz], [[1]])
    >>> E, px, py, pz = symbols('E p_x p_y p_z', positive=True)
    >>> P.data = [E, px, py, pz]

    The contravariant and covariant components are, respectively:

    >>> print(str(P(i0).data))
    [E p_x p_y p_z]
    >>> print(str(P(-i0).data))
    [E -p_x -p_y -p_z]

    The contraction of a 1-index tensor by itself is usually indicated by a
    power by two:

    >>> P(i0)**2
    E**2 - p_x**2 - p_y**2 - p_z**2

    As the power by two is clearly identical to `P_\mu P^\mu`, it is possible to
    simply contract the ``TensorHead`` object, without specifying the indices

    >>> P**2
    E**2 - p_x**2 - p_y**2 - p_z**2

    """

    is_commutative = False

    def __new__(cls, name, typ, comm=0, matrix_behavior=0, **kw_args):
        if isinstance(name, str):
            name_symbol = Symbol(name)
        elif isinstance(name, Symbol):
            name_symbol = name
        else:
            raise ValueError('invalid name')

        comm2i = TensorManager.comm_symbols2i(comm)

        obj = Basic.__new__(cls, name_symbol, typ, **kw_args)

        obj._matrix_behavior = matrix_behavior

        obj._name = obj.args[0].name
        obj._rank = len(obj.index_types)
        obj._types = typ.types
        obj._symmetry = typ.symmetry
        obj._comm = comm2i
        return obj

    @property
    def name(self):
        return self._name

    @property
    def rank(self):
        return self._rank

    @property
    def types(self):
        return self._types[:]

    @property
    def symmetry(self):
        return self._symmetry

    @property
    def typ(self):
        return self.args[1]

    @property
    def comm(self):
        return self._comm

    @property
    def index_types(self):
        return self.args[1].index_types[:]

    def __lt__(self, other):
        return (self.name, self.index_types) < (other.name, other.index_types)

    def commutes_with(self, other):
        """
        Returns ``0`` if ``self`` and ``other`` commute, ``1`` if they anticommute.

        Returns ``None`` if ``self`` and ``other`` neither commute nor anticommute.

        """
        r = TensorManager.get_comm(self._comm, other._comm)
        return r

    def _print(self):
        return f"{self.name}({','.join([str(x) for x in self.index_types])})"

    def _check_auto_matrix_indices_in_call(self, *indices):
        matrix_behavior_kinds = {}

        if len(indices) != len(self.index_types):
            if not self._matrix_behavior:
                raise ValueError('wrong number of indices')

            # Take the last one or two missing
            # indices as auto-matrix indices:
            ldiff = len(self.index_types) - len(indices)
            if ldiff > 2:
                raise ValueError('wrong number of indices')
            if ldiff == 2:
                mat_ind = [len(indices), len(indices) + 1]
            elif ldiff == 1:
                mat_ind = [len(indices)]
            not_equal = True
        else:
            not_equal = False
            mat_ind = [i for i, e in enumerate(indices) if e is True]
            if mat_ind:
                not_equal = True
            indices = tuple(_ for _ in indices if _ is not True)

            for i, el in enumerate(indices):
                if not isinstance(el, TensorIndex):
                    not_equal = True
                    break
                if el._tensortype != self.index_types[i]:
                    not_equal = True
                    break

        if not_equal:
            for el in mat_ind:
                eltyp = self.index_types[el]
                if eltyp in matrix_behavior_kinds:
                    elind = -self.index_types[el].auto_right
                    matrix_behavior_kinds[eltyp].append(elind)
                else:
                    elind = self.index_types[el].auto_left
                    matrix_behavior_kinds[eltyp] = [elind]
                indices = indices[:el] + (elind,) + indices[el:]

        return indices, matrix_behavior_kinds

    def __call__(self, *indices, **kw_args):
        """
        Returns a tensor with indices.

        There is a special behavior in case of indices denoted by ``True``,
        they are considered auto-matrix indices, their slots are automatically
        filled, and confer to the tensor the behavior of a matrix or vector
        upon multiplication with another tensor containing auto-matrix indices
        of the same ``TensorIndexType``. This means indices get summed over the
        same way as in matrix multiplication. For matrix behavior, define two
        auto-matrix indices, for vector behavior define just one.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b = tensor_indices('a b', Lorentz)
        >>> A = tensorhead('A', [Lorentz]*2, [[1]*2])
        >>> t = A(a, -b)
        >>> t
        A(a, -b)

        To use the auto-matrix index behavior, just put a ``True`` on the
        desired index position.

        >>> r = A(True, True)
        >>> r
        A(auto_left, -auto_right)

        Here ``auto_left`` and ``auto_right`` are automatically generated
        tensor indices, they are only two for every ``TensorIndexType`` and
        can be assigned to just one or two indices of a given type.

        Auto-matrix indices can be assigned many times in a tensor, if indices
        are of different ``TensorIndexType``

        >>> Spinor = TensorIndexType('Spinor', dummy_fmt='S')
        >>> B = tensorhead('B', [Lorentz, Lorentz, Spinor, Spinor], [[1]*4])
        >>> s = B(True, True, True, True)
        >>> s
        B(auto_left, -auto_right, auto_left, -auto_right)

        Here, ``auto_left`` and ``auto_right`` are repeated twice, but they are
        not the same indices, as they refer to different ``TensorIndexType``s.

        Auto-matrix indices are automatically contracted upon multiplication,

        >>> r*s
        A(auto_left, L_0)*B(-L_0, -auto_right, auto_left, -auto_right)

        The multiplication algorithm has found an ``auto_right`` index in ``A``
        and an ``auto_left`` index in ``B`` referring to the same
        ``TensorIndexType`` (``Lorentz``), so they have been contracted.

        Auto-matrix indices can be accessed from the ``TensorIndexType``:

        >>> Lorentz.auto_right
        auto_right
        >>> Lorentz.auto_left
        auto_left

        There is a special case, in which the ``True`` parameter is not needed
        to declare an auto-matrix index, i.e. when the matrix behavior has been
        declared upon ``TensorHead`` construction, in that case the last one or
        two tensor indices may be omitted, so that they automatically become
        auto-matrix indices:

        >>> C = tensorhead('C', [Lorentz, Lorentz], [[1]*2], matrix_behavior=True)
        >>> C()
        C(auto_left, -auto_right)

        """
        indices, _ = self._check_auto_matrix_indices_in_call(*indices)
        tensor = Tensor._new_with_dummy_replacement(self, indices, **kw_args)
        return tensor

    def __pow__(self, other):
        if self.data is None:
            raise ValueError('No power on abstract tensors.')
        numpy = import_module('numpy')
        metrics = [_.data for _ in self.args[1].args[0]]

        marray = self.data
        for metric in metrics:
            marray = numpy.tensordot(marray, numpy.tensordot(metric, marray, (1, 0)), (0, 0))
        pow2 = marray[()]
        return pow2 ** (Rational(1, 2) * other)

    @property
    def data(self):
        return _tensor_data_substitution_dict[self]

    @data.setter
    def data(self, data):
        _tensor_data_substitution_dict[self] = data

    @data.deleter
    def data(self):
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]

    def __iter__(self):
        return self.data.flatten().__iter__()


@doctest_depends_on(modules=('numpy',))
class TensExpr(Basic):
    """
    Abstract base class for tensor expressions

    Notes
    =====

    A tensor expression is an expression formed by tensors;
    currently the sums of tensors are distributed.

    A ``TensExpr`` can be a ``TensAdd`` or a ``TensMul``.

    ``TensAdd`` objects are put in canonical form using the Butler-Portugal
    algorithm for canonicalization under monoterm symmetries.

    ``TensMul`` objects are formed by products of component tensors,
    and include a coefficient, which is a Diofant expression.


    In the internal representation contracted indices are represented
    by ``(ipos1, ipos2, icomp1, icomp2)``, where ``icomp1`` is the position
    of the component tensor with contravariant index, ``ipos1`` is the
    slot which the index occupies in that component tensor.

    Contracted indices are therefore nameless in the internal representation.

    """

    _op_priority = 11.0
    is_commutative = False

    def __neg__(self):
        return self*Integer(-1)

    def __abs__(self):
        raise NotImplementedError

    def __add__(self, other):
        raise NotImplementedError

    def __radd__(self, other):
        raise NotImplementedError

    def __sub__(self, other):
        raise NotImplementedError

    def __rsub__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    def __pow__(self, other):
        if self.data is None:
            raise ValueError('No power without ndarray data.')
        numpy = import_module('numpy')
        free = self.free

        marray = self.data
        for metric in free:
            marray = numpy.tensordot(
                marray,
                numpy.tensordot(
                    metric[0]._tensortype.data,
                    marray,
                    (1, 0)
                ),
                (0, 0)
            )
        pow2 = marray[()]
        return pow2 ** (Rational(1, 2) * other)

    def __rpow__(self, other):
        raise NotImplementedError

    def __truediv__(self, other):
        raise NotImplementedError

    def __rtruediv__(self, other):
        raise NotImplementedError()

    @doctest_depends_on(modules=('numpy',))
    def get_matrix(self):
        """
        Returns ndarray components data as a matrix, if components data are
        available and ndarray dimension does not exceed 2.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> sym2 = tensorsymmetry([1]*2)
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> A = S2('A')

        The tensor ``A`` is symmetric in its indices, as can be deduced by the
        ``[1, 1]`` Young tableau when constructing `sym2`. One has to be
        careful to assign symmetric component data to ``A``, as the symmetry
        properties of data are currently not checked to be compatible with the
        defined tensor symmetry.

        >>> Lorentz.data = [1, -1, -1, -1]
        >>> i0, i1 = tensor_indices('i0:2', Lorentz)
        >>> A.data = [[j+i for j in range(4)] for i in range(4)]
        >>> A(i0, i1).get_matrix()
        Matrix([
        [0, 1, 2, 3],
        [1, 2, 3, 4],
        [2, 3, 4, 5],
        [3, 4, 5, 6]])

        It is possible to perform usual operation on matrices, such as the
        matrix multiplication:

        >>> A(i0, i1).get_matrix()*ones(4, 1)
        Matrix([
        [ 6],
        [10],
        [14],
        [18]])

        >>> del A.data

        """
        if 0 < self.rank <= 2:
            rows = self.data.shape[0]
            columns = self.data.shape[1] if self.rank == 2 else 1
            if self.rank == 2:
                mat_list = [] * rows
                for i in range(rows):
                    mat_list.append([])
                    for j in range(columns):
                        mat_list[i].append(self[i, j])
            else:
                mat_list = [None] * rows
                for i in range(rows):
                    mat_list[i] = self[i]
            return Matrix(mat_list)
        raise NotImplementedError('missing multidimensional reduction to matrix.')


@doctest_depends_on(modules=('numpy',))
class TensAdd(TensExpr):
    """
    Sum of tensors

    Parameters
    ==========

    free_args : list of the free indices

    Attributes
    ==========

    ``args`` : tuple
        of addends
    ``rank`` : tuple
        rank of the tensor
    ``free_args`` : list
        of the free indices in sorted order

    Notes
    =====

    Sum of more than one tensor are put automatically in canonical form.

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b = tensor_indices('a b', Lorentz)
    >>> p, q = tensorhead('p q', [Lorentz], [[1]])
    >>> t = p(a) + q(a)
    >>> t
    p(a) + q(a)
    >>> t(b)
    p(b) + q(b)

    Examples with components data added to the tensor expression:

    >>> Lorentz.data = [1, -1, -1, -1]
    >>> a, b = tensor_indices('a, b', Lorentz)
    >>> p.data = [2, 3, -2, 7]
    >>> q.data = [2, 3, -2, 7]
    >>> t = p(a) + q(a)
    >>> t
    p(a) + q(a)
    >>> t(b)
    p(b) + q(b)

    The following are: 2**2 - 3**2 - 2**2 - 7**2 ==> -58

    >>> (p(a)*p(-a)).data
    -58
    >>> p(a)**2
    -58

    """

    def __new__(cls, *args, **kw_args):
        args = [sympify(x) for x in args if x]
        args = TensAdd._tensAdd_flatten(args)

        if not args:
            return Integer(0)

        if len(args) == 1 and not isinstance(args[0], TensExpr):
            return args[0]

        # replace auto-matrix indices so that they are the same in all addends
        args = TensAdd._tensAdd_check_automatrix(args)

        # now check that all addends have the same indices:
        TensAdd._tensAdd_check(args)

        # if TensAdd has only 1 TensMul element in its `args`:
        if len(args) == 1 and isinstance(args[0], TensMul):
            obj = Basic.__new__(cls, *args, **kw_args)
            return obj

        # TODO: do not or do canonicalize by default?
        # Technically, one may wish to have additions of non-canonicalized
        # tensors. This feature should be removed in the future.
        # Unfortunately this would require to rewrite a lot of tests.
        # canonicalize all TensMul
        args = [canon_bp(x) for x in args if x]
        args = [x for x in args if x]

        # if there are no more args (i.e. have cancelled out),
        # just return zero:
        if not args:
            return Integer(0)

        if len(args) == 1:
            return args[0]

        # collect canonicalized terms
        def sort_key(t):
            x = get_tids(t)
            return x.components, x.free, x.dum
        args.sort(key=sort_key)
        args = TensAdd._tensAdd_collect_terms(args)
        if not args:
            return Integer(0)
        # it there is only a component tensor return it
        if len(args) == 1:
            return args[0]

        obj = Basic.__new__(cls, *args, **kw_args)
        return obj

    @staticmethod
    def _tensAdd_flatten(args):
        # flatten TensAdd, coerce terms which are not tensors to tensors

        if not all(isinstance(x, TensExpr) for x in args):
            args_expanded = []
            for x in args:
                if isinstance(x, TensAdd):
                    args_expanded.extend(list(x.args))
                else:
                    args_expanded.append(x)
            args_tensor = []
            args_scalar = []
            for x in args_expanded:
                if isinstance(x, TensExpr) and x.coeff:
                    args_tensor.append(x)
                if not isinstance(x, TensExpr):
                    args_scalar.append(x)
            t1 = TensMul.from_data(Add(*args_scalar), [], [], [])
            args = [t1] + args_tensor
        a = []
        for x in args:
            if isinstance(x, TensAdd):
                a.extend(list(x.args))
            else:
                a.append(x)
        args = [x for x in a if x.coeff]
        return args

    @staticmethod
    def _tensAdd_check_automatrix(args):
        # check that all automatrix indices are the same.

        # if there are no addends, just return.
        if not args:
            return args

        # @type auto_left_types: set
        auto_left_types = set()
        auto_right_types = set()
        args_auto_left_types = []
        args_auto_right_types = []
        for i, arg in enumerate(args):
            arg_auto_left_types = set()
            arg_auto_right_types = set()
            for index in get_indices(arg):
                # @type index: TensorIndex
                if index in (index._tensortype.auto_left, -index._tensortype.auto_left):
                    auto_left_types.add(index._tensortype)
                    arg_auto_left_types.add(index._tensortype)
                if index in (index._tensortype.auto_right, -index._tensortype.auto_right):
                    auto_right_types.add(index._tensortype)
                    arg_auto_right_types.add(index._tensortype)
            args_auto_left_types.append(arg_auto_left_types)
            args_auto_right_types.append(arg_auto_right_types)
        for arg, aas_left, aas_right in zip(args, args_auto_left_types, args_auto_right_types):
            missing_left = auto_left_types - aas_left
            missing_right = auto_right_types - aas_right
            missing_intersection = missing_left & missing_right
            for j in missing_intersection:
                args[i] *= j.delta(j.auto_left, -j.auto_right)
            if missing_left != missing_right:
                raise ValueError('cannot determine how to add auto-matrix indices on some args')

        return args

    @staticmethod
    def _tensAdd_check(args):
        # check that all addends have the same free indices
        indices0 = {x[0] for x in get_tids(args[0]).free}
        list_indices = [{y[0] for y in get_tids(x).free} for x in args[1:]]
        if not all(x == indices0 for x in list_indices):
            raise ValueError('all tensors must have the same indices')

    @staticmethod
    def _tensAdd_collect_terms(args):
        # collect TensMul terms differing at most by their coefficient
        a = []
        prev = args[0]
        prev_coeff = get_coeff(prev)
        changed = False

        for x in args[1:]:
            # if x and prev have the same tensor, update the coeff of prev
            x_tids = get_tids(x)
            prev_tids = get_tids(prev)
            if x_tids.components == prev_tids.components \
                    and x_tids.free == prev_tids.free and x_tids.dum == prev_tids.dum:
                prev_coeff = prev_coeff + get_coeff(x)
                changed = True
                op = 0
            else:
                # x and prev are different; if not changed, prev has not
                # been updated; store it
                if not changed:
                    a.append(prev)
                else:
                    # get a tensor from prev with coeff=prev_coeff and store it
                    if prev_coeff:
                        t = TensMul.from_data(prev_coeff, prev_tids.components,
                                              prev_tids.free, prev_tids.dum)
                        a.append(t)
                # move x to prev
                op = 1
                prev = x
                prev_coeff = get_coeff(x)
                changed = False
        # if the case op=0 prev was not stored; store it now
        # in the case op=1 x was not stored; store it now (as prev)
        if op == 0 and prev_coeff:
            prev = TensMul.from_data(prev_coeff, prev_tids.components, prev_tids.free, prev_tids.dum)
            a.append(prev)
        elif op == 1:
            a.append(prev)
        return a

    @property
    def rank(self):
        return self.args[0].rank

    @property
    def free_args(self):
        return self.args[0].free_args

    def __call__(self, *indices):
        """Returns tensor with ordered free indices replaced by ``indices``

        Parameters
        ==========

        indices

        Examples
        ========

        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> i0, i1, i2, i3, i4 = tensor_indices('i0:5', Lorentz)
        >>> p, q = tensorhead('p q', [Lorentz], [[1]])
        >>> g = Lorentz.metric
        >>> t = p(i0)*p(i1) + g(i0, i1)*q(i2)*q(-i2)
        >>> t(i0, i2)
        metric(i0, i2)*q(L_0)*q(-L_0) + p(i0)*p(i2)
        >>> t(i0, i1) - t(i1, i0)
        0

        """
        free_args = self.free_args
        indices = list(indices)
        if [x._tensortype for x in indices] != [x._tensortype for x in free_args]:
            raise ValueError('incompatible types')
        if indices == free_args:
            return self
        index_tuples = list(zip(free_args, indices))
        a = [x.func(*x.fun_eval(*index_tuples).args) for x in self.args]
        res = TensAdd(*a)

        return res

    def canon_bp(self):
        """
        Canonicalize using the Butler-Portugal algorithm for canonicalization
        under monoterm symmetries.

        """
        args = [x.canon_bp() for x in self.args]
        res = TensAdd(*args)
        return res

    def equals(self, other):
        other = sympify(other)
        if isinstance(other, TensMul) and other._coeff == 0:
            return all(x._coeff == 0 for x in self.args)
        if isinstance(other, TensExpr):
            if self.rank != other.rank:
                return False
        if isinstance(other, TensAdd):
            if set(self.args) != set(other.args):
                return False
            return True
        t = self - other
        if not isinstance(t, TensExpr):
            return t == 0
        if isinstance(t, TensMul):
            return t._coeff == 0
        return all(x._coeff == 0 for x in t.args)

    def __add__(self, other):
        return TensAdd(self, other)

    def __radd__(self, other):
        return TensAdd(other, self)

    def __sub__(self, other):
        return TensAdd(self, -other)

    def __rsub__(self, other):
        return TensAdd(other, -self)

    def __mul__(self, other):
        return TensAdd(*(x*other for x in self.args))

    def __rmul__(self, other):
        return self*other

    def __truediv__(self, other):
        other = sympify(other)
        if isinstance(other, TensExpr):
            raise ValueError('cannot divide by a tensor')
        return TensAdd(*(x/other for x in self.args))

    def __rtruediv__(self, other):
        raise ValueError('cannot divide by a tensor')

    def __getitem__(self, item):
        # pylint: disable=unsubscriptable-object
        return self.data[item]

    def contract_delta(self, delta):
        args = [x.contract_delta(delta) for x in self.args]
        t = TensAdd(*args)
        return canon_bp(t)

    def contract_metric(self, g):
        """
        Raise or lower indices with the metric ``g``

        Parameters
        ==========

        g :  metric

        contract_all : if True, eliminate all ``g`` which are contracted

        Notes
        =====

        See Also
        ========

        TensorIndexType

        """
        args = [contract_metric(x, g) for x in self.args]
        t = TensAdd(*args)
        return canon_bp(t)

    def fun_eval(self, *index_tuples):
        """
        Return a tensor with free indices substituted according to ``index_tuples``

        Parameters
        ==========

        index_types : list of tuples ``(old_index, new_index)``

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i, j, k, l = tensor_indices('i j k l', Lorentz)
        >>> A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
        >>> t = A(i, k)*B(-k, -j) + A(i, -j)
        >>> t.fun_eval((i, k), (-j, l))
        A(k, L_0)*B(l, -L_0) + A(k, l)

        """
        args = self.args
        args1 = []
        for x in args:
            y = x.fun_eval(*index_tuples)
            args1.append(y)
        return TensAdd(*args1)

    def substitute_indices(self, *index_tuples):
        """
        Return a tensor with free indices substituted according to ``index_tuples``

        Parameters
        ==========

        index_types : list of tuples ``(old_index, new_index)``

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i, j, k, l = tensor_indices('i j k l', Lorentz)
        >>> A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
        >>> t = A(i, k)*B(-k, -j)
        >>> t
        A(i, L_0)*B(-L_0, -j)
        >>> t.substitute_indices((i, j), (j, k))
        A(j, L_0)*B(-L_0, -k)

        """
        args = self.args
        args1 = []
        for x in args:
            y = x.substitute_indices(*index_tuples)
            args1.append(y)
        return TensAdd(*args1)

    def _print(self):
        a = []
        args = self.args
        for x in args:
            a.append(str(x))
        a.sort()
        s = ' + '.join(a)
        s = s.replace('+ -', '- ')
        return s

    @staticmethod
    def from_TIDS_list(coeff, tids_list):
        """
        Given a list of coefficients and a list of ``TIDS`` objects, construct
        a ``TensAdd`` instance, equivalent to the one that would result from
        creating single instances of ``TensMul`` and then adding them.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i, j = tensor_indices('i j', Lorentz)
        >>> A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
        >>> ea = 3*A(i, j)
        >>> eb = 2*B(j, i)
        >>> t1 = ea._tids
        >>> t2 = eb._tids
        >>> c1 = ea.coeff
        >>> c2 = eb.coeff
        >>> TensAdd.from_TIDS_list([c1, c2], [t1, t2])
        2*B(i, j) + 3*A(i, j)

        If the coefficient parameter is a scalar, then it will be applied
        as a coefficient on all ``TIDS`` objects.

        >>> TensAdd.from_TIDS_list(4, [t1, t2])
        4*A(i, j) + 4*B(i, j)

        """
        if not isinstance(coeff, (list, tuple, Tuple)):
            coeff = [coeff] * len(tids_list)
        tensmul_list = [TensMul.from_TIDS(c, t) for c, t in zip(coeff, tids_list)]
        return TensAdd(*tensmul_list)

    @property
    def data(self):
        return _tensor_data_substitution_dict[self]


@doctest_depends_on(modules=('numpy',))
class Tensor(TensExpr):
    """
    Base tensor class, i.e. this represents a tensor, the single unit to be
    put into an expression.

    This object is usually created from a ``TensorHead``, by attaching indices
    to it. Indices preceded by a minus sign are considered contravariant,
    otherwise covariant.

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> mu, nu = tensor_indices('mu nu', Lorentz)
    >>> A = tensorhead('A', [Lorentz, Lorentz], [[1], [1]])
    >>> A(mu, -nu)
    A(mu, -nu)
    >>> A(mu, -mu)
    A(L_0, -L_0)

    """

    is_commutative = False

    def __new__(cls, tensor_head, indices, **kw_args):
        tids = TIDS.from_components_and_indices((tensor_head,), indices)
        obj = Basic.__new__(cls, tensor_head, Tuple(*indices), **kw_args)
        obj._tids = tids
        obj._indices = indices
        obj._is_canon_bp = kw_args.get('is_canon_bp', False)
        return obj

    @staticmethod
    def _new_with_dummy_replacement(tensor_head, indices, **kw_args):
        tids = TIDS.from_components_and_indices((tensor_head,), indices)
        indices = tids.get_indices()
        return Tensor(tensor_head, indices, **kw_args)

    @property
    def is_canon_bp(self):
        return self._is_canon_bp

    @property
    def indices(self):
        return self._indices

    @property
    def free(self):
        return self._tids.free

    @property
    def dum(self):
        return self._tids.dum

    @property
    def rank(self):
        return len(self.free)

    @property
    def free_args(self):
        return sorted(x[0] for x in self.free)

    def perm2tensor(self, g, canon_bp=False):
        """
        Returns the tensor corresponding to the permutation ``g``

        For further details, see the method in ``TIDS`` with the same name.

        """
        return perm2tensor(self, g, canon_bp)

    def canon_bp(self):
        if self._is_canon_bp:
            return self
        g, dummies, msym, v = self._tids.canon_args()
        can = canonicalize(g, dummies, msym, *v)
        if can == 0:
            return Integer(0)
        tensor = self.perm2tensor(can, True)
        return tensor

    @property
    def types(self):
        return get_tids(self).components[0]._types

    @property
    def coeff(self):
        return Integer(1)

    @property
    def component(self):
        return self.args[0]

    @property
    def components(self):
        return [self.args[0]]

    def split(self):
        return [self]

    def get_indices(self):
        """Get a list of indices, corresponding to those of the tensor."""
        return self._tids.get_indices()

    def substitute_indices(self, *index_tuples):
        return substitute_indices(self, *index_tuples)

    def __call__(self, *indices):
        """Returns tensor with ordered free indices replaced by ``indices``

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i0, i1, i2, i3, i4 = tensor_indices('i0:5', Lorentz)
        >>> A = tensorhead('A', [Lorentz]*5, [[1]*5])
        >>> t = A(i2, i1, -i2, -i3, i4)
        >>> t
        A(L_0, i1, -L_0, -i3, i4)
        >>> t(i1, i2, i3)
        A(L_0, i1, -L_0, i2, i3)

        """
        free_args = self.free_args
        indices = list(indices)
        if [x._tensortype for x in indices] != [x._tensortype for x in free_args]:
            raise ValueError('incompatible types')
        if indices == free_args:
            return self
        t = self.fun_eval(*list(zip(free_args, indices)))

        # object is rebuilt in order to make sure that all contracted indices
        # get recognized as dummies, but only if there are contracted indices.
        if len({i if i.is_up else -i for i in indices}) != len(indices):
            return t.func(*t.args)
        return t

    def fun_eval(self, *index_tuples):
        free = self.free
        free1 = []
        for j, ipos, cpos in free:
            # search j in index_tuples
            for i, v in index_tuples:
                if i == j:
                    free1.append((v, ipos, cpos))
                    break
            else:
                free1.append((j, ipos, cpos))
        return TensMul.from_data(self.coeff, self.components, free1, self.dum)

    # TODO: put this into TensExpr?
    def __iter__(self):
        return self.data.flatten().__iter__()

    # TODO: put this into TensExpr?
    def __getitem__(self, item):
        return self.data[item]

    @property
    def data(self):
        return _tensor_data_substitution_dict[self]

    @data.setter
    def data(self, data):
        # TODO: check data compatibility with properties of tensor.
        _tensor_data_substitution_dict[self] = data

    def __mul__(self, other):
        if isinstance(other, TensAdd):
            return TensAdd(*[self*arg for arg in other.args])
        tmul = TensMul(self, other)
        return tmul

    def __rmul__(self, other):
        return TensMul(other, self)

    def __truediv__(self, other):
        if isinstance(other, TensExpr):
            raise ValueError('cannot divide by a tensor')
        return TensMul(self, Integer(1)/other, is_canon_bp=self.is_canon_bp)

    def __rtruediv__(self, other):
        raise ValueError('cannot divide by a tensor')

    def __add__(self, other):
        return TensAdd(self, other)

    def __radd__(self, other):
        return TensAdd(other, self)

    def __sub__(self, other):
        return TensAdd(self, -other)

    def __neg__(self):
        return TensMul(Integer(-1), self)

    def _print(self):
        indices = [str(ind) for ind in self.indices]
        component = self.component
        if component.rank > 0:
            return f"{component.name}({', '.join(indices)})"
        return f'{component.name}'

    def equals(self, other):
        if other == 0:
            return self.coeff == 0
        other = sympify(other)
        if not isinstance(other, TensExpr):
            assert not self.components
            return Integer(1) == other

        def _get_compar_comp(self):
            t = self.canon_bp()
            r = (t.coeff, tuple(t.components),
                 tuple(sorted(t.free)), tuple(sorted(t.dum)))
            return r

        return _get_compar_comp(self) == _get_compar_comp(other)

    def contract_metric(self, metric):
        tids, sign = get_tids(self).contract_metric(metric)
        return TensMul.from_TIDS(sign, tids)


@doctest_depends_on(modules=('numpy',))
class TensMul(TensExpr):
    """
    Product of tensors

    Parameters
    ==========

    coeff : Diofant coefficient of the tensor
    args

    Attributes
    ==========

    ``components`` : list of ``TensorHead`` of the component tensors
    ``types`` : list of nonrepeated ``TensorIndexType``
    ``free`` : list of ``(ind, ipos, icomp)``, see Notes
    ``dum`` : list of ``(ipos1, ipos2, icomp1, icomp2)``, see Notes
    ``ext_rank`` : tuple
        rank of the tensor counting the dummy indices
    ``rank`` : tuple
        rank of the tensor
    ``coeff`` : Expr
        Diofant coefficient of the tensor
    ``free_args`` : list
        list of the free indices in sorted order
    ``is_canon_bp`` : ``True`` if the tensor in in canonical form

    Notes
    =====

    ``args[0]``   list of ``TensorHead`` of the component tensors.

    ``args[1]``   list of ``(ind, ipos, icomp)``
    where ``ind`` is a free index, ``ipos`` is the slot position
    of ``ind`` in the ``icomp``-th component tensor.

    ``args[2]`` list of tuples representing dummy indices.
    ``(ipos1, ipos2, icomp1, icomp2)`` indicates that the contravariant
    dummy index is the ``ipos1``-th slot position in the ``icomp1``-th
    component tensor; the corresponding covariant index is
    in the ``ipos2`` slot position in the ``icomp2``-th component tensor.

    """

    def __new__(cls, *args, **kw_args):
        # make sure everything is sympified:
        args = [sympify(arg) for arg in args]

        # flatten:
        args = TensMul._flatten(args)

        is_canon_bp = kw_args.get('is_canon_bp', False)
        if not any(isinstance(arg, TensExpr) for arg in args):
            tids = TIDS([], [], [])
        else:
            tids_list = [arg._tids for arg in args if isinstance(arg, (Tensor, TensMul))]
            if len(tids_list) == 1:
                for arg in args:
                    if not isinstance(arg, Tensor):
                        continue
                    is_canon_bp = kw_args.get('is_canon_bp', arg._is_canon_bp)
            tids = functools.reduce(lambda a, b: a*b, tids_list)

        coeff = functools.reduce(lambda a, b: a*b, [Integer(1)] + [arg for arg in args if not isinstance(arg, TensExpr)])
        args = tids.get_tensors()
        if coeff != 1:
            args = [coeff] + args
        if len(args) == 1:
            return args[0]

        obj = Basic.__new__(cls, *args)
        obj._types = []
        for t in tids.components:
            obj._types.extend(t._types)
        obj._tids = tids
        obj._ext_rank = len(obj._tids.free) + 2*len(obj._tids.dum)
        obj._coeff = coeff
        obj._is_canon_bp = is_canon_bp
        return obj

    @staticmethod
    def _flatten(args):
        a = []
        for arg in args:
            if isinstance(arg, TensMul):
                a.extend(arg.args)
            else:
                a.append(arg)
        return a

    @staticmethod
    def from_data(coeff, components, free, dum, **kw_args):
        tids = TIDS(components, free, dum)
        return TensMul.from_TIDS(coeff, tids, **kw_args)

    @staticmethod
    def from_TIDS(coeff, tids, **kw_args):
        return TensMul(coeff, *tids.get_tensors(), **kw_args)

    @property
    def free_args(self):
        return sorted(x[0] for x in self.free)

    @property
    def components(self):
        return self._tids.components[:]

    @property
    def free(self):
        return self._tids.free[:]

    @property
    def coeff(self):
        return self._coeff

    @property
    def dum(self):
        return self._tids.dum[:]

    @property
    def rank(self):
        return len(self.free)

    @property
    def types(self):
        return self._types[:]

    def equals(self, other):
        if other == 0:
            return self.coeff == 0
        other = sympify(other)
        if not isinstance(other, TensExpr):
            assert not self.components
            return self._coeff == other

        def _get_compar_comp(self):
            t = self.canon_bp()
            r = (get_coeff(t), tuple(t.components),
                 tuple(sorted(t.free)), tuple(sorted(t.dum)))
            return r

        return _get_compar_comp(self) == _get_compar_comp(other)

    def get_indices(self):
        """
        Returns the list of indices of the tensor

        The indices are listed in the order in which they appear in the
        component tensors.
        The dummy indices are given a name which does not collide with
        the names of the free indices.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0 m1 m2', Lorentz)
        >>> g = Lorentz.metric
        >>> p, q = tensorhead('p q', [Lorentz], [[1]])
        >>> t = p(m1)*g(m0, m2)
        >>> t.get_indices()
        [m1, m0, m2]

        """
        return self._tids.get_indices()

    def split(self):
        """
        Returns a list of tensors, whose product is ``self``

        Dummy indices contracted among different tensor components
        become free indices with the same name as the one used to
        represent the dummy indices.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b, c, d = tensor_indices('a b c d', Lorentz)
        >>> A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
        >>> t = A(a, b)*B(-b, c)
        >>> t
        A(a, L_0)*B(-L_0, c)
        >>> t.split()
        [A(a, L_0), B(-L_0, c)]

        """
        if not self.args:
            return [self]
        splitp = []
        res = 1
        for arg in self.args:
            if isinstance(arg, Tensor):
                splitp.append(res*arg)
                res = 1
            else:
                res *= arg
        return splitp

    def __add__(self, other):
        return TensAdd(self, other)

    def __radd__(self, other):
        return TensAdd(other, self)

    def __sub__(self, other):
        return TensAdd(self, -other)

    def __rsub__(self, other):
        return TensAdd(other, -self)

    def __mul__(self, other):
        """
        Multiply two tensors using Einstein summation convention.

        If the two tensors have an index in common, one contravariant
        and the other covariant, in their product the indices are summed

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0 m1 m2', Lorentz)
        >>> g = Lorentz.metric
        >>> p, q = tensorhead('p q', [Lorentz], [[1]])
        >>> t1 = p(m0)
        >>> t2 = q(-m0)
        >>> t1*t2
        p(L_0)*q(-L_0)

        """
        other = sympify(other)
        if not isinstance(other, TensExpr):
            coeff = self.coeff*other
            tmul = TensMul.from_TIDS(coeff, self._tids, is_canon_bp=self._is_canon_bp)
            return tmul
        if isinstance(other, TensAdd):
            return TensAdd(*[self*x for x in other.args])

        new_tids = self._tids*other._tids
        coeff = self.coeff*other.coeff
        tmul = TensMul.from_TIDS(coeff, new_tids)
        return tmul

    def __rmul__(self, other):
        other = sympify(other)
        coeff = other*self._coeff
        tmul = TensMul.from_TIDS(coeff, self._tids)
        return tmul

    def __truediv__(self, other):
        other = sympify(other)
        if isinstance(other, TensExpr):
            raise ValueError('cannot divide by a tensor')
        coeff = self._coeff/other
        tmul = TensMul.from_TIDS(coeff, self._tids, is_canon_bp=self._is_canon_bp)
        return tmul

    def __getitem__(self, item):
        return self.data[item]

    def sorted_components(self):
        """
        Returns a tensor with sorted components
        calling the corresponding method in a ``TIDS`` object.

        """
        new_tids, sign = self._tids.sorted_components()
        coeff = -self.coeff if sign == -1 else self.coeff
        t = TensMul.from_TIDS(coeff, new_tids)
        return t

    def perm2tensor(self, g, canon_bp=False):
        """
        Returns the tensor corresponding to the permutation ``g``

        For further details, see the method in ``TIDS`` with the same name.

        """
        return perm2tensor(self, g, canon_bp)

    def canon_bp(self):
        """
        Canonicalize using the Butler-Portugal algorithm for canonicalization
        under monoterm symmetries.

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0 m1 m2', Lorentz)
        >>> A = tensorhead('A', [Lorentz]*2, [[2]])
        >>> t = A(m0, -m1)*A(m1, -m0)
        >>> t.canon_bp()
        -A(L_0, L_1)*A(-L_0, -L_1)
        >>> t = A(m0, -m1)*A(m1, -m2)*A(m2, -m0)
        >>> t.canon_bp()
        0

        """
        if self._is_canon_bp:
            return self
        if not self.components:
            return self
        t = self.sorted_components()
        g, dummies, msym, v = t._tids.canon_args()
        can = canonicalize(g, dummies, msym, *v)
        if can == 0:
            return Integer(0)
        tmul = t.perm2tensor(can, True)
        return tmul

    def contract_delta(self, delta):
        t = self.contract_metric(delta)
        return t

    def contract_metric(self, g):
        """
        Raise or lower indices with the metric ``g``

        Parameters
        ==========

        g : metric

        Notes
        =====

        See Also
        ========

        TensorIndexType

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0 m1 m2', Lorentz)
        >>> g = Lorentz.metric
        >>> p, q = tensorhead('p q', [Lorentz], [[1]])
        >>> t = p(m0)*q(m1)*g(-m0, -m1)
        >>> t.canon_bp()
        metric(L_0, L_1)*p(-L_0)*q(-L_1)
        >>> t.contract_metric(g).canon_bp()
        p(L_0)*q(-L_0)

        """
        tids, sign = get_tids(self).contract_metric(g)
        res = TensMul.from_TIDS(sign*self.coeff, tids)
        return res

    def substitute_indices(self, *index_tuples):
        return substitute_indices(self, *index_tuples)

    def fun_eval(self, *index_tuples):
        """
        Return a tensor with free indices substituted according to ``index_tuples``

        ``index_types`` list of tuples ``(old_index, new_index)``

        Examples
        ========

        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i, j, k, l = tensor_indices('i j k l', Lorentz)
        >>> A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
        >>> t = A(i, k)*B(-k, -j)
        >>> t
        A(i, L_0)*B(-L_0, -j)
        >>> t.fun_eval((i, k), (-j, l))
        A(k, L_0)*B(-L_0, l)

        """
        free = self.free
        free1 = []
        for j, ipos, cpos in free:
            # search j in index_tuples
            for i, v in index_tuples:
                if i == j:
                    free1.append((v, ipos, cpos))
                    break
            else:
                free1.append((j, ipos, cpos))
        return TensMul.from_data(self.coeff, self.components, free1, self.dum)

    def __call__(self, *indices):
        """Returns tensor product with ordered free indices replaced by ``indices``

        Examples
        ========

        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> i0, i1, i2, i3, i4 = tensor_indices('i0:5', Lorentz)
        >>> g = Lorentz.metric
        >>> p, q = tensorhead('p q', [Lorentz], [[1]])
        >>> t = p(i0)*q(i1)*q(-i1)
        >>> t(i1)
        p(i1)*q(L_0)*q(-L_0)

        """
        free_args = self.free_args
        indices = list(indices)
        if [x._tensortype for x in indices] != [x._tensortype for x in free_args]:
            raise ValueError('incompatible types')
        if indices == free_args:
            return self
        t = self.fun_eval(*list(zip(free_args, indices)))

        # object is rebuilt in order to make sure that all contracted indices
        # get recognized as dummies, but only if there are contracted indices.
        if len({i if i.is_up else -i for i in indices}) != len(indices):
            return t.func(*t.args)
        return t

    def _print(self):
        args = self.args

        def get_str(arg):
            return str(arg) if arg.is_Atom or isinstance(arg, TensExpr) else f'({arg!s})'

        if not args:
            # no arguments is equivalent to "1", i.e. TensMul().
            # If tensors are constructed correctly, this should never occur.
            return '1'
        if self.coeff == -1:
            # expressions like "-A(a)"
            return '-'+'*'.join([get_str(arg) for arg in args[1:]])

        # prints expressions like "A(a)", "3*A(a)", "(1+x)*A(a)"
        return '*'.join([get_str(arg) for arg in self.args])

    @property
    def data(self):
        dat = _tensor_data_substitution_dict[self]
        if dat is not None:
            return self.coeff * dat

    def __iter__(self):
        if self.data is None:
            raise ValueError('No iteration on abstract tensors')
        return self.data.flatten().__iter__()


def canon_bp(p):
    """Butler-Portugal canonicalization."""
    if isinstance(p, TensExpr):
        return p.canon_bp()
    return p


def tensor_mul(*a):
    """Product of tensors."""
    if not a:
        return TensMul.from_data(Integer(1), [], [], [])
    t = a[0]
    for tx in a[1:]:
        t = t*tx
    return t


def riemann_cyclic_replace(t_r):
    """Replace Riemann tensor with an equivalent expression.

    ``R(m,n,p,q) -> 2/3*R(m,n,p,q) - 1/3*R(m,q,n,p) + 1/3*R(m,p,n,q)``

    """
    free = sorted(t_r.free, key=lambda x: x[1])
    m, n, p, q = (x[0] for x in free)
    t0 = Rational(2, 3)*t_r
    t1 = - Rational(1, 3)*t_r.substitute_indices((m, m), (n, q), (p, n), (q, p))
    t2 = Rational(1, 3)*t_r.substitute_indices((m, m), (n, p), (p, n), (q, q))
    t3 = t0 + t1 + t2
    return t3


def riemann_cyclic(t2):
    """
    Replace each Riemann tensor with an equivalent expression
    satisfying the cyclic identity.

    This trick is discussed in the reference guide to Cadabra.

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> i, j, k, l = tensor_indices('i j k l', Lorentz)
    >>> R = tensorhead('R', [Lorentz]*4, [[2, 2]])
    >>> t = R(i, j, k, l)*(R(-i, -j, -k, -l) - 2*R(-i, -k, -j, -l))
    >>> riemann_cyclic(t)
    0

    """
    if isinstance(t2, (TensMul, Tensor)):
        args = [t2]
    else:
        args = t2.args
    a1 = [x.split() for x in args]
    a2 = [[riemann_cyclic_replace(tx) for tx in y] for y in a1]
    a3 = [tensor_mul(*v) for v in a2]
    t3 = TensAdd(*a3)
    if not t3:
        return t3
    return canon_bp(t3)


def get_indices(t):
    if not isinstance(t, TensExpr):
        return ()
    return t.get_indices()


def get_tids(t):
    if isinstance(t, TensExpr):
        return t._tids
    return TIDS([], [], [])


def get_coeff(t):
    if isinstance(t, Tensor):
        return Integer(1)
    if isinstance(t, TensMul):
        return t.coeff
    if isinstance(t, TensExpr):
        raise ValueError('no coefficient associated to this tensor expression')
    return t


def contract_metric(t, g):
    if isinstance(t, TensExpr):
        return t.contract_metric(g)
    return t


def perm2tensor(t, g, canon_bp=False):
    """
    Returns the tensor corresponding to the permutation ``g``

    For further details, see the method in ``TIDS`` with the same name.

    """
    if not isinstance(t, TensExpr):
        return t
    new_tids = get_tids(t).perm2tensor(g, canon_bp)
    coeff = get_coeff(t)
    if g[-1] != len(g) - 1:
        coeff = -coeff
    res = TensMul.from_TIDS(coeff, new_tids, is_canon_bp=canon_bp)
    return res


def substitute_indices(t, *index_tuples):
    """
    Return a tensor with free indices substituted according to ``index_tuples``

    ``index_types`` list of tuples ``(old_index, new_index)``

    Note: this method will neither raise or lower the indices, it will just replace their symbol.

    Examples
    ========

    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> i, j, k, l = tensor_indices('i j k l', Lorentz)
    >>> A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
    >>> t = A(i, k)*B(-k, -j)
    >>> t
    A(i, L_0)*B(-L_0, -j)
    >>> t.substitute_indices((i, j), (j, k))
    A(j, L_0)*B(-L_0, -k)

    """
    if not isinstance(t, TensExpr):
        return t
    free = t.free
    free1 = []
    for j, ipos, cpos in free:
        for i, v in index_tuples:
            if i._name == j._name and i._tensortype == j._tensortype:
                if i._is_up == j._is_up:
                    free1.append((v, ipos, cpos))
                else:
                    free1.append((-v, ipos, cpos))
                break
        else:
            free1.append((j, ipos, cpos))

    t = TensMul.from_data(t.coeff, t.components, free1, t.dum)
    return t
