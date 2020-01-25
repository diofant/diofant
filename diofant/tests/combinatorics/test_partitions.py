import random

import pytest

from diofant import Set, default_sort_key, ordered
from diofant.abc import x
from diofant.combinatorics import (IntegerPartition, Partition, RGS_enum,
                                   RGS_rank, RGS_unrank)
from diofant.combinatorics.partitions import random_integer_partition
from diofant.utilities.iterables import partitions


__all__ = ()


def test_partition():
    pytest.raises(ValueError, lambda: Partition(*range(3)))
    pytest.raises(ValueError, lambda: Partition([1, 1, 2]))

    a = Partition([1, 2, 3], [4])
    b = Partition([1, 2], [3, 4])
    c = Partition([x])
    l = [a, b, c]
    l.sort(key=default_sort_key)
    assert l == [c, a, b]
    l.sort(key=lambda w: default_sort_key(w, order='rev-lex'))
    assert l == [c, a, b]

    assert (a == b) is False
    assert a <= b
    assert a < b
    assert (a > b) is False
    assert a != b

    assert (a + 2).partition == [[1, 2], [3, 4]]
    assert (b - 1).partition == [[1, 2, 4], [3]]

    assert (a - 1).partition == [[1, 2, 3, 4]]
    assert (a + 1).partition == [[1, 2, 4], [3]]
    assert (b + 1).partition == [[1, 2], [3], [4]]

    assert a.rank == 1
    assert b.rank == 3

    assert a.RGS == (0, 0, 0, 1)
    assert b.RGS == (0, 0, 1, 1)


def test_integer_partition():
    # no zeros in partition
    pytest.raises(ValueError, lambda: IntegerPartition(list(range(3))))
    # check fails since 1 + 2 != 100
    pytest.raises(ValueError, lambda: IntegerPartition(100, list(range(1, 3))))
    a = IntegerPartition(8, [1, 3, 4])
    assert a.as_dict() == {1: 1, 3: 1, 4: 1}
    b = a.next_lex()
    c = IntegerPartition([1, 3, 4])
    d = IntegerPartition(8, {1: 3, 3: 1, 2: 1})
    assert a == c
    assert a.integer == d.integer
    assert a.conjugate == [3, 2, 2, 1]
    assert (a == b) is False
    assert a <= b
    assert (a > b) is False
    assert a != b

    for i in range(1, 11):
        next = set()
        prev = set()
        a = IntegerPartition([i])
        ans = {IntegerPartition(p) for p in partitions(i)}
        n = len(ans)
        for j in range(n):
            next.add(a)
            a = a.next_lex()
            IntegerPartition(i, a.partition)  # check it by giving i
        for j in range(n):
            prev.add(a)
            a = a.prev_lex()
            IntegerPartition(i, a.partition)  # check it by giving i
        assert next == ans
        assert prev == ans

    assert IntegerPartition([1, 2, 3]).as_ferrers() == '###\n##\n#'
    assert IntegerPartition([1, 1, 3]).as_ferrers('o') == 'ooo\no\no'
    assert str(IntegerPartition([1, 1, 3])) == '[3, 1, 1]'
    assert IntegerPartition([1, 1, 3]).partition == [3, 1, 1]

    pytest.raises(ValueError, lambda: random_integer_partition(-1))
    assert random_integer_partition(1) == [1]
    assert random_integer_partition(10,
                                    seed=[1, 3, 2,
                                          1, 5, 1]) == [5, 2, 1, 1, 1]

    random.seed(0)
    assert random_integer_partition(5, seed=1) == [2, 2, 1]


def test_rgs():
    pytest.raises(ValueError, lambda: RGS_unrank(-1, 3))
    pytest.raises(ValueError, lambda: RGS_unrank(3, 0))
    pytest.raises(ValueError, lambda: RGS_unrank(10, 1))

    pytest.raises(ValueError, lambda: Partition.from_rgs(list(range(3)), list(range(2))))
    pytest.raises(ValueError, lambda: Partition.from_rgs(list(range(1, 3)), list(range(2))))
    assert RGS_enum(-1) == 0
    assert RGS_enum(1) == 1
    assert RGS_unrank(7, 5) == [0, 0, 1, 0, 2]
    assert RGS_unrank(23, 14) == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2]
    assert RGS_rank(RGS_unrank(40, 100)) == 40


def test_sympyissue_9608():
    a = Partition([1, 2, 3], [4])
    b = Partition([1, 2], [3, 4])
    assert list(ordered([a, b], Set._infimum_key))  # does not raise an error
