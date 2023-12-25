"""Combinatorics package."""

from .generators import alternating, cyclic, dihedral, symmetric
from .graycode import GrayCode
from .group_constructs import DirectProduct
from .named_groups import (AbelianGroup, AlternatingGroup, CyclicGroup,
                           DihedralGroup, RubikGroup, SymmetricGroup)
from .partitions import (IntegerPartition, Partition, RGS_enum, RGS_rank,
                         RGS_unrank)
from .perm_groups import PermutationGroup
from .permutations import Cycle, Permutation
from .polyhedron import (Polyhedron, cube, dodecahedron, icosahedron,
                         octahedron, tetrahedron)
from .prufer import Prufer
from .subsets import Subset


__all__ = ('alternating', 'cyclic', 'dihedral', 'symmetric', 'GrayCode',
           'DirectProduct', 'AbelianGroup', 'AlternatingGroup', 'CyclicGroup',
           'DihedralGroup', 'RubikGroup', 'SymmetricGroup', 'IntegerPartition',
           'Partition', 'RGS_enum', 'RGS_rank', 'RGS_unrank',
           'PermutationGroup', 'Cycle', 'Permutation', 'Polyhedron', 'cube',
           'dodecahedron', 'icosahedron', 'octahedron', 'tetrahedron',
           'Prufer', 'Subset')
