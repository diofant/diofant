"""
Combinatorics package.
"""

from .permutations import Permutation, Cycle
from .prufer import Prufer
from .generators import cyclic, alternating, symmetric, dihedral
from .subsets import Subset
from .partitions import (Partition, IntegerPartition, RGS_rank,
                         RGS_unrank, RGS_enum)
from .polyhedron import (Polyhedron, tetrahedron, cube, octahedron,
                         dodecahedron, icosahedron)
from .perm_groups import PermutationGroup
from .group_constructs import DirectProduct
from .graycode import GrayCode
from .named_groups import (SymmetricGroup, DihedralGroup, CyclicGroup,
                           AlternatingGroup, AbelianGroup, RubikGroup)
