"""
Combinatorics package.
"""

from .permutations import Permutation, Cycle  # noqa: F401
from .prufer import Prufer  # noqa: F401
from .generators import cyclic, alternating, symmetric, dihedral  # noqa: F401
from .subsets import Subset  # noqa: F401
from .partitions import (Partition, IntegerPartition, RGS_rank,  # noqa: F401
                         RGS_unrank, RGS_enum)
from .polyhedron import (Polyhedron, tetrahedron, cube, octahedron,  # noqa: F401
                         dodecahedron, icosahedron)
from .perm_groups import PermutationGroup  # noqa: F401
from .group_constructs import DirectProduct  # noqa: F401
from .graycode import GrayCode  # noqa: F401
from .named_groups import (SymmetricGroup, DihedralGroup, CyclicGroup,  # noqa: F401
                           AlternatingGroup, AbelianGroup, RubikGroup)
