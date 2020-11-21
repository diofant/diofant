"""Implementation of :class:`SimpleDomain` class."""

from .domain import Domain


class SimpleDomain(Domain):
    """Base class for simple domains, e.g. ZZ, QQ."""

    def inject(self, *gens):
        """Inject generators into this domain."""
        return self.poly_ring(*gens)
