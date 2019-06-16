"""Implementation of :class:`CompositeDomain` class. """

from ..polys.polyerrors import GeneratorsError
from .domain import Domain


__all__ = 'CompositeDomain',


class CompositeDomain(Domain):
    """Base class for composite domains, e.g. ZZ[x], ZZ(X). """

    is_Composite = True

    gens, ngens, symbols, domain = [None]*4

    def inject(self, *symbols):
        """Inject generators into this domain. """
        if not (set(self.symbols) & set(symbols)):
            return self.__class__(self.domain, self.symbols + symbols, self.order)
        else:
            raise GeneratorsError("common generators in %s and %s" % (self.symbols, symbols))
