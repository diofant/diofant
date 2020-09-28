"""Implementation of :class:`CompositeDomain` class."""

from ..polys.polyerrors import GeneratorsError
from .domain import Domain


class CompositeDomain(Domain):
    """Base class for composite domains, e.g. ZZ[x], ZZ(X)."""

    gens, ngens, symbols, domain = [None]*4

    def inject(self, *symbols, front=False):
        """Inject generators into this domain."""
        if not (set(self.symbols) & set(symbols)):
            if front:
                symbols += self.symbols
            else:
                symbols = self.symbols + symbols

            return self.__class__(self.domain, symbols, self.order)
        else:
            raise GeneratorsError(f'common generators in {self.symbols} and {symbols}')
