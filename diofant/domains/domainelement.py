"""Trait for implementing domain elements. """

from ..printing.defaults import DefaultPrinting


__all__ = 'DomainElement',


class DomainElement(DefaultPrinting):
    """
    Represents an element of a domain.

    Mix in this trait into a class which instances should be recognized as
    elements of a domain.  Property ``parent`` gives that domain.
    """

    @property
    def parent(self):
        raise TypeError
