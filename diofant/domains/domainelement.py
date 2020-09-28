"""Trait for implementing domain elements."""

import abc

from ..printing.defaults import DefaultPrinting


class DomainElement(abc.ABC, DefaultPrinting):
    """
    Represents an element of a domain.

    Mix in this trait into a class which instances should be recognized as
    elements of a domain.  Property ``parent`` gives that domain.
    """

    @property
    @abc.abstractmethod
    def parent(self):
        raise NotImplementedError
