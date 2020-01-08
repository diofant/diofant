from .ndim_array import NDimArray


class MutableNDimArray(NDimArray):
    """A mutable version of the N-dim array."""

    def _subs(self, old, new, **hints):
        return self
