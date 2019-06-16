from .ndim_array import NDimArray


class MutableNDimArray(NDimArray):

    def _subs(self, old, new, **hints):
        return self
