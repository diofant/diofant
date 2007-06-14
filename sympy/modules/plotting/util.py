def getkwarg(kwargs, kw, default):
    if kw in kwargs:
        return kwargs[kw]
    return default