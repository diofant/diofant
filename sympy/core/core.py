""" The core's core. """


class Registry(object):
    """
    Base class for registry objects.

    Registries map a name to an object using attribute notation. Registry
    classes behave singletonically: all their instances share the same state,
    which is stored in the class object.

    All subclasses should set `__slots__ = []`.
    """
    __slots__ = []

    def __setattr__(self, name, obj):
        setattr(self.__class__, name, obj)

    def __delattr__(self, name):
        delattr(self.__class__, name)

#: A set containing all sympy class objects
all_classes = set()


class BasicMeta(type):
    def __init__(cls, *args, **kwargs):
        all_classes.add(cls)
