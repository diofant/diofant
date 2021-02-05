"""Singleton mechanism"""

from __future__ import annotations

import typing

from .assumptions import ManagedProperties


class SingletonRegistry:
    """
    The registry for the singleton classes.

    Several classes in Diofant appear so often that they are
    singletonized, that is, using some metaprogramming they are made
    so that they can only be instantiated once (see the
    :class:`diofant.core.singleton.Singleton` class for details).  For
    instance, every time you create ``Integer(0)``, this will return
    the same instance, :class:`diofant.core.numbers.Zero`.

    >>> a = Integer(0)
    >>> a is S.Zero
    True

    All singleton instances are attributes of the
    :data:`~diofant.core.singleton.S` object, so ``Integer(0)`` can
    also be accessed as ``S.Zero``.

    Notes
    =====

    For the most part, the fact that certain objects are singletonized
    is an implementation detail that users shouldn't need to worry
    about.  In Diofant library code, :keyword:`is` comparison is often
    used for performance purposes.  The primary advantage of
    :data:`~diofant.core.singleton.S` for end users is the convenient
    access to certain instances that are otherwise difficult to type,
    like ``S.Half`` (instead of ``Rational(1, 2)``).

    When using ``is`` comparison, make sure the argument is a
    :class:`~diofant.core.basic.Basic` instance.  For example,

    >>> int(0) is S.Zero
    False

    """

    def __init__(self):
        self._classes_to_install = {}
        # Dict of classes that have been registered, but that have not have been
        # installed as an attribute of this SingletonRegistry.
        # Installation automatically happens at the first attempt to access the
        # attribute.
        # The purpose of this is to allow registration during class
        # initialization during import, but not trigger object creation until
        # actual use (which should not happen until after all imports are
        # finished).

    def register(self, cls):
        self._classes_to_install[cls.__name__] = cls

    def __setattr__(self, name, obj):
        setattr(self.__class__, name, obj)

    def __getattr__(self, name):
        """Python calls __getattr__ if no attribute of that name was installed
        yet.

        This __getattr__ checks whether a class with the requested name was
        already registered but not installed; if no, raises an AttributeError.
        Otherwise, retrieves the class, calculates its singleton value, installs
        it as an attribute of the given name, and unregisters the class.

        """
        if name not in self._classes_to_install:
            raise AttributeError(
                f"Attribute '{name}' was not installed on Diofant registry {self}")
        class_to_install = self._classes_to_install[name]
        value_to_install = class_to_install()
        self.__setattr__(name, value_to_install)
        del self._classes_to_install[name]
        return value_to_install

    def __repr__(self):
        return 'S'


#: Alias for instance of :class:`SingletonRegistry`.
S: SingletonRegistry = SingletonRegistry()


class Singleton(type):
    """
    Metaclass for singleton classes.

    A singleton class has only one instance which is returned every
    time the class is instantiated. Additionally, this instance can be
    accessed through the global registry object
    :data:`~diofant.core.singleton.S` as ``S.<class_name>``.

    Examples
    ========

    >>> class MySingleton(Basic, metaclass=Singleton):
    ...     pass
    >>> Basic() is Basic()
    False
    >>> MySingleton() is MySingleton()
    True
    >>> S.MySingleton is MySingleton()
    True

    Notes
    =====

    Instance creation is delayed until the first time the value is accessed.

    This metaclass is a subclass of ManagedProperties because that is the
    metaclass of many classes that need to be Singletons (Python does not allow
    subclasses to have a different metaclass than the superclass, except the
    subclass may use a subclassed metaclass).

    """

    _instances: dict[type[typing.Any], typing.Any] = {}
    'Maps singleton classes to their instances.'

    def __new__(cls, *args, **kwargs):
        result = super().__new__(cls, *args, **kwargs)
        S.register(result)
        return result

    def __call__(cls, *args, **kwargs):
        # Called when application code says SomeClass(), where SomeClass is a
        # class of which Singleton is the metaclass.
        # __call__ is invoked first, before __new__() and __init__().
        if cls not in cls.__class__._instances:
            # Invokes the standard constructor of SomeClass.
            cls.__class__._instances[cls] = super().__call__(*args, **kwargs)
        return cls.__class__._instances[cls]


class SingletonWithManagedProperties(Singleton, ManagedProperties):
    """Metaclass for singleton classes with managed properties."""
