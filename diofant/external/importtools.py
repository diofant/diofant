"""Tools to assist importing optional external modules."""

import sys
import warnings
from distutils.version import LooseVersion

from .. import DIOFANT_DEBUG


def import_module(module, min_module_version=None, min_python_version=None,
                  warn_not_installed=False, warn_old_version=True,
                  module_version_attr='__version__',
                  module_version_attr_call_args=None,
                  __import__kwargs={}, catch=()):
    """
    Import and return a module if it is installed.

    If the module is not installed, it returns None.

    A minimum version for the module can be given as the keyword argument
    min_module_version.  This should be comparable against the module version.
    By default, module.__version__ is used to get the module version.  To
    override this, set the module_version_attr keyword argument.  If the
    attribute of the module to get the version should be called (e.g.,
    module.version()), then set module_version_attr_call_args to the args such
    that module.module_version_attr(*module_version_attr_call_args) returns the
    module's version.

    If the module version is less than min_module_version using the Python <
    comparison, None will be returned, even if the module is installed. You can
    use this to keep from importing an incompatible older version of a module.

    You can also specify a minimum Python version by using the
    min_python_version keyword argument.  This should be comparable against
    sys.version_info.

    If the keyword argument warn_not_installed is set to True (or
    DIOFANT_DEBUG is True), the function will emit a UserWarning when
    the module is not installed.

    By default, or if the keyword argument warn_old_version is set to True (or
    DIOFANT_DEBUG is True), the function will emit a UserWarning when
    the library is installed, but cannot be imported because of the
    min_module_version or min_python_version options.

    This function uses __import__() to import the module.  To pass additional
    options to __import__(), use the __import__kwargs keyword argument.  For
    example, to import a submodule A.B, you must pass a nonempty fromlist option
    to __import__.

    This catches ImportError to determine if the module is not installed.  To
    catch additional errors, pass them as a tuple to the catch keyword
    argument.

    Examples
    ========

    >>> numpy = import_module('numpy')

    >>> numpy = import_module('numpy', min_python_version=(2, 7),
    ...                       warn_old_version=False)

    numpy.__version__ is a string

    >>> numpy = import_module('numpy', min_module_version='1.5',
    ...                       warn_old_version=False)

    gmpy does not have __version__, but it does have gmpy.version()

    >>> gmpy = import_module('gmpy2', min_module_version='2.0.0',
    ...                      module_version_attr='version',
    ...                      module_version_attr_call_args=(),
    ...                      warn_old_version=False)

    To import a submodule, you must pass a nonempty fromlist to
    __import__().  The values do not matter.

    >>> p3 = import_module('mpl_toolkits.mplot3d',
    ...                     __import__kwargs={'fromlist':['something']})

    matplotlib.pyplot can raise RuntimeError when the display cannot be opened

    >>> matplotlib = import_module('matplotlib',
    ...                            __import__kwargs={'fromlist':['pyplot']},
    ...                            catch=(RuntimeError,))

    See Also
    ========

    __import__

    """
    warn_old_version = True if DIOFANT_DEBUG else warn_old_version
    warn_not_installed = True if DIOFANT_DEBUG else warn_not_installed

    # Check Python first so we don't waste time importing a module we can't use
    if min_python_version:
        if sys.version_info < min_python_version:
            min_python_version = '.'.join(map(str, min_python_version))
            if warn_old_version:
                warnings.warn("Python version is too old to use %s "
                              "(%s or newer required)" % (module,
                                                          min_python_version),
                              UserWarning)
            return

    try:
        mod = __import__(module, **__import__kwargs)
    except ImportError:
        if warn_not_installed:
            warnings.warn("%s module is not installed" % module, UserWarning)
        return

    if min_module_version:
        modversion = getattr(mod, module_version_attr)
        if module_version_attr_call_args is not None:
            modversion = modversion(*module_version_attr_call_args)
        if LooseVersion(modversion) < LooseVersion(min_module_version):
            if warn_old_version:
                warnings.warn("%s version is too old to use "
                              "(%s or newer required)" % (module,
                                                          min_module_version),
                              UserWarning)
            return

    return mod
