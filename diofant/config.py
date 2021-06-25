"""Configuration utilities."""

import ast
import contextlib
import os


__all__ = 'setup',


_default_config = {
    'USE_COLLINS_RESULTANT':      False,
    'USE_HEU_GCD':                True,
    'HEU_GCD_MAX':                6,

    'FALLBACK_GCD_ZZ_METHOD':     'prs',
    'GCD_AA_METHOD':              'prs',

    'USE_IRREDUCIBLE_IN_FACTOR':  False,
    'USE_CYCLOTOMIC_FACTOR':      True,

    'EEZ_RESTART_IF_NEEDED':      True,
    'EEZ_NUMBER_OF_CONFIGS':      3,
    'EEZ_NUMBER_OF_TRIES':        5,
    'EEZ_MODULUS_STEP':           2,

    'GF_IRRED_METHOD':            'rabin',
    'GF_FACTOR_METHOD':           'zassenhaus',

    'AA_FACTOR_METHOD':           'modular',

    'GROEBNER':                   'buchberger',
    'MINPOLY_METHOD':             'compose',

    'KARATSUBA_CUTOFF':           100,
}

_current_config = {}


@contextlib.contextmanager
def using(**kwargs):
    for k, v in kwargs.items():
        setup(k, v)

    yield

    for k in kwargs:
        setup(k)


def setup(key, value=None):
    """Assign a value to (or reset) a configuration item."""
    key = key.upper()

    if value is not None:
        _current_config[key] = value
    else:
        _current_config[key] = _default_config[key]


def query(key):
    """Ask for a value of the given configuration item."""
    return _current_config.get(key.upper(), None)


def configure():
    """Initialized configuration of polys module."""
    for key, default in _default_config.items():
        _current_config[key] = default

        value = os.getenv('DIOFANT_' + key)
        if value is not None:
            try:
                value = ast.literal_eval(value)
            except (SyntaxError, ValueError):
                pass
            else:
                if type(value) is type(default):
                    _current_config[key] = value


configure()
