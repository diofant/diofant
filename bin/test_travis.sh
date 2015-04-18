#! /usr/bin/env bash

set -e -x # exit on error and echo each command

if [[ "${TEST_SPHINX}" == "true" ]]; then
    make -C doc html-errors man latex
    LATEXOPTIONS="-interaction=nonstopmode" make -C doc/_build/latex
elif [[ "${TEST_SAGE}" == "true" ]]; then
    sage -v
    sage -python py.test sympy/external/tests/test_sage.py
else
    if [[ "${TEST_DOCTESTS}" == "true" ]]; then
        cat << EOF | python
import sympy
if not sympy.doctest():
    raise Exception('Tests failed')
EOF
        bin/doctest doc/
    elif [[ "${TEST_SLOW}" == "true" ]]; then
        py.test -m 'slow' --duration=100 --split="${SPLIT}" sympy/
    elif [[ "${TEST_EXTRA}" == "true" ]]; then
        if [[ "${TRAVIS_PYTHON_VERSION}" == "2.7" ]]; then
            py.test --duration=100 --cov sympy \
                sympy/printing/tests/test_theanocode.py \
                sympy/external/tests/test_autowrap.py \
                sympy/polys/ sympy/plotting/
        else
            py.test --duration=100 \
                sympy/printing/tests/test_theanocode.py \
                sympy/external/tests/test_autowrap.py \
                sympy/polys/ sympy/plotting/
        fi
    elif [[ "${TRAVIS_PYTHON_VERSION}" == "2.7" ]]; then
        py.test -m 'not slow' --duration=100 --cov sympy --split="${SPLIT}" sympy/
    else
        py.test -m 'not slow' --duration=100 --split="${SPLIT}" sympy/
    fi
fi
