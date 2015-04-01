#! /usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_SPHINX}" == "true" ]]; then
    cd doc
    make html-errors
    make man
    make latex
    cd _build/latex
    export LATEXOPTIONS="-interaction=nonstopmode"
    make all
elif [[ "${TEST_SAGE}" == "true" ]]; then
    sage -v
    sage -python bin/test sympy/external/tests/test_sage.py
elif [[ "${TEST_ASCII}" == "true" ]]; then
    export LANG=c
    mkdir empty
    cd empty
    cat <<EOF | python
import sympy
sympy.test('print')
EOF
    cd ..
    bin/doctest
else
    # We change directories to make sure that we test the installed version of
    # sympy.
    mkdir empty
    cd empty

    if [[ "${TEST_DOCTESTS}" == "true" ]]; then
        cat << EOF | python
import sympy
if not sympy.doctest():
    raise Exception('Tests failed')
EOF
        cd ..
        bin/doctest doc/
    elif [[ "${TEST_SLOW}" == "true" ]]; then
        cd ..
        py.test -m 'slow' --duration=100 --split="${SPLIT}" sympy/
    elif [[ "${TEST_THEANO}" == "true" ]]; then
        cat << EOF | python
import sympy
if not sympy.test('*theano*'):
    raise Exception('Tests failed')
EOF
    elif [[ "${TEST_GMPY}" == "true" ]] && [[ "${TEST_MATPLOTLIB}" == "true" ]]; then
        cd ..
        py.test --duration=100 sympy/polys/ sympy/plotting/
    elif [[ "${TEST_AUTOWRAP}" == "true" ]]; then
        cd ..
        py.test sympy/external/tests/test_autowrap.py
    else
        cd ..
        py.test -m 'not slow' --duration=100 --split="${SPLIT}" sympy/
    fi
fi
