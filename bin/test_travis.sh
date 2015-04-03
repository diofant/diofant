#! /usr/bin/env bash

set -e -x # exit on error and echo each command

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
    export LANG=C
    py.test -k 'print' sympy/
    bin/doctest
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
        py.test sympy/printing/tests/test_theanocode.py
        py.test sympy/external/tests/test_autowrap.py
        py.test --duration=100 sympy/polys/ sympy/plotting/
    else
        py.test -m 'not slow' --duration=100 --split="${SPLIT}" sympy/
    fi
fi
