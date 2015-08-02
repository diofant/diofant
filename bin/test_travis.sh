#! /usr/bin/env bash

set -e -x # exit on error and echo each command

if [[ "${TEST_SLOW}" == "true" ]]; then
    py.test -m 'slow' --duration=100 --split="${SPLIT}" \
        --ignore sympy/utilities/autowrap.py \
        --ignore sympy/plotting/plot.py \
        --ignore sympy/plotting/plot_implicit.py sympy/
elif [[ "${TEST_EXTRA}" == "true" ]]; then
    if [[ "${TRAVIS_PYTHON_VERSION}" == "2.7" ]]; then
        py.test --duration=100 --cov sympy --doctest-modules \
            sympy/printing/tests/test_theanocode.py \
            sympy/external/tests/test_autowrap.py \
            sympy/polys/ sympy/plotting/ docs \
            --ignore docs/tutorial/gotchas.rst # XXX: workaround __future__ imports!
    else
        py.test --duration=100 --doctest-modules \
            sympy/printing/tests/test_theanocode.py \
            sympy/external/tests/test_autowrap.py \
            sympy/polys/ sympy/plotting/ docs
    fi
    make -C docs html-errors man latex
    LATEXOPTIONS="-interaction=nonstopmode" make -C docs/_build/latex
    python examples/all.py -q
elif [[ "${TRAVIS_PYTHON_VERSION}" == "2.7" ]]; then
    py.test -m 'not slow' --duration=100 --cov sympy --split="${SPLIT}" \
        --ignore sympy/utilities/autowrap.py \
        --ignore sympy/plotting/plot.py \
        --ignore sympy/plotting/plot_implicit.py \
        --doctest-modules sympy/
else
    if [[ "${TRAVIS_PYTHON_VERSION}" == "pypy3" ]]; then
        EXTRA_IGNORE="--ignore sympy/matrices/dense.py --ignore sympy/tensor/tensor.py --ignore sympy/utilities/lambdify.py"
    else
        EXTRA_IGNORE=""
    fi
    py.test -m 'not slow' --duration=100 --split="${SPLIT}" \
        --ignore sympy/utilities/autowrap.py \
        --ignore sympy/plotting/plot.py \
        --ignore sympy/plotting/plot_implicit.py \
        ${EXTRA_IGNORE} --doctest-modules sympy/
fi
