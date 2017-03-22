===========
SymPy 0.6.5
===========

17 Jul 2009

* Geometric Algebra Improvements

  - Upgrade GA module with left and right contraction operations
  - Add intersection test for the vertical segment, reimplementation of convex_hull

* Implement ``series()`` as function
* Core improvements

  - Refactor ``Number.eval_power``
  - fix bugs in ``Number.eval_power``

* Matrix improvements:

  - Improve jacobian function, introduce vec and vech

* Solver improvements:

  - solutions past linear factor found in tsolve
  - refactor sympy.solvers.guess_solve_strategy
  - Small cleanups to the ODE solver and tests
  - Fix corner case for Bernoulli equation

* Improvements on partial differential equations solvers

  - Added separation of variables for PDEs

* Expand improvements

  - Refactoring
  - ``exp(x) exp(y)`` is no longer automatically combined into ``exp(x+y)``, use ``powsimp`` for that

* Documentation improvements:

  - Test also documentation under doc/
  - Added many docstrings
  - Fix Sphinx complaints/warnings/errors
  - Doctest coverage

* New logic module

  - Efficient DPLL algorithm

* LaTeX printer improvements:

  - Handle standard form in the LaTeX printer correctly
  - Latex: print_Mul fix (:sympyissue:`4381`)
  - Robust printing of latex sub and superscripts
  - sorting print_Add output using a main variable
  - Matrix printing improvements

* MathML printing improvements:

  - MathML's printer extended

* Testing framework improvements

  - Make tests pass without the "py" module

* Polynomial module improvements:

  - Fixed subresultant PRS computation and ``ratint()``
  - Removed old module ``sympy.polynomials``

* limit fixes:

  - Compute the finite parts of the limit of a sum by direct substitution

* Test coverage script
* Code quality improvements (remove string exceptions, code quality test improvements)
* C code generation
* Update mpmath
