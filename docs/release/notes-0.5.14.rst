============
SymPy 0.5.14
============

26 Apr 2008

* SymPy is now 25% faster on average compared to the previous release

  - ``__eq__``/``__ne__``/``__nonzero__`` returns ``True``/``False`` directly so dict lookups are not expensive anymore
  - ``sum(x**i/i,i=1..400)`` is now 4.8x faster
  - ``isinstance(term, C.Mul)`` was replaced by ``term.is_Mul`` and similarly for other basic classes

* Documentation was improved a lot. See https://docs.sympy.org/
* rsolve_poly & rsolve_hyper fixed
* ``subs`` and ``subs_dict`` unified to ``.subs()``
* faster and more robust polynomials module
* improved ``Matrix.det()``, implemented Berkowitz algorithm
* improved isympy (interactive shell for SymPy)
* pretty-printing improved
* ``Rel``, ``Eq``, ``Ne``, ``Lt``, ``Le``, ``Gt``, ``Ge`` implemented
* ``Limit`` class represents unevaluated limits now
* Bailey-Borwein-Plouffe algorithm (finds the nth hexidecimal digit of pi without calculating the previous digits) implemented
* solver for transcendental equations added
* ``.nseries()`` methods implemented (more robust/faster than .oseries)
* multivariate Lambdas implemented
