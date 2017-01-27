===========
SymPy 0.6.2
===========

17 Aug 2008

* SymPy is now 50% faster on average (cache:on) and 130% (cache:off) compared to previous release.
* adaptive and faster ``evalf()``
* evalf: numerical summation of hypergeometric series
* evalf: fast and accurate numerical summation
* evalf: oscillatory quadrature
* integrals now support variable transformation
* we can now ``integrate(f(x)⋅diff(f(x),x), x)``
* we can now solve ``a⋅cos(x)=y`` and ``exp(x)+exp(-x)=y``
* printing system refactored
* pprint: new symbol for multiply in unicode mode(x*y -> x⋅y)
* pprint: matrices now look much better
* printing of dicts and sets are now more human-friendly
* latex: now supports sub- and superscripts in symbol names
* ``RootSum.doit()``, now works on all roots
* Wild can now have additional predicates
* numpy-like zeros and ones functions
* ``var('x,y,z')`` now works
* ``((x+y+z)**50).expand()`` is now 4.8x faster
* big assumptions cleanup and rewrite
* access to all object attributes is now ~2.5 times faster
* we try not to let 'is_commutative' to go through (slow) assumptions path
* Add/Mul were optimized (for some cases significantly)
* isympy and sympy.interactive code were merged
* multiple inheritance removed (NoArithMeths, NoRelMeths, RelMeths, ArithMeths are gone)
* ``.nseries()`` is now used as default in ``.series()``
* doctesting was made more robust
