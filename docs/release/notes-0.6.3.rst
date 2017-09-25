===========
SymPy 0.6.3
===========

19 Nov 2008

* port to python2.6 (all tests pass)
* port to jython (all tests pass except those depending on the "ast" module)
* true division fixed (all tests pass with "-Qnew" Python option)
* ``http://buildbot.sympy.org`` created, sympy is now regularly tested on python2.4, 2.5, 2.6 on both i386 and amd64 architectures.
* py.bench -- py.test based benchmarking added
* bin/test -- simple py.test like testing framework, without external dependencies, nice colored output
* most limits now work
* factorization over `Z[x]` greatly improved
* Piecewise function added
* ``nsimplify()`` implemented
* ``symbols`` and ``var`` syntax unified
* C code printing
* many bugfixes
