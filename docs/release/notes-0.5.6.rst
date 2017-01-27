===========
SymPy 0.5.6
===========

30 Oct 2007

* ``_sage_()`` methods implemented to convert any SymPy expression to a SAGE expression.
* ``isympy`` fixed so that it always tries the local unpacked sympy first (the one in the directory where isympy sits) and only then the system wide installation of sympy (Debian package for example).
