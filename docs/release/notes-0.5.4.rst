===========
SymPy 0.5.4
===========

5 Oct 2007

* ``Log`` and ``ApplyLog`` classes were simplified to ``log``, as was in the 0.4.3 version (the same for all other classes, like ``sin`` or ``cos``).
* Limits algorithm was fixed and it works very reliably (there are some bugs in the series facility though that make some limits fail), see `this post <https://groups.google.com/forum/#!topic/sympy/mKBEvVrFN8o>`_ for more details.
* All functions arguments are now accessed using the ``sin(x)[:]`` idiom again, as in the 0.4.3 version (instead of the old ``sin(x)._args`` or ``sin(x).args`` which was briefly introduced in the 0.5.x series).
