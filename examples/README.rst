This directory contains Diofant example programs.

-------------------
DIRECTORY STRUCTURE
-------------------

The examples are broken up into three categories based on difficulty of both
the mathematics and programming concepts.  They roughly follow the following
guide:

intermediate :
  Demonstrations of more complex mathematical concepts, but still for
  someone with little programming experience.

advanced :
  Larger demonstrations of advanced mathematical topics.

----------------
RUNNING EXAMPLES
----------------

All the working examples can be run by executing the "all.py" script, use
``./all.py -h`` for usage, if an example is known to be broken it will be
commented out in this script.

To run the individual examples one needs to have Python version >= 3.4
installed and Diofant must be in your PYTHONPATH environment variable.  Most
examples can be run from the command line python and the name of the example::

    $ python examples/intermediate/trees.py
    -258121*x**20 - 205634*x**19 - 123440*x**18 - 67402*x**17 - 34663*x**16 - 17711*x**15 - 8878*x**14 - 4577*x**13 - 2373*x**12 - 1607*x**11 + 106*x**10 + 47*x**9 + 23*x**8 + 11*x**7 + 6*x**6 + 3*x**5 + 2*x**4 + x**3 + x**2 + x + 1
    [Integer(1), Integer(1), Integer(1), Integer(1), Integer(2), Integer(3), Integer(6), Integer(11), Integer(23), Integer(47), Integer(106)]

Note, that on most systems, the current directory is searched by Python
automatically, so "python examples/intermediate/trees.py" works from the diofant root
directory, however there are systems (Ubuntu Intrepid) where this
doesn't work by default, unless you put "PYTHONPATH=." into your
.bashrc for example.
