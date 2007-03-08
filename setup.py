#!/usr/bin/env python
"""Distutils based setup script for Sympy.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other useful commands

    python setup.py test  -> will run the complete test suite
    python setup.py test_core -> will run only tests concerning core features
    
Or, if all else fails, feel free to write to the sympy list at
sympy@googlegroups.com and ask for help.
"""

from distutils.core import setup
from distutils.core import Command
import sys

# Make sure I have the right Python version.
if sys.version_info[2] < 4:
    print "Sympy requires Python 2.4 or newer.  Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)

class test_sympy(Command):
    """Runs all tests under the tests/ folder
    """
    
    description = "Automatically run the test suite for Sympy."
    user_options = []  # distutils complains if this is not here.


    def initialize_options(self):  # distutils wants this
        pass
    
    def finalize_options(self):    # this too
        pass
    
    def run(self):
        try:
            import py
        except ImportError:
            print """In order to run the tests, you need codespeak's py.lib
            web page: http://codespeak.net/py/dist/
            If you are on debian systems, the package is named python-codespeak-lib
            """
            sys.exit(-1)
        py.test.cmdline.main(args=["tests"])
        # change back to the current directory

class test_sympy_core(Command):
    """Run only the tests concerning features of sympy.core.
    It's a lot faster than running the complete test suite.
    """
    
    description = "Automatically run the core test suite for Sympy."
    user_options = []  # distutils complains if this is not here.

    tests_to_run = ["tests/test_arit.py", "tests/test_basic.py", 
                   "tests/test_diff.py", "tests/test_equal.py", 
                   "tests/test_eval.py", "tests/test_evalf.py", 
                   "tests/test_functions.py", "tests/test_hashing.py", 
                   "tests/test_numbers.py", "tests/test_series.py", 
                   "tests/test_str.py", "tests/test_subs.py", 
                   "tests/test_symbol.py", "tests/test_util.py" 
                   ]
    
    def initialize_options(self):  # distutils wants this
        pass
    
    def finalize_options(self):    # this too
        pass
    
    def run(self):
        try:
            import py
        except ImportError:
            print """In order to run the tests, you need codespeak's py.lib
            web page: http://codespeak.net/py/dist/
            If you are on debian systems, the package is named python-codespeak-lib
            """
            sys.exit(-1)
        #for x in tests_ro_run
        py.test.cmdline.main(args=self.tests_to_run)


setup(
      name='Sympy', 
      version='1.0-pre', 
      description='Computer algebra system (CAS) in Python', 
      url='http://code.google.com/p/sympy', 
      packages=['sympy', 'sympy.core', 'sympy.modules'],
      ext_modules = [],
      cmdclass    = {'test': test_sympy, 
                     'test_core' : test_sympy_core,
                     },
      )

