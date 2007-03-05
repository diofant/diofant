#!/usr/bin/env python

#from distutils.command.build_py import build_py
from distutils.core import setup
from distutils.core import Command
from distutils.extension import Extension
import os, sys


class test_sympy(Command):
    description = "Automatically run the test suite for Biopython."
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
            return
        py.test.cmdline.main(args=["tests"])
        # change back to the current directory

setup(
      name='Sympy', 
      version='1.0-pre', 
      description='Computer algebra system (CAS) in Python', 
      url='http://code.google.com/p/sympy', 
      packages=['sym', 'sym.core', 'sym.modules'],
      ext_modules = [],
      cmdclass    = {'test': test_sympy},
      )

