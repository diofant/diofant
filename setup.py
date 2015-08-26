#!/usr/bin/env python
"""
Setuptools-based setup script for SymPy.

This uses Setuptools (http://pythonhosted.org/setuptools/setuptools.html),
a collection of enhancements to the standard Python distutils.
For the easiest installation just type the command (you'll probably
need root privileges for that):

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other commands:

    python setup.py clean -> will clean all trash (*.pyc and stuff)
    python setup.py test  -> will run the complete test suite
    python setup.py audit -> will run pyflakes checker on source code

To get a full list of avaiable commands, read the output of:

    python setup.py --help-commands

Or, if all else fails, feel free to write to the sympy list at
sympy@googlegroups.com and ask for help.
"""

import sys
import subprocess
import os
import shutil
import glob
from setuptools import setup, Command, find_packages
from setuptools.command.test import test as TestCommand


# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 7):
    print("SymPy requires Python 2.7 or newer. Python %d.%d detected" % sys.version_info[:2])
    sys.exit(-1)


class audit(Command):
    """Audits SymPy's source code for following issues:
        - Names which are used but not defined or used before they are defined.
        - Names which are redefined without having been used.
    """

    description = "Audit SymPy source with PyFlakes"
    user_options = []

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        import os
        try:
            import pyflakes.scripts.pyflakes as flakes
        except ImportError:
            print("In order to run the audit, you need to have PyFlakes installed.")
            sys.exit(-1)
        dirs = (os.path.join(*d) for d in (m.split('.') for m in find_packages()))
        warns = 0
        for dir in dirs:
            for filename in os.listdir(dir):
                if filename.endswith('.py') and filename != '__init__.py':
                    warns += flakes.checkPath(os.path.join(dir, filename))
        if warns > 0:
            print("Audit finished with total %d warnings" % warns)


class clean(Command):
    """Cleans *.pyc and debian trashs, so you should get the same copy as
    is in the VCS.
    """

    description = "remove build files"
    user_options = [("all", "a", "the same")]

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        dir_setup = os.path.dirname(os.path.realpath(__file__))
        curr_dir = os.getcwd()
        for root, dirs, files in os.walk(dir_setup):
            for file in files:
                if file.endswith('.pyc') and os.path.isfile:
                    os.remove(os.path.join(root, file))

        os.chdir(dir_setup)
        names = ["python-build-stamp-2.4", "MANIFEST", "build", "dist", "doc/_build", ".coverage"]

        for f in names:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)

        os.chdir(curr_dir)


class test_sympy(TestCommand):
    """Runs all tests."""

    description = "run all tests and doctests"

    def finalize_options(self):
        TestCommand.finalize_options(self)
        _test_args = [
            '--ignore=setup.py',
            '--verbose',
            '--durations=100',
            '--doctest-modules',
        ]
        extra_args = os.environ.get('PYTEST_EXTRA_ARGS')
        if extra_args is not None:
            _test_args.extend(extra_args.split())
        self.test_args = _test_args
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


exec(open('sympy/release.py').read())
with open('sympy/__init__.py') as f:
    long_description = f.read().split('"""')[1]

setup(name='sympy',
      version=__version__,
      description='Computer algebra system (CAS) in Python',
      long_description=long_description,
      maintainer='Sergey B Kirpichev',
      maintainer_email='skirpichev@gmail.com',
      license='BSD',
      keywords="Math CAS",
      url='http://omg.rtfd.org',
      packages=find_packages(),
      ext_modules=[],
      package_data={
          'sympy.logic.benchmarks': ['input/*.cnf'],
      },
      cmdclass={'test': test_sympy,
                'clean': clean,
                'audit': audit},
      classifiers=[
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
      ],
      tests_require=['pytest'],
      install_requires=['mpmath>=0.19', 'decorator', 'strategies>=0.2.3']
      )
