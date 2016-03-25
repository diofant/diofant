#!/usr/bin/env python
"""
This script uses Setuptools (http://pythonhosted.org/setuptools/), a
collection of enhancements to the standard Python distutils.
"""

import sys
import os
import shutil
from setuptools import setup, Command, find_packages
from setuptools.command.test import test as TestCommand


# Make sure I have the right Python version.
if sys.version_info[:2] < (3, 2):
    print("SymPy requires Python 3.2 or newer. Python %d.%d detected" % sys.version_info[:2])
    sys.exit(-1)


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
        names = ["python-build-stamp-2.4", "MANIFEST", "build", "dist", "doc/_build", ".coverage", ".cache"]

        for f in names:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)

        os.chdir(curr_dir)


class test(TestCommand):
    """Runs all tests."""

    description = "run all tests and doctests"
    user_options = [('cov', None, "gatter coverage information"),
                    ('mark=', "m", "run tests matching given mark expression")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.cov = None
        self.mark = None

    def finalize_options(self):
        TestCommand.finalize_options(self)
        try:
            # New setuptools don't need this anymore
            self.test_args = []
            self.test_suite = True
        except AttributeError:
            pass
        self.pytest_args = []
        if self.cov is not None:
            self.pytest_args.extend(["--cov", "sympy"])
        if self.mark is not None:
            self.pytest_args.extend(["-m", self.mark])

    def run_tests(self):
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


exec(open('sympy/release.py').read())
with open('sympy/__init__.py') as f:
    long_description = f.read().split('"""')[1]

setup(name='sympy',
      version=__version__,
      description='Computer algebra system (CAS) in Python',
      long_description=long_description,
      maintainer='Sergey B Kirpichev',
      license='BSD',
      keywords="Math CAS",
      url='http://omg.rtfd.org',
      packages=find_packages(),
      ext_modules=[],
      cmdclass={'test': test,
                'clean': clean},
      classifiers=[
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.2',
          'Programming Language :: Python :: 3.4',
      ],
      tests_require=['pytest>=2.7.0', 'flake8', 'pep8-naming', 'pytest-cov'],
      install_requires=['mpmath>=0.19', 'strategies>=0.2.3', 'cachetools'],
      setup_requires=['setuptools>=5.5.1,<=19.4', 'pip>=6.0'],
      extras_require={
          'exports': ["numpy", "scipy", "Theano"],
          'gmpy': ["gmpy>=1.16"],
          'plot': ["matplotlib"],
          'interactive': ["ipython>=2.3.0"],
      }
)
