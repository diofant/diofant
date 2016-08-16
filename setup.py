#!/usr/bin/env python3
"""
This script uses Setuptools (http://pythonhosted.org/setuptools/), a
collection of enhancements to the standard Python distutils.
"""

import re
import sys
import os
import shutil
from setuptools import setup, Command, find_packages


# Make sure I have the right Python version.  We can drop this
# when setuptools 24.2.1 enter into the Debian stable.
if sys.version_info[:2] < (3, 4):
    print('Diofant requires Python 3.4 or newer. '
          'Python %d.%d detected' % sys.version_info[:2])
    sys.exit(-1)


class clean(Command):
    """Cleans *.pyc and debian trashs, so you should get the same copy as
    is in the VCS.
    """

    description = 'remove build files'
    user_options = [('all', 'a', 'the same')]

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
        names = ['python-build-stamp-2.4', 'MANIFEST', 'build', 'dist',
                 'doc/_build', '.coverage', '.cache']

        for f in names:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)

        os.chdir(curr_dir)


with open('diofant/__init__.py') as f:
    source = f.read()
    long_description = source.split('"""')[1]

    for line in source.splitlines():
        m = re.match('^__version__ = "([0-9a.]+)".*$', line)
        if m:
            __version__ = m.group(1)

setup_reqs = ['setuptools>=5.5.1', 'pip>=6.0', 'pytest-runner']
extra_reqs = {'exports': ['numpy', 'scipy', 'Theano'],
              'gmpy': ['gmpy2'],
              'plot': ['pyparsing!=2.1.2', 'matplotlib'],
              'interactive': ['ipython>=2.3.0'],
              'docs': ['sphinx>=1.2.3', 'numpydoc', 'sphinx_rtd_theme'],
              }
extra_reqs['develop'] = ['pytest>=3.0', 'flake8>=2.5.5', 'pep8-naming',
                         'pytest-cov', 'coverage'] + setup_reqs

setup(name='Diofant',
      version=__version__,
      description='Computer algebra system (CAS) in Python',
      long_description=long_description,
      maintainer='Sergey B Kirpichev',
      license='BSD',
      keywords='Math CAS',
      url='http://diofant.rtfd.io',
      packages=find_packages(),
      ext_modules=[],
      cmdclass={'clean': clean},
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
      ],
      python_requires='>=3.4',
      tests_require=extra_reqs['develop'],
      install_requires=['mpmath>=0.19', 'strategies>=0.2.3', 'cachetools'],
      setup_requires=setup_reqs,
      extras_require=extra_reqs,
)
