#!/usr/bin/env python3
"""
Setup script for Diofant.

This script uses Setuptools (https://setuptools.readthedocs.io/en/latest/),
a collection of enhancements to the standard Python distutils.
"""

import re

from setuptools import find_packages, setup


with open('diofant/__init__.py') as f:
    source = f.read()
    long_description = source.split('"""')[1]

    for line in source.splitlines():
        m = re.match('^__version__ = "([0-9ab.]+(dev[0-9]+)?)".*$', line)
        if m:
            __version__ = m.group(1)

setup_reqs = ['setuptools>=5.5.1', 'pip>=6.0', 'pytest-runner', 'isort']
extra_reqs = {'exports': ['numpy>=1.12.1', 'scipy', 'Theano>=0.9.0'],
              'gmpy': ['gmpy2>=2.0.8'],
              'plot': ['pyparsing!=2.1.2', 'matplotlib!=2.1.1'],
              'interactive': ['ipython>=2.3.0'],
              'docs': ['docutils!=0.13.1', 'sphinx>=1.6.7', 'numpydoc',
                       'sphinx_rtd_theme>=0.2.4'],
              }
extra_reqs['develop'] = ['pytest>=3.0', 'flake8>=2.5.5,!=3.1.0',
                         'flake8-docstrings>=1.2.0', 'pydocstyle', 'pep8-naming',
                         'flake8-comprehensions', 'flake8-isort', 'hypothesis',
                         'pytest-xdist', 'pytest-cov', 'pytest-timeout',
                         'coverage'] + setup_reqs

setup(name='Diofant',
      version=__version__,
      description='Computer algebra system (CAS) in Python',
      long_description=long_description,
      maintainer='Sergey B Kirpichev',
      license='BSD',
      keywords='Math CAS',
      url='https://diofant.readthedocs.io/',
      packages=find_packages(),
      ext_modules=[],
      classifiers=[
          'Development Status :: 4 - Beta',
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
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],
      python_requires='>=3.5',
      tests_require=extra_reqs['develop'],
      install_requires=['mpmath>=0.19', 'strategies>=0.2.3', 'cachetools'],
      setup_requires=setup_reqs,
      extras_require=extra_reqs,
      zip_safe=True)
