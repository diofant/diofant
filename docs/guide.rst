.. _guide:

=================
Developer's Guide
=================

.. role:: input(strong)

.. warning::

    If you are new to Diofant, start with the :ref:`Tutorial <tutorial>`.  If
    you are willing to contribute - it's assumed you know the Python
    programming language and the Git Version Control System.

.. include:: ../CONTRIBUTING.rst

Versioning and Release procedure
================================

We use standard `Semantic Versioning`_ MAJOR.MINOR.MAINTENANCE
numbering scheme, but adopt `PEP 440`_ for alpha ("aN" suffix),
beta ("bN" suffix) and release candidates ("rcN" suffix).

Releasing a new version is done as follows:

1. Increase ``__version__`` in ``diofant/__init__.py``.

2. Commit version update::

    $ git commit -am 'Bump version to X.Y.Z'

3. Merge commit to the master branch.

4. Tag merge commit and push release tag::

    $ git pull
    $ git tag -s vX.Y.Z
    $ git push origin vX.Y.Z

.. _Semantic Versioning: http://semver.org/
.. _PEP 440: https://www.python.org/dev/peps/pep-0440/

.. XXX adopt http://www.contribution-guide.org

.. include:: README.rst
