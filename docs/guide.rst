.. _guide:

=================
Developer's Guide
=================

.. role:: input(strong)

.. warning::

    If you are new to Diofant, start with the :ref:`Tutorial <tutorial>`.  If
    you are willing to contribute - it's assumed you know the Python
    programming language and the Git Version Control System.

.. include:: CONTRIBUTING.rst

Rosetta Stone
=============

The Diofant project is a `SymPy`_'s fork, so it could be handy to collect here
some facts about SymPy and explain historical code conventions.

First, the SymPy project was hosted in SVN repository on the Google Code and
our master branch include only commits, that added after moving project on the
Github.  But it's not a problem for us - we keep old history on the branch
sympy-svn-history.  Also, you can see this history as part of master's, if you
:ref:`clone our repo <installation-src>` and simply do this::

    $ git fetch origin 'refs/replace/*:refs/replace/*'

Please note, that we have dozens of references to SymPy issues in our
codebase.  Such reference must be either a direct URL of the issue, or
a fully qualified reference in the Github format, like
``sympy/sympy#123``.  Unqualified references like ``#123`` or ``issue
123`` --- are reserved for `Diofant's issues`_.  Functions for
regression tests should be named like ``test_sympyissue_123`` and
``test_diofantissue_123``, respectively.

However, in the old Git history, before commit :commit:`cbdd072`,
please expect that ``#123``, ``issue #123`` or ``issue 123`` --- are
references to the SymPy's issues.  The whole story is a little worse,
because before commit :commit:`6f68fa1` - such unqualified references
assume issues on the Google Code, not Github, unless other clearly
stated.  SymPy issues from the Google Code were moved to the Github in
March 2014 (see :sympyissue:`7235`).  Transfered issue numbers were
shifted by 3099.  I.e. ``issue 123`` in the history - does mean issue
``sympy/sympy#3222`` on Github.

.. _SymPy : http://www.sympy.org/
.. _Diofant's issues : https://github.com/diofant/diofant/issues

Versioning and Release Procedure
================================

We use standard `Semantic Versioning`_ MAJOR.MINOR.MAINTENANCE
numbering scheme, but adopt `PEP 440`_ for alpha ("aN" suffix), beta
("bN"), release candidates ("rcN") and development (".devN").

Releasing a new version is done as follows:

1. Increase ``__version__`` in ``diofant/__init__.py``.

2. Commit version update::

    $ git commit -am 'Bump version to X.Y.Z'

3. Merge commit to the master branch.

4. Tag merge commit and push release tag::

    $ git pull
    $ git tag -s vX.Y.Z
    $ git push origin vX.Y.Z

.. _Semantic Versioning: https://semver.org/
.. _PEP 440: https://www.python.org/dev/peps/pep-0440/

.. XXX adopt http://www.contribution-guide.org

Labeling Issues and Pull Requests
=================================

Following table lists meanings of labels, including
`provided by Github
<https://help.github.com/articles/about-labels/#using-default-labels>`_:

+------------------+----------------------------------------------------+
| Label            | Description                                        |
+==================+====================================================+
| bug              | an unexpected problem or unintended behavior       |
+------------------+----------------------------------------------------+
| wrong answer     | if mathematically wrong result was obtained        |
+------------------+----------------------------------------------------+
| duplicate        | indicates similar issues or pull requests          |
+------------------+----------------------------------------------------+
| enhancement      | new feature requests (or implementation)           |
+------------------+----------------------------------------------------+
| good first issue | indicates a good issue for first-time contributors |
+------------------+----------------------------------------------------+
| help wanted      | indicates that a maintainer wants help here        |
+------------------+----------------------------------------------------+
| invalid          | mark that a problem is no longer relevant          |
+------------------+----------------------------------------------------+
| question         | support request or need for more information       |
+------------------+----------------------------------------------------+
| wontfix          | indicates that work won't continue on this issue   |
+------------------+----------------------------------------------------+

.. include:: README.rst
