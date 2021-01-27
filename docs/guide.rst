.. _guide:

===========
Development
===========

.. role:: input(strong)

.. warning::

    If you are new to Diofant, start with the :ref:`Tutorial <tutorial>`.  If
    you are willing to contribute - it's assumed you know the Python
    programming language and the Git Version Control System.

This project adheres to `No Code of Conduct`_.  Contributions will
be judged by their technical merit.  Nothing else matters.

.. _reporting-issues:

Reporting Issues
================

When opening a new issue, please take the following steps:

1. Please search `GitHub issues`_ to avoid duplicate reports.

2. If possible, try updating to master and reproducing your issue.

3. Try to include a minimal reproducible code example that demonstrates the
   problem.  Avoid screenshots, if possible.

4. Include any relevant details of your local setup (i.e. Python
   version, installed libraries).

.. note::

    Please avoid changing your messages on the GitHub, unless you
    want fix a typo and so on.  Just add new comments.

Contributing Code
=================

All work should be submitted via `Pull Requests (PR)`_.

1. PR can be submitted as soon as there is code worth discussing.
   Please make a draft PR, if one is not intended to be merged
   in its present shape even if all checks pass.

2. Please put your work on the branch of your fork, not in the
   master branch.  PR should generally be made against master.

3. One logical change per commit.  Make good commit messages: short
   (<= 78 characters) one-line summary, then newline followed by
   verbose description of your changes.  Please `mention closed
   issues`_ with commit message.

4. Please conform to `PEP 8`_ and `PEP 257`_, enable `flake8 git hook
   <http://flake8.pycqa.org/en/stable/user/using-hooks.html>`_ to
   prevent badly formatted commits.

5. PR should include tests:

   1. Bugfixes should include regression tests.  Please format
      them accordingly, to include references for fixed
      issues (e.g. by naming test like ``test_diofantissue_123``
      or adding comment with issue number).
   2. All new functionality should be tested, every new line
      should be covered by tests.  Please use in tests only
      public interfaces.  Regression tests are not accounted in
      the coverage statistics.
   3. Optionally, provide doctests to illustrate usage.  But keep in
      mind, doctests are not tests.  Think of them as examples that
      happen to be tested.

6. It's good idea to be sure that **all** existing tests
   pass and you don't break anything, so please run::

       pytest

   To check also doctests, run::

       pytest --doctest-modules

7. Please also check for potential flaws in your Python code with::

       pylint diofant

   and do type checking::

       mypy diofant

8. If your change affects documentation, please build it by::

       python setup.py build_sphinx -W

   and check that it looks as expected.

Rosetta Stone
=============

The Diofant project is a `SymPy`_'s fork, so it could be handy to collect here
some facts about SymPy and explain historical code conventions.

First, the SymPy project was hosted in SVN repository on the Google Code and
our master branch include only commits, that added after moving project on the
Github.  But it's not a problem for us - we keep old history on the branch
sympy-svn-history.  Also, you can see this history as part of master's, if you
:ref:`clone our repo <installation-src>` and simply do this::

    git fetch origin 'refs/replace/*:refs/replace/*'

Please note, that we have dozens of references to SymPy issues in our
codebase.  Such reference must be either a direct URL of the issue, or
a fully qualified reference in the Github format, like
``sympy/sympy#123``.  Unqualified references like ``#123`` or ``issue
123`` --- are reserved for our `Github issues`_.  Functions for
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

Versioning and Release Procedure
================================

We use standard `Semantic Versioning`_ numbering scheme, but adopt
`PEP 440`_ for alpha ("aN" suffix), beta ("bN") and development
(".devN") releases.

To release a new version, tag latest commit to the master branch
and publish this release tag::

    git pull
    git tag -s vX.Y.Z
    git push origin vX.Y.Z

.. XXX adopt http://www.contribution-guide.org

.. _SymPy : https://www.sympy.org/
.. _Semantic Versioning: https://semver.org/
.. _PEP 440: https://www.python.org/dev/peps/pep-0440/
.. _GitHub issues: https://github.com/diofant/diofant/issues
.. _Pull Requests (PR): https://github.com/diofant/diofant/pulls
.. _PEP 8: https://www.python.org/dev/peps/pep-0008/
.. _PEP 257: https://www.python.org/dev/peps/pep-0257/
.. _flake8: http://flake8.rtfd.io/
.. _No Code of Conduct: https://github.com/domgetter/NCoC
.. _mention closed issues: https://help.github.com/en/github/managing-your-work-on-github/linking-a-pull-request-to-an-issue
