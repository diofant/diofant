Guidelines for contributing
===========================

.. _reporting-issues:

Reporting issues
----------------

When opening a new issue, please take the following steps:

1. Search `GitHub issues`_ for your issue to avoid duplicate
   reports.  Keyword searches for your error messages are most helpful.

2. If possible, try updating to master and reproducing your issue,
   because we may have already fixed it.

3. Try to include a minimal reproducible test case as an example.

4. Include relevant system information (python version,
   program version, any relevant package versions e.g. mpmath's).

Contributing code
-----------------

All work should be submitted via `Pull Requests (PR)`_.

1. PR can be submitted as soon as there is code worth discussing.

2. Please put your work on the branch of your fork, not
   in the master branch.

3. PR should generally be made against master.

4. One logical change per commit.  Make good commit messages: short
   (<= 78 characters) one-line summary, then newline followed by
   verbose description of your changes.

5. PR should include tests:

   1. Bugfixes should include regression tests.
   2. All new functionality should be tested.
   3. All new public interfaces (methods, functions or classes) should
      have doctests showing how to use them.
   4. Keep in mind, doctests are not tests.  Think of them as
      examples that happen to be tested.

6. Usually, it's good idea to be sure that all tests
   pass and you don't break anything, so please run::

       $ ./setup.py test

License
-------

By submitting a PR, you agree to license your code under the SymPy's
`BSD license`_ and `LGPL`_.


.. _GitHub issues: https://github.com/skirpichev/omg/issues
.. _Pull Requests (PR): https://github.com/skirpichev/omg/pulls
.. _BSD license: LICENSE
.. _LGPL: https://www.gnu.org/copyleft/lesser.html
