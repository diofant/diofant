Guidelines for Contributing
===========================

This project adheres to `No Code of Conduct`_.  Contributions will
be judged by their technical merit.  Nothing else matters.

.. _reporting-issues:

Reporting Issues
----------------

When opening a new issue, please take the following steps:

1. Please search `GitHub issues`_ to avoid duplicate reports.

2. If possible, try updating to master and reproducing your issue.

3. Try to include a minimal reproducible test case as an example.

4. Include any relevant details of your local setup (i.e. Python
   version, installed libraries).

Contributing Code
-----------------

All work should be submitted via `Pull Requests (PR)`_.

1. PR can be submitted as soon as there is code worth discussing.

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

   1. Bugfixes should include regression tests.
   2. All new functionality should be tested, every new line
      should be covered by tests.
   3. Optionally, provide doctests to illustrate usage.  But keep in
      mind, doctests are not tests.  Think of them as examples that
      happen to be tested.

6. Usually, it's good idea to be sure that all existing tests
   pass and you don't break anything, so please run::

       $ python setup.py test

License
-------

By submitting a PR, you agree to license your code under the
`BSD license`_ and `LGPL`_.

.. _GitHub issues: https://github.com/diofant/diofant/issues
.. _Pull Requests (PR): https://github.com/diofant/diofant/pulls
.. _PEP 8: https://www.python.org/dev/peps/pep-0008/
.. _PEP 257: https://www.python.org/dev/peps/pep-0257/
.. _flake8: http://flake8.rtfd.io/
.. _BSD license: https://github.com/diofant/diofant/blob/master/LICENSE
.. _LGPL: https://www.gnu.org/copyleft/lesser.html
.. _No Code of Conduct: https://github.com/domgetter/NCoC
.. _mention closed issues: https://help.github.com/articles/closing-issues-via-commit-messages
