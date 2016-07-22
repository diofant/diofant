===============
 Preliminaries
===============

This tutorial assumes that the reader already knows the basics of the Python programming
language.  If you do not, the `official Python
tutorial <http://docs.python.org/3/tutorial/index.html>`_ is excellent.

This tutorial assumes a decent mathematical background.  Most examples require
knowledge lower than a calculus level, and some require knowledge at a
calculus level.  Some of the advanced features require more than this. If you
come across a section that uses some mathematical function you are not
familiar with, you can probably skip over it, or replace it with a similar one
that you are more familiar with.  Or look up the function on Wikipedia and
learn something new.  Some important mathematical concepts that are not common
knowledge will be introduced as necessary.

Installation
============

You will need to install Diofant first.  See the :ref:`installation guide
<installation>`.

Exercises
=========

This tutorial was the basis for a tutorial given at the 2013 SciPy conference
in Austin, TX.  The website for that tutorial is `here
<http://certik.github.io/scipy-2013-tutorial/html/index.html>`_. It has links
to videos, materials, and IPython notebook exercises.  The IPython notebook
exercises in particular are recommended to anyone going through this tutorial.

About This Tutorial
===================

This tutorial aims to give an introduction to Diofant for someone who has not
used the library before.  Many features of Diofant will be introduced in this
tutorial, but they will not be exhaustive. In fact, virtually every
functionality shown in this tutorial will have more options or capabilities
than what will be shown.  The rest of the Diofant documentation serves as API
documentation, which extensively lists every feature and option of each
function.

These are the goals of this tutorial:

.. NB: This is mainly here for you, the person who is editing and adding to
   this tutorial. Try to keep these principles in mind.

- To give a guide, suitable for someone who has never used Diofant (but who has
  used Python and knows the necessary mathematics).

- To be written in a narrative format, which is both easy and fun to follow.
  It should read like a book.

- To give insightful examples and exercises, to help the reader learn and to
  make it entertaining to work through.

- To introduce concepts in a logical order.

.. In other words, don't try to get ahead of yourself.

- To use good practices and idioms, and avoid antipatterns.  Functions or
  methodologies that tend to lead to antipatterns are avoided. Features that
  are only useful to advanced users are not shown.

- To be consistent.  If there are multiple ways to do it, only the best way is
  shown.

.. For example, there are at least five different ways to create Symbols.
   ``symbols`` is the only one that is general and doesn't lead to
   antipatterns, so it is the only one used.

- To avoid unnecessary duplication, it is assumed that previous sections of
  the tutorial have already been read.

Feedback on this tutorial, or on Diofant in general is always welcome.
