#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A tool to help keep .mailmap and aboutus.rst up-to-date.
"""

import os
import re
import sys

from fabric.api import local, env
from fabric.colors import yellow, blue, green, red
from fabric.utils import error

mailmap_update_path = os.path.abspath(__file__)
mailmap_update_dir = os.path.dirname(mailmap_update_path)
sympy_top = os.path.split(mailmap_update_dir)[0]
sympy_dir = os.path.join(sympy_top, 'sympy')

if os.path.isdir(sympy_dir):
    sys.path.insert(0, sympy_top)

from sympy.utilities.misc import filldedent

try:
    # Only works in newer versions of fabric
    env.colorize_errors = True
except AttributeError:
    pass

git_command = 'git log --format="%aN" | sort -u'

git_people = unicode(local(git_command, capture=True), 'utf-8').strip().split("\n")

from distutils.version import LooseVersion

git_ver = local('git --version', capture=True)[12:]
if LooseVersion(git_ver) < LooseVersion('1.8.4.2'):
    print(yellow("Please use a newer git version >= 1.8.4.2"))

with open(os.path.realpath(os.path.join(__file__, os.path.pardir,
                           os.path.pardir, "docs/aboutus.rst"))) as fd:
    AUTHORS = unicode(fd.read(), 'utf-8')

authors = []

for l in AUTHORS.splitlines():
    if l.startswith("#. "):
        authors.append(re.sub(r'^#\. ([^:]+):.*', r'\1', l))

# People who don't want to be listed in aboutus.rst
authors_skip = ["Kirill Smelkov"]

predate_git = 0

exit1 = False

print(blue("""All people who contributed by sending at least a
patch or more (in the order of the date of their first contribution) should
appear in the docs/aboutus.rst (except those who explicitly didn't want to
be mentioned) with an entry like this:

#. Jonh Smith: did this and that.

People with additional space after "#." are not found in the
metadata of the git history."""))

print()
print(yellow("People who are in aboutus.rst but not in git:"))
print()

for name in sorted(set(authors) - set(git_people)):
    if name.startswith(" "):
        # People who are in aboutus.rst but predate git
        predate_git += 1
        continue
    exit1 = True
    print(name)

print()
print(yellow("People who are in git but not in aboutus.rst:"))
print()

for name in sorted(set(git_people) - set(authors) - set(authors_skip)):
    exit1 = True
    print(name)

authors_count = len(authors)
adjusted_authors_count = (
    authors_count
    - predate_git
    + len(authors_skip)
)
git_count = len(git_people)

print()
print(yellow("There are {git_count} people in git, and {adjusted_authors_count} "
    "(adjusted) people from aboutus.rst".format(git_count=git_count,
    adjusted_authors_count=adjusted_authors_count)))

if git_count != adjusted_authors_count:
    error("These two numbers are not the same!")
else:
    print()
    print(green(filldedent("""Congratulations. The aboutus.rst and .mailmap
files appear to be up to date.""")))

if exit1:
    print()
    print(red("There were errors. Please fix them."))
    sys.exit(1)
