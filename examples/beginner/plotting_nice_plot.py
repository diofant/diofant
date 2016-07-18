#!/usr/bin/env python

"""Plotting example

Demonstrates simple plotting.
"""

from diofant import cos, sin, log, tan
from diofant.plotting.plot import plot3d

from diofant.abc import x, y


def main():
    fun1 = cos(x)*sin(y)
    fun2 = sin(x)*sin(y)
    fun3 = cos(y) + log(tan(y/2)) + 0.2*x

    plot3d(fun1, fun2, fun3, (x, -0.00, 12.4), (y, 0.1, 2))

if __name__ == "__main__":
    main()
