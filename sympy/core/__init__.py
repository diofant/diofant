"""Core module. 

Provides the basic operations needed in SymPy. It can be viewed as the smallest
possible self-consistent part (core) of SymPy. It should have a well-defined
interface. Allowing us to possibly rewrite the core in C++ (optionally), or
just use some clever new algorithms to speed it up and the sympy.modules should
work without touching them.
"""
