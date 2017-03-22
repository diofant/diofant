============
SymPy 0.5.15
============

24 May 2008

* all SymPy functions support vector arguments, e.g. ``sin([1, 2, 3])``
* lambdify can now use numpy/math/mpmath
* the order of lambdify arguments has changed
* all SymPy objects are pickable
* simplify improved and made more robust
* broken limit_series was removed, we now have just one limit implementation
* limits now use ``.nseries()``
* ``.nseries()`` improved a lot
* Polys improved
* Basic kronecker delta and Levi-Civita implementation
