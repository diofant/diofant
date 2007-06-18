import sys
sys.path.append(".")

try:
    from OpenGL.GL import *
    disabled = False
except:
    disabled = True

from sympy import Symbol, log
x,y = Symbol('x'), Symbol('y')

class TestPlotting:
    
    def __init__(self):
        global disabled
        self.disabled = disabled

    def test_import(self):
        from sympy.modules.plotting import Plot
    
    def test_plot_2d(self):
        from sympy.modules.plotting import Plot
        Plot(x, [x, -5, 5, 10], show=False)
        
    def test_plot_2d_discontinuous(self):
        from sympy.modules.plotting import Plot
        Plot(1/x, [x, -1, 1, 2], show=False)
    
    def test_plot_3d(self):
        from sympy.modules.plotting import Plot
        Plot(x*y, [x, -5, 5, 10], [y, -5, 5, 10], show=False)

    def test_plot_3d_discontinuous(self):
        from sympy.modules.plotting import Plot
        Plot(1/x, [x, -1, 1, 2], [y, -1, 1, 1], show=False)
