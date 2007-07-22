from pyglet.gl import *
from plot_function import PlotFunction, PlotFunctionRegistry, get_vars, count_vars, vrange

class PolarFunction(PlotFunction):

    def __new__(cls, f, intervals, options):
        raise NotImplementedError("Polar functions not implemented.")
        #d = max( [count_vars(f), len(intervals)] )
        #if d == 1:
            #return object.__new__(PolarFunction2d, f, intervals, options)
        #elif d == 2:
            #return object.__new__(PolarFunction3d, f, intervals, options)
        #else:
            #raise ValueError("Cannot plot a polar function with %i variables." % d)

class PolarFunction2d(PlotFunction):
    
    def __init__(self, f, intervals, options):
        self.f = f
        self.intervals = intervals
        self.options = options

    def render(self):
        pass

class PolarFunction3d(PlotFunction):
    
    def __init__(self, f, intervals, options):
        self.f = f
        self.intervals = intervals
        self.options = options

    def render(self):
        pass

PlotFunctionRegistry.register('polar', PolarFunction)
