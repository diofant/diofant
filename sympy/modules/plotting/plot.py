from sympy import Basic, Symbol
from scene import Scene
from sympy.modules.plotting.scene import Renderable
from renderables import CartesianCurve, CartesianSurface, ParametricCurve, Triangle, Simple3dAxes
from util import getkwarg

class Plot(object):
    """
    Plot is the primary user interface for plotting in sympy.

    Usage Examples
    =======
    see examples\plotting.py
    """

    window_title = "SymPy Plot - %s"
    f_str_list = []
    def __init__(self, *fargs, **kwargs):
        """
        
        """
        dist = getkwarg(kwargs, 'dist', 3.0)

        self._scene = Scene(dist = dist)
        self.append(*fargs, **kwargs)
        self._scene.title = self.window_title % (", ".join(self.f_str_list))

        if getkwarg(kwargs, 'axes', True):
            self._scene.append(Simple3dAxes())
        if getkwarg(kwargs, 'show', len(fargs) > 0):
            self.show()

    def append(self, *fargs, **kwargs):

        for functions, intervals in self.parse_args(*fargs):
            for f in functions:
                vars = f.atoms(type=Symbol)
                var_count = len(vars)
                i_count = len(intervals)
                if var_count not in (0, 1, 2):
                    raise ValueError( "Cannot plot %d variables." % (var_count) )
                for v in vars:
                    var_found = False
                    for i in intervals:
                        if i[0] == v:
                            var_found = True
                            break
                    if not var_found:
                        raise ValueError( "No interval given for variable %s." % (str(v)) )
                self.f_str_list.append(str(f))
                if i_count == 2:
                    self._scene.append( CartesianSurface(f, intervals[0], intervals[1]) )
                else:
                    self._scene.append( CartesianCurve(f, intervals[0]) )

    def parse_args(self, *fargs):
        """
        Generator which unpacks arguments to append() (and hence to __init__).
        
        Grammar
        =======
        *fargs := (Renderable|cluster)*
        cluster := sympy-function+, interval+
        interval := [var, var_min, var_max, var_steps]
        
        Ex. In Grammar
        --------------
        sympy-function, interval
        Renderable, sympy-function, interval, interval, Renderable, sympy-function, sympy-function, interval
        
        Ex. Not In Grammar
        ------------------
        interval, sympy-function
        """

        def error_check(functions, intervals):
            """
            Error checking helper which removes redundancy.
            """
            f_error = "No functions specified for interval(s) '%s'"
            i_error = "No interval specified for function(s) '%s'"
            if len(functions) == 0: raise ValueError( f_error % (intervals) )
            if len(intervals) == 0: raise ValueError( i_error % (functions) )
            return (functions, intervals)

        functions = []; intervals = []
        for token in fargs:
            if isinstance(token, Renderable):
                self._scene.append(token)
            else:
                try:
                    if len(token) == 4: #>= 3 and len(token) <= 4:
                        if isinstance(token[0], Symbol): intervals.append(token)
                except:
                    try: f = Basic.sympify(token)
                    except: raise ValueError( "Could not interpret token '%s'" % (token) )

                    if len(intervals) != 0:
                        yield error_check(functions, intervals)

                        functions = []; intervals = []
                    functions.append(f)

        if len(functions) != 0 or len(intervals) != 0:
            yield error_check(functions, intervals)

    def show(self):
        self._scene.show()
