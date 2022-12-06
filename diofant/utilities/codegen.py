"""
module for generating C, C++, Fortran77, Fortran90 and Octave/Matlab routines
that evaluate diofant expressions.  This module is work in progress.  Only the
milestones with a '+' character in the list below have been completed.

--- How is diofant.utilities.codegen different from diofant.printing.ccode? ---

We considered the idea to extend the printing routines for diofant functions in
such a way that it prints complete compilable code, but this leads to a few
unsurmountable issues that can only be tackled with dedicated code generator:

- For C, one needs both a code and a header file, while the printing routines
  generate just one string. This code generator can be extended to support
  .pyf files for f2py.

- Diofant functions are not concerned with programming-technical issues, such
  as input, output and input-output arguments. Other examples are contiguous
  or non-contiguous arrays, including headers of other libraries such as gsl
  or others.

- It is highly interesting to evaluate several diofant functions in one C
  routine, eventually sharing common intermediate results with the help
  of the cse routine. This is more than just printing.

- From the programming perspective, expressions with constants should be
  evaluated in the code generator as much as possible. This is different
  for printing.

--- Basic assumptions ---

* A generic Routine data structure describes the routine that must be
  translated into C/Fortran/... code. This data structure covers all
  features present in one or more of the supported languages.

* Descendants from the CodeGen class transform multiple Routine instances
  into compilable code. Each derived class translates into a specific
  language.

* In many cases, one wants a simple workflow. The friendly functions in the
  last part are a simple api on top of the Routine/CodeGen stuff. They are
  easier to use, but are less powerful.

--- Milestones ---

+ First working version with scalar input arguments, generating C code,
  tests
+ Friendly functions that are easier to use than the rigorous
  Routine/CodeGen workflow.
+ Integer and Real numbers as input and output
+ Output arguments
+ InputOutput arguments
+ Sort input/output arguments properly
+ Contiguous array arguments (numpy matrices)
+ Also generate .pyf code for f2py (in autowrap module)
+ Isolate constants and evaluate them beforehand in double precision
+ Fortran 90
+ Octave/Matlab

- Common Subexpression Elimination
- User defined comments in the generated code
- Optional extra include lines for libraries/objects that can eval special
  functions
- Test other C compilers and libraries: gcc, tcc, libtcc, gcc+gsl, ...
- Contiguous array arguments (diofant matrices)
- Non-contiguous array arguments (diofant matrices)
- ccode must raise an error when it encounters something that can not be
  translated into c. ccode(integrate(sin(x)/x, x)) does not make sense.
- Complex numbers as input and output
- A default complex datatype
- Include extra information in the header: date, user, hostname, sha1
  hash, ...
- Fortran 77
- C++
- Python
- ...

"""

import os
import textwrap
from io import StringIO

from .. import __version__ as diofant_version
from ..core import Dummy, Equality, Expr, Function, Integer, Symbol, Tuple
from ..core.compatibility import is_sequence
from ..matrices import (ImmutableMatrix, MatrixBase, MatrixExpr, MatrixSlice,
                        MatrixSymbol)
from ..printing.ccode import CCodePrinter, ccode
from ..printing.fcode import FCodePrinter, fcode
from ..printing.octave import OctaveCodePrinter, octave_code
from ..tensor import Idx, Indexed, IndexedBase


__all__ = (
    # description of routines
    'Routine', 'DataType', 'default_datatypes', 'get_default_datatype',
    'Argument', 'InputArgument', 'OutputArgument', 'Result',
    # routines -> code
    'CodeGen', 'CCodeGen', 'FCodeGen', 'OctaveCodeGen',
    # friendly functions
    'codegen', 'make_routine',
)


#
# Description of routines
#


class Routine:
    """Generic description of evaluation routine for set of expressions.

    A CodeGen class can translate instances of this class into code in a
    particular language.  The routine specification covers all the features
    present in these languages.  The CodeGen part must raise an exception
    when certain features are not present in the target language.  For
    example, multiple return values are possible in Python, but not in C or
    Fortran.  Another example: Fortran and Python support complex numbers,
    while C does not.

    """

    def __init__(self, name, arguments, results, local_vars, global_vars):
        """Initialize a Routine instance.

        Parameters
        ==========

        name : string
            Name of the routine.

        arguments : list of Arguments
            These are things that appear in arguments of a routine, often
            appearing on the right-hand side of a function call.  These are
            commonly InputArguments but in some languages, they can also be
            OutputArguments or InOutArguments (e.g., pass-by-reference in C
            code).

        results : list of Results
            These are the return values of the routine, often appearing on
            the left-hand side of a function call.  The difference between
            Results and OutputArguments and when you should use each is
            language-specific.

        local_vars : list of Results
            These are variables that will be defined at the beginning of the
            function.

        global_vars : list of Symbols
            Variables which will not be passed into the function.

        """
        # extract all input symbols and all symbols appearing in an expression
        input_symbols = set()
        symbols = set()
        for arg in arguments:
            if isinstance(arg, OutputArgument):
                symbols.update(arg.expr.free_symbols)
            elif isinstance(arg, InputArgument):
                input_symbols.add(arg.name)
            elif isinstance(arg, InOutArgument):
                input_symbols.add(arg.name)
                symbols.update(arg.expr.free_symbols)
            else:
                raise ValueError(f'Unknown Routine argument: {arg}')

        for r in results:
            if not isinstance(r, Result):
                raise ValueError(f'Unknown Routine result: {r}')
            symbols.update(r.expr.free_symbols)

        local_symbols = set()
        for r in local_vars:
            if isinstance(r, Result):
                symbols.update(r.expr.free_symbols)
                local_symbols.add(r.name)
            else:
                local_symbols.add(r)

        # Check that all symbols in the expressions are covered by
        # InputArguments/InOutArguments---subset because user could
        # specify additional (unused) InputArguments or local_vars.
        notcovered = symbols.difference(input_symbols | local_symbols | global_vars)
        if notcovered != set():
            raise ValueError('Symbols needed for output are not in input ' +
                             ', '.join([str(x) for x in notcovered]))

        self.name = name
        self.arguments = arguments
        self.results = results
        self.local_vars = local_vars
        self.global_vars = global_vars

    @property
    def variables(self):
        """Returns a set of all variables possibly used in the routine.

        For routines with unnamed return values, the dummies that may or
        may not be used will be included in the set.

        """
        v = set(self.local_vars)
        for arg in self.arguments:
            v.add(arg.name)
        for res in self.results:
            v.add(res.result_var)
        return v

    @property
    def result_variables(self):
        """Returns a list of OutputArgument, InOutArgument and Result.

        If return values are present, they are at the end ot the list.

        """
        args = [arg for arg in self.arguments if isinstance(
            arg, (OutputArgument, InOutArgument))]
        args.extend(self.results)
        return args


class DataType:
    """Holds strings for a certain datatype in different languages."""

    def __init__(self, cname, fname, pyname, octname):
        self.cname = cname
        self.fname = fname
        self.pyname = pyname
        self.octname = octname


default_datatypes = {
    'int': DataType('int', 'INTEGER*4', 'int', ''),
    'float': DataType('double', 'REAL*8', 'float', '')
}


def get_default_datatype(expr):
    """Derives an appropriate datatype based on the expression."""
    if expr.is_integer:
        return default_datatypes['int']
    else:
        return default_datatypes['float']


class Variable:
    """Represents a typed variable."""

    def __init__(self, name, datatype=None, dimensions=None, precision=None):
        """Return a new variable.

        Parameters
        ==========

        name : Symbol or MatrixSymbol

        datatype : optional
            When not given, the data type will be guessed based on the
            assumptions on the symbol argument.

        dimension : sequence containing tupes, optional
            If present, the argument is interpreted as an array, where this
            sequence of tuples specifies (lower, upper) bounds for each
            index of the array.

        precision : int, optional
            Controls the precision of floating point constants.

        """
        if not isinstance(name, (Dummy, Symbol, MatrixSymbol)):
            raise TypeError('The first argument must be a diofant symbol.')
        if datatype is None:
            datatype = get_default_datatype(name)
        elif not isinstance(datatype, DataType):
            raise TypeError("The (optional) `datatype' argument must be an "
                            'instance of the DataType class.')
        if dimensions and not isinstance(dimensions, (tuple, list)):
            raise TypeError(
                'The dimension argument must be a sequence of tuples')

        self._name = name
        self._datatype = {
            'C': datatype.cname,
            'FORTRAN': datatype.fname,
            'OCTAVE': datatype.octname,
            'PYTHON': datatype.pyname
        }
        self.dimensions = dimensions
        self.precision = precision

    @property
    def name(self):
        return self._name

    def get_datatype(self, language):
        """Returns the datatype string for the requested language.

        Examples
        ========

        >>> x = Variable(Symbol('x'))
        >>> x.get_datatype('c')
        'double'
        >>> x.get_datatype('fortran')
        'REAL*8'

        """
        try:
            return self._datatype[language.upper()]
        except KeyError as exc:
            raise CodeGenError(f"Has datatypes for languages: {', '.join(self._datatype)}") from exc


class Argument(Variable):
    """An abstract Argument data structure: a name and a data type.

    This structure is refined in the descendants below.

    """


class InputArgument(Argument):
    """Input argument class."""


class ResultBase:
    """Base class for all "outgoing" information from a routine.

    Objects of this class stores a diofant expression, and a diofant object
    representing a result variable that will be used in the generated code
    only if necessary.

    """

    def __init__(self, expr, result_var):
        self.expr = expr
        self.result_var = result_var


class OutputArgument(Argument, ResultBase):
    """OutputArgument are always initialized in the routine."""

    def __init__(self, name, result_var, expr, datatype=None, dimensions=None, precision=None):
        """Return a new variable.

        Parameters
        ==========

        name : Symbol, MatrixSymbol
            The name of this variable.  When used for code generation, this
            might appear, for example, in the prototype of function in the
            argument list.

        result_var : Symbol, Indexed
            Something that can be used to assign a value to this variable.
            Typically the same as `name` but for Indexed this should be e.g.,
            "y[i]" whereas `name` should be the Symbol "y".

        expr : object
            The expression that should be output, typically a Diofant
            expression.

        datatype : optional
            When not given, the data type will be guessed based on the
            assumptions on the symbol argument.

        dimension : sequence containing tupes, optional
            If present, the argument is interpreted as an array, where this
            sequence of tuples specifies (lower, upper) bounds for each
            index of the array.

        precision : int, optional
            Controls the precision of floating point constants.

        """
        Argument.__init__(self, name, datatype, dimensions, precision)
        ResultBase.__init__(self, expr, result_var)


class InOutArgument(Argument, ResultBase):
    """InOutArgument are never initialized in the routine."""

    def __init__(self, name, result_var, expr, datatype=None, dimensions=None, precision=None):
        if not datatype:
            datatype = get_default_datatype(expr)
        Argument.__init__(self, name, datatype, dimensions, precision)
        ResultBase.__init__(self, expr, result_var)
    __init__.__doc__ = OutputArgument.__init__.__doc__


class Result(Variable, ResultBase):
    """An expression for a return value.

    The name result is used to avoid conflicts with the reserved word
    "return" in the python language.  It is also shorter than ReturnValue.

    These may or may not need a name in the destination (e.g., "return(x*y)"
    might return a value without ever naming it).

    """

    def __init__(self, expr, name=None, result_var=None, datatype=None,
                 dimensions=None, precision=None):
        """Initialize a return value.

        Parameters
        ==========

        expr : Diofant expression

        name : Symbol, MatrixSymbol, optional
            The name of this return variable.  When used for code generation,
            this might appear, for example, in the prototype of function in a
            list of return values.  A dummy name is generated if omitted.

        result_var : Symbol, Indexed, optional
            Something that can be used to assign a value to this variable.
            Typically the same as `name` but for Indexed this should be e.g.,
            "y[i]" whereas `name` should be the Symbol "y".  Defaults to
            `name` if omitted.

        datatype : optional
            When not given, the data type will be guessed based on the
            assumptions on the symbol argument.

        dimension : sequence containing tupes, optional
            If present, this variable is interpreted as an array,
            where this sequence of tuples specifies (lower, upper)
            bounds for each index of the array.

        precision : int, optional
            Controls the precision of floating point constants.

        """
        if not isinstance(expr, (Expr, MatrixBase)):
            raise TypeError('The first argument must be a diofant expression.')

        if name is None:
            name = f'result_{abs(hash(expr))}'

        if isinstance(name, str):
            if isinstance(expr, (MatrixBase, MatrixExpr)):
                name = MatrixSymbol(name, *expr.shape)
            else:
                name = Symbol(name)

        if result_var is None:
            result_var = name

        Variable.__init__(self, name, datatype=datatype,
                          dimensions=dimensions, precision=precision)
        ResultBase.__init__(self, expr, result_var)


#
# Transformation of routine objects into code
#

class CodeGen:
    """Abstract class for the code generators."""

    printer = None  # will be set to an instance of a CodePrinter subclass

    def __init__(self, project='project', cse=False):
        """Initialize a code generator.

        Derived classes will offer more options that affect the generated
        code.

        """
        self.project = project
        self.cse = cse

    def routine(self, name, expr, argument_sequence, global_vars=None):
        """Creates an Routine object that is appropriate for this language.

        This implementation is appropriate for at least C/Fortran.  Subclasses
        can override this if necessary.

        Here, we assume at most one return value (the l-value) which must be
        scalar.  Additional outputs are OutputArguments (e.g., pointers on
        right-hand-side or pass-by-reference).  Matrices are always returned
        via OutputArguments.  If ``argument_sequence`` is None, arguments will
        be ordered alphabetically, but with all InputArguments first, and then
        OutputArgument and InOutArguments.

        """
        if self.cse:
            from ..simplify import cse

            if is_sequence(expr) and not isinstance(expr, (MatrixBase, MatrixExpr)):
                if not expr:
                    raise ValueError('No expression given')
                for e in expr:
                    if not e.is_Equality:
                        raise CodeGenError(f'Lists of expressions must all be Equalities. {e} is not.')

                # create a list of right hand sides and simplify them
                rhs = [e.rhs for e in expr]
                common, simplified = cse(rhs)

                # pack the simplified expressions back up with their left hand sides
                expr = [Equality(e.lhs, rhs) for e, rhs in zip(expr, simplified)]
            else:
                if isinstance(expr, Equality):
                    common, simplified = cse(expr.rhs)
                    expr = Equality(expr.lhs, simplified[0])
                else:
                    common, simplified = cse(expr)
                    expr = simplified

            local_vars = [Result(b, a) for a, b in common]
            local_symbols = {a for a, _ in common}
            local_expressions = Tuple(*[b for _, b in common])
        else:
            local_expressions = Tuple()

        if is_sequence(expr) and not isinstance(expr, (MatrixBase, MatrixExpr)):
            if not expr:
                raise ValueError('No expression given')
            expressions = Tuple(*expr)
        else:
            expressions = Tuple(expr)

        # local variables
        if not self.cse:
            # local variables for indexed expressions
            local_vars = {i.label for i in expressions.atoms(Idx)}
            local_symbols = local_vars

        # global variables
        global_vars = set() if global_vars is None else set(global_vars)

        # symbols that should be arguments
        symbols = (expressions.free_symbols | local_expressions.free_symbols) - local_symbols - global_vars
        # Decide whether to use output argument or return value
        return_val = []
        output_args = []
        for expr in expressions:
            if isinstance(expr, Equality):
                out_arg = expr.lhs
                expr = expr.rhs
                if isinstance(out_arg, Indexed):
                    dims = tuple((Integer(0), dim - 1) for dim in out_arg.shape)
                    symbol = out_arg.base.label
                elif isinstance(out_arg, Symbol):
                    dims = []
                    symbol = out_arg
                elif isinstance(out_arg, MatrixSymbol):
                    dims = tuple((Integer(0), dim - 1) for dim in out_arg.shape)
                    symbol = out_arg
                else:
                    raise CodeGenError('Only Indexed, Symbol, or MatrixSymbol '
                                       'can define output arguments.')

                if expr.has(symbol):
                    output_args.append(
                        InOutArgument(symbol, out_arg, expr, dimensions=dims))
                else:
                    output_args.append(
                        OutputArgument(symbol, out_arg, expr, dimensions=dims))

                # avoid duplicate arguments
                symbols.remove(symbol)
            elif isinstance(expr, (ImmutableMatrix, MatrixSlice)):
                # Create a "dummy" MatrixSymbol to use as the Output arg
                out_arg = MatrixSymbol(f'out_{abs(hash(expr))}', *expr.shape)
                dims = tuple((Integer(0), dim - 1) for dim in out_arg.shape)
                output_args.append(
                    OutputArgument(out_arg, out_arg, expr, dimensions=dims))
            else:
                return_val.append(Result(expr))

        arg_list = []

        # setup input argument list
        array_symbols = {}
        for array in expressions.atoms(Indexed) | local_expressions.atoms(Indexed):
            array_symbols[array.base.label] = array
        for array in expressions.atoms(MatrixSymbol) | local_expressions.atoms(MatrixSymbol):
            array_symbols[array] = array

        for symbol in sorted(symbols, key=str):
            if symbol in array_symbols:
                dims = []
                array = array_symbols[symbol]
                for dim in array.shape:
                    dims.append((Integer(0), dim - 1))
                metadata = {'dimensions': dims}
            else:
                metadata = {}

            arg_list.append(InputArgument(symbol, **metadata))

        output_args.sort(key=lambda x: str(x.name))
        arg_list.extend(output_args)

        if argument_sequence is not None:
            # if the user has supplied IndexedBase instances, we'll accept that
            new_sequence = []
            for arg in argument_sequence:
                if isinstance(arg, IndexedBase):
                    new_sequence.append(arg.label)
                else:
                    new_sequence.append(arg)
            argument_sequence = new_sequence

            missing = [x for x in arg_list if x.name not in argument_sequence]
            if missing:
                msg = "Argument list didn't specify: {0} "
                msg = msg.format(', '.join([str(m.name) for m in missing]))
                raise CodeGenArgumentListError(msg, missing)

            # create redundant arguments to produce the requested sequence
            name_arg_dict = {x.name: x for x in arg_list}
            new_args = []
            for symbol in argument_sequence:
                try:
                    new_args.append(name_arg_dict[symbol])
                except KeyError:
                    dims = None
                    if isinstance(symbol, MatrixSymbol):
                        dims = tuple((Integer(0), dim - 1) for dim in symbol.shape)
                    new_args.append(InputArgument(symbol, dimensions=dims))
            arg_list = new_args

        return Routine(name, arg_list, return_val, local_vars, global_vars)

    def write(self, routines, prefix, to_files=False, header=True, empty=True):
        """Writes all the source code files for the given routines.

        The generated source is returned as a list of (filename, contents)
        tuples, or is written to files (see below).  Each filename consists
        of the given prefix, appended with an appropriate extension.

        Parameters
        ==========

        routines : list
            A list of Routine instances to be written

        prefix : string
            The prefix for the output files

        to_files : bool, optional
            When True, the output is written to files.  Otherwise, a list
            of (filename, contents) tuples is returned.  [default: False]

        header : bool, optional
            When True, a header comment is included on top of each source
            file. [default: True]

        empty : bool, optional
            When True, empty lines are included to structure the source
            files. [default: True]

        """
        for routine in routines:
            if not isinstance(routine, Routine):
                raise CodeGenError(f'Routine expected, got {routine}')

        if to_files:
            for dump_fn in self.dump_fns:
                filename = f'{prefix}.{dump_fn.extension}'
                with open(filename, 'w', encoding='utf-8') as f:
                    dump_fn(self, routines, f, prefix, header, empty)
        else:
            result = []
            for dump_fn in self.dump_fns:
                filename = f'{prefix}.{dump_fn.extension}'
                contents = StringIO()
                dump_fn(self, routines, contents, prefix, header, empty)
                result.append((filename, contents.getvalue()))
            return result

    def dump_code(self, routines, f, prefix, header=True, empty=True):
        """Write the code by calling language specific methods.

        The generated file contains all the definitions of the routines in
        low-level code and refers to the header file if appropriate.

        Parameters
        ==========

        routines : list
            A list of Routine instances.

        f : file-like
            Where to write the file.

        prefix : string
            The filename prefix, used to refer to the proper header file.
            Only the basename of the prefix is used.

        header : bool, optional
            When True, a header comment is included on top of each source
            file.  [default : True]

        empty : bool, optional
            When True, empty lines are included to structure the source
            files.  [default : True]

        """
        code_lines = self._preprocessor_statements(prefix)

        for routine in routines:
            if empty:
                code_lines.append('\n')
            code_lines.extend(self._get_routine_opening(routine))
            code_lines.extend(self._declare_arguments(routine))
            code_lines.extend(self._declare_globals(routine))
            code_lines.extend(self._declare_locals(routine))
            if empty:
                code_lines.append('\n')
            code_lines.extend(self._call_printer(routine))
            if empty:
                code_lines.append('\n')
            code_lines.extend(self._get_routine_ending(routine))

        code_lines = self._indent_code(''.join(code_lines))

        if header:
            code_lines = ''.join(self._get_header() + [code_lines])

        if code_lines:
            f.write(code_lines)

    def _printer_method_with_settings(self, method, settings=None, *args, **kwargs):
        settings = settings or {}
        ori = {k: self.printer._settings[k] for k in settings}
        for k, v in settings.items():
            self.printer._settings[k] = v
        result = getattr(self.printer, method)(*args, **kwargs)
        for k, v in ori.items():
            self.printer._settings[k] = v
        return result


class CodeGenError(Exception):
    pass


class CodeGenArgumentListError(Exception):
    @property
    def missing_args(self):
        return self.args[1]


header_comment = """Code generated with diofant %(version)s

See https://diofant.readthedocs.io/ for more information.

This file is part of '%(project)s'
"""


class CCodeGen(CodeGen):
    """Generator for C code.

    The .write() method inherited from CodeGen will output a code file and
    an interface file, <prefix>.c and <prefix>.h respectively.

    """

    code_extension = 'c'
    interface_extension = 'h'

    def __init__(self, project='project', printer=None,
                 preprocessor_statements=None, cse=False):
        super().__init__(project=project, cse=cse)
        self.printer = printer or CCodePrinter()

        self.preprocessor_statements = preprocessor_statements
        if preprocessor_statements is None:
            self.preprocessor_statements = ['#include <math.h>\n']

    def _get_header(self):
        """Writes a common header for the generated files."""
        code_lines = []
        code_lines.append('/' + '*'*78 + '\n')
        tmp = header_comment % {'version': diofant_version,
                                'project': self.project}
        for line in tmp.splitlines():
            code_lines.append(f' *{line.center(76)}*\n')
        code_lines.append(' ' + '*'*78 + '/\n')
        return code_lines

    def get_prototype(self, routine):
        """Returns a string for the function prototype of the routine.

        If the routine has multiple result objects, an CodeGenError is
        raised.

        See: https://en.wikipedia.org/wiki/Function_prototype

        """
        if len(routine.results) > 1:
            raise CodeGenError('C only supports a single or no return value.')
        if len(routine.results) == 1:
            ctype = routine.results[0].get_datatype('C')
        else:
            ctype = 'void'

        type_args = []
        for arg in routine.arguments:
            name = ccode(arg.name)
            if arg.dimensions or isinstance(arg, ResultBase):
                type_args.append((arg.get_datatype('C'), f'*{name}'))
            else:
                type_args.append((arg.get_datatype('C'), name))
        arguments = ', '.join([f"{t[0] + ' ' + t[1]}" for t in type_args])
        return f'{ctype} {routine.name}({arguments})'

    def _preprocessor_statements(self, prefix):
        code_lines = []
        code_lines.append(f'#include "{os.path.basename(prefix)}.h\"\n')
        code_lines.extend(self.preprocessor_statements)
        return code_lines

    def _get_routine_opening(self, routine):
        prototype = self.get_prototype(routine)
        return [f'{prototype} {{\n']

    def _declare_arguments(self, routine):
        # arguments are declared in prototype
        return []

    def _declare_globals(self, routine):
        # global variables are not explicitly declared within C functions
        return []

    def _declare_locals(self, routine):
        # Compose a list of symbols to be dereferenced in the function
        # body. These are the arguments that were passed by a reference
        # pointer, excluding arrays.
        dereference = []
        for arg in routine.arguments:
            if isinstance(arg, ResultBase) and not arg.dimensions:
                dereference.append(arg.name)

        code_lines = []
        for result in routine.local_vars:

            # local variables that are simple symbols such as those used as indices into
            # for loops are defined declared elsewhere.
            if not isinstance(result, Result):
                continue

            assign_to = result.name
            t = result.get_datatype('c')
            prefix = f'const {t} '

            *_, c_expr = self._printer_method_with_settings(
                'doprint', {'human': False, 'dereference': dereference},
                result.expr, assign_to=assign_to)

            code_lines.append(f'{prefix}{c_expr}\n')

        return code_lines

    def _call_printer(self, routine):
        code_lines = []

        # Compose a list of symbols to be dereferenced in the function
        # body. These are the arguments that were passed by a reference
        # pointer, excluding arrays.
        dereference = []
        for arg in routine.arguments:
            if isinstance(arg, ResultBase) and not arg.dimensions:
                dereference.append(arg.name)

        return_val = None
        for result in routine.result_variables:
            if isinstance(result, Result):
                assign_to = routine.name + '_result'
                t = result.get_datatype('c')
                code_lines.append(f'{t} {assign_to!s};\n')
                return_val = assign_to
            else:
                assign_to = result.result_var

            constants, _, c_expr = ccode(result.expr, human=False,
                                         assign_to=assign_to, dereference=dereference)

            for name, value in sorted(constants, key=str):
                code_lines.append(f'double const {name} = {value};\n')
            code_lines.append(f'{c_expr}\n')

        if return_val:
            code_lines.append(f'   return {return_val};\n')
        return code_lines

    def _indent_code(self, codelines):
        p = CCodePrinter()
        return p.indent_code(codelines)

    def _get_routine_ending(self, routine):
        return ['}\n']

    def dump_c(self, routines, f, prefix, header=True, empty=True):
        self.dump_code(routines, f, prefix, header, empty)
    dump_c.extension = code_extension  # type: ignore[attr-defined]
    dump_c.__doc__ = CodeGen.dump_code.__doc__

    def dump_h(self, routines, f, prefix, header=True, empty=True):
        """Writes the C header file.

        This file contains all the function declarations.

        Parameters
        ==========

        routines : list
            A list of Routine instances.

        f : file-like
            Where to write the file.

        prefix : string
            The filename prefix, used to construct the include guards.
            Only the basename of the prefix is used.

        header : bool, optional
            When True, a header comment is included on top of each source
            file.  [default : True]

        empty : bool, optional
            When True, empty lines are included to structure the source
            files.  [default : True]

        """
        if header:
            print(''.join(self._get_header()), file=f)
        guard_name = f"{self.project.replace(' ', '_').upper()}__{prefix.replace('/', '_').upper()}__H"
        # include guards
        if empty:
            print(file=f)
        print(f'#ifndef {guard_name}', file=f)
        print(f'#define {guard_name}', file=f)
        if empty:
            print(file=f)
        # declaration of the function prototypes
        for routine in routines:
            prototype = self.get_prototype(routine)
            print(f'{prototype};', file=f)
        # end if include guards
        if empty:
            print(file=f)
        print('#endif', file=f)
        if empty:
            print(file=f)
    dump_h.extension = interface_extension  # type: ignore[attr-defined]

    # This list of dump functions is used by CodeGen.write to know which dump
    # functions it has to call.
    dump_fns = [dump_c, dump_h]


class FCodeGen(CodeGen):
    """Generator for Fortran 95 code

    The .write() method inherited from CodeGen will output a code file and
    an interface file, <prefix>.f90 and <prefix>.h respectively.

    """

    code_extension = 'f90'
    interface_extension = 'h'

    def __init__(self, project='project'):
        CodeGen.__init__(self, project)

    def _get_symbol(self, s):
        """Returns the symbol as fcode prints it."""
        return fcode(s).strip()

    def _get_header(self):
        """Writes a common header for the generated files."""
        code_lines = []
        code_lines.append('!' + '*'*78 + '\n')
        tmp = header_comment % {'version': diofant_version,
                                'project': self.project}
        for line in tmp.splitlines():
            code_lines.append(f'!*{line.center(76)}*\n')
        code_lines.append('!' + '*'*78 + '\n')
        return code_lines

    def _preprocessor_statements(self, prefix):
        return []

    def _get_routine_opening(self, routine):
        """Returns the opening statements of the fortran routine."""
        code_list = []
        if len(routine.results) > 1:
            raise CodeGenError('Fortran only supports a single '
                               'or no return value.')
        if len(routine.results) == 1:
            result = routine.results[0]
            code_list.append(result.get_datatype('fortran'))
            code_list.append('function')
        else:
            code_list.append('subroutine')

        args = ', '.join(f'{self._get_symbol(arg.name)}'
                         for arg in routine.arguments)

        call_sig = f'{routine.name}({args})\n'
        # Fortran 95 requires all lines be less than 132 characters, so wrap
        # this line before appending.
        call_sig = ' &\n'.join(textwrap.wrap(call_sig,
                                             width=60,
                                             break_long_words=False)) + '\n'
        code_list.append(call_sig)
        code_list = [' '.join(code_list)]
        code_list.append('implicit none\n')
        return code_list

    def _declare_arguments(self, routine):
        # argument type declarations
        code_list = []
        array_list = []
        scalar_list = []
        for arg in routine.arguments:

            if isinstance(arg, InputArgument):
                typeinfo = f"{arg.get_datatype('fortran')}, intent(in)"
            elif isinstance(arg, InOutArgument):
                typeinfo = f"{arg.get_datatype('fortran')}, intent(inout)"
            else:
                typeinfo = f"{arg.get_datatype('fortran')}, intent(out)"

            fprint = self._get_symbol

            if arg.dimensions:
                # fortran arrays start at 1
                dimstr = ', '.join([f'{fprint(dim[0] + 1)}:{fprint(dim[1] + 1)}'
                                    for dim in arg.dimensions])
                typeinfo += f', dimension({dimstr})'
                array_list.append(f'{typeinfo} :: {fprint(arg.name)}\n')
            else:
                scalar_list.append(f'{typeinfo} :: {fprint(arg.name)}\n')

        # scalars first, because they can be used in array declarations
        code_list.extend(scalar_list)
        code_list.extend(array_list)

        return code_list

    def _declare_globals(self, routine):
        # Global variables not explicitly declared within Fortran 90 functions.
        # Note: a future F77 mode may need to generate "common" blocks.
        return []

    def _declare_locals(self, routine):
        code_list = []
        for var in sorted(routine.local_vars, key=str):
            typeinfo = get_default_datatype(var)
            code_list.append(f'{typeinfo.fname} :: {self._get_symbol(var)}\n')
        return code_list

    def _get_routine_ending(self, routine):
        """Returns the closing statements of the fortran routine."""
        if len(routine.results) == 1:
            return ['end function\n']
        else:
            return ['end subroutine\n']

    def get_interface(self, routine):
        """Returns a string for the function interface.

        The routine should have a single result object, which can be None.
        If the routine has multiple result objects, a CodeGenError is
        raised.

        See: https://en.wikipedia.org/wiki/Function_prototype

        """
        prototype = ['interface\n']
        prototype.extend(self._get_routine_opening(routine))
        prototype.extend(self._declare_arguments(routine))
        prototype.extend(self._get_routine_ending(routine))
        prototype.append('end interface\n')

        return ''.join(prototype)

    def _call_printer(self, routine):
        declarations = []
        code_lines = []
        for result in routine.result_variables:
            if isinstance(result, Result):
                assign_to = routine.name
            else:
                assign_to = result.result_var

            constants, not_fortran, f_expr = fcode(result.expr,
                                                   assign_to=assign_to, source_format='free', human=False)

            for obj, v in sorted(constants, key=str):
                t = get_default_datatype(obj)
                declarations.append(
                    f'{t.fname}, parameter :: {obj} = {v}\n')
            for obj in sorted(not_fortran, key=str):
                t = get_default_datatype(obj)
                if isinstance(obj, Function):
                    name = obj.func
                else:
                    name = obj
                declarations.append(f'{t.fname} :: {name}\n')

            code_lines.append(f'{f_expr}\n')
        return declarations + code_lines

    def _indent_code(self, codelines):
        p = FCodePrinter({'source_format': 'free', 'human': False})
        return p.indent_code(codelines)

    def dump_f95(self, routines, f, prefix, header=True, empty=True):
        # check that symbols are unique with ignorecase
        for r in routines:
            lowercase = {str(x).lower() for x in r.variables}
            orig_case = {str(x) for x in r.variables}
            if len(lowercase) < len(orig_case):
                raise CodeGenError(f"Fortran ignores case. Got symbols: {', '.join([str(var) for var in r.variables])}")
        self.dump_code(routines, f, prefix, header, empty)
    dump_f95.extension = code_extension  # type: ignore[attr-defined]
    dump_f95.__doc__ = CodeGen.dump_code.__doc__

    def dump_h(self, routines, f, prefix, header=True, empty=True):
        """Writes the interface to a header file.

        This file contains all the function declarations.

        Parameters
        ==========

        routines : list
            A list of Routine instances.

        f : file-like
            Where to write the file.

        prefix : string
            The filename prefix.

        header : bool, optional
            When True, a header comment is included on top of each source
            file.  [default : True]

        empty : bool, optional
            When True, empty lines are included to structure the source
            files.  [default : True]

        """
        if header:
            print(''.join(self._get_header()), file=f)
        if empty:
            print(file=f)
        # declaration of the function prototypes
        for routine in routines:
            prototype = self.get_interface(routine)
            f.write(prototype)
        if empty:
            print(file=f)
    dump_h.extension = interface_extension  # type: ignore[attr-defined]

    # This list of dump functions is used by CodeGen.write to know which dump
    # functions it has to call.
    dump_fns = [dump_f95, dump_h]


class OctaveCodeGen(CodeGen):
    """Generator for Octave code.

    The .write() method inherited from CodeGen will output a code file
    <prefix>.m.

    Octave .m files usually contain one function.  That function name should
    match the filename (``prefix``).  If you pass multiple ``name_expr`` pairs,
    the latter ones are presumed to be private functions accessed by the
    primary function.

    You should only pass inputs to ``argument_sequence``: outputs are ordered
    according to their order in ``name_expr``.

    """

    code_extension = 'm'

    def routine(self, name, expr, argument_sequence, global_vars=None):
        """Specialized Routine creation for Octave."""
        # FIXME: this is probably general enough for other high-level
        # languages, perhaps its the C/Fortran one that is specialized!

        if is_sequence(expr) and not isinstance(expr, (MatrixBase, MatrixExpr)):
            if not expr:
                raise ValueError('No expression given')
            expressions = Tuple(*expr)
        else:
            expressions = Tuple(expr)

        # local variables
        local_vars = {i.label for i in expressions.atoms(Idx)}

        # global variables
        global_vars = set() if global_vars is None else set(global_vars)

        # symbols that should be arguments
        symbols = expressions.free_symbols - local_vars - global_vars

        # Octave supports multiple return values
        return_vals = []
        for (i, expr) in enumerate(expressions):
            if isinstance(expr, Equality):
                out_arg = expr.lhs
                expr = expr.rhs
                symbol = out_arg
                if isinstance(out_arg, Indexed):
                    symbol = out_arg.base.label
                if not isinstance(out_arg, (Indexed, Symbol, MatrixSymbol)):
                    raise CodeGenError('Only Indexed, Symbol, or MatrixSymbol '
                                       'can define output arguments.')

                return_vals.append(Result(expr, name=symbol, result_var=out_arg))
                if not expr.has(symbol):
                    # this is a pure output: remove from the symbols list, so
                    # it doesn't become an input.
                    symbols.remove(symbol)

            else:
                # we have no name for this output
                return_vals.append(Result(expr, name=f'out{i + 1}'))

        # setup input argument list
        arg_list = []
        array_symbols = {}
        for array in expressions.atoms(Indexed):
            array_symbols[array.base.label] = array
        for array in expressions.atoms(MatrixSymbol):
            array_symbols[array] = array

        for symbol in sorted(symbols, key=str):
            arg_list.append(InputArgument(symbol))

        if argument_sequence is not None:
            # if the user has supplied IndexedBase instances, we'll accept that
            new_sequence = []
            for arg in argument_sequence:
                if isinstance(arg, IndexedBase):
                    new_sequence.append(arg.label)
                else:
                    new_sequence.append(arg)
            argument_sequence = new_sequence

            missing = [x for x in arg_list if x.name not in argument_sequence]
            if missing:
                msg = "Argument list didn't specify: {0} "
                msg = msg.format(', '.join([str(m.name) for m in missing]))
                raise CodeGenArgumentListError(msg, missing)

            # create redundant arguments to produce the requested sequence
            name_arg_dict = {x.name: x for x in arg_list}
            new_args = []
            for symbol in argument_sequence:
                try:
                    new_args.append(name_arg_dict[symbol])
                except KeyError:
                    new_args.append(InputArgument(symbol))
            arg_list = new_args

        return Routine(name, arg_list, return_vals, local_vars, global_vars)

    def _get_symbol(self, s):
        """Print the symbol appropriately."""
        return octave_code(s).strip()

    def _get_header(self):
        """Writes a common header for the generated files."""
        code_lines = []
        tmp = header_comment % {'version': diofant_version,
                                'project': self.project}
        for line in tmp.splitlines():
            if line == '':
                code_lines.append('%\n')
            else:
                code_lines.append(f'%   {line}\n')
        return code_lines

    def _preprocessor_statements(self, prefix):
        return []

    def _get_routine_opening(self, routine):
        """Returns the opening statements of the routine."""
        code_list = []
        code_list.append('function ')

        # Outputs
        outs = []
        for result in routine.results:
            # Note: name not result_var; want `y` not `y(i)` for Indexed
            s = self._get_symbol(result.name)
            outs.append(s)
        if len(outs) > 1:
            code_list.append('[' + (', '.join(outs)) + ']')
        else:
            code_list.append(''.join(outs))
        code_list.append(' = ')

        # Inputs
        args = []
        for arg in routine.arguments:
            if isinstance(arg, (OutputArgument, InOutArgument)):
                raise CodeGenError(f'Octave: invalid argument of type {type(arg)!s}')
            args.append(f'{self._get_symbol(arg.name)}')
        args = ', '.join(args)
        code_list.append(f'{routine.name}({args})\n')
        code_list = [''.join(code_list)]

        return code_list

    def _declare_arguments(self, routine):
        return []

    def _declare_globals(self, routine):
        if not routine.global_vars:
            return []
        s = ' '.join(sorted(self._get_symbol(g) for g in routine.global_vars))
        return ['global ' + s + '\n']

    def _declare_locals(self, routine):
        return []

    def _get_routine_ending(self, routine):
        return ['end\n']

    def _call_printer(self, routine):
        declarations = []
        code_lines = []
        for result in routine.results:
            assign_to = result.result_var

            constants, not_supported, oct_expr = octave_code(result.expr,
                                                             assign_to=assign_to,
                                                             human=False, inline=False)

            for obj, v in sorted(constants, key=str):
                declarations.append(
                    f'  {obj} = {v};  % constant\n')
            for obj in sorted(not_supported, key=str):
                if isinstance(obj, Function):
                    name = obj.func
                else:
                    name = obj
                declarations.append(
                    f'  % unsupported: {name}\n')
            code_lines.append(f'{oct_expr}\n')
        return declarations + code_lines

    def _indent_code(self, codelines):
        # Note that indenting seems to happen twice, first
        # statement-by-statement by OctavePrinter then again here.
        p = OctaveCodePrinter({'human': False})
        return p.indent_code(codelines)

    def dump_m(self, routines, f, prefix, header=True, empty=True, inline=True):
        # Note used to call self.dump_code() but we need more control for header

        code_lines = self._preprocessor_statements(prefix)

        for i, routine in enumerate(routines):
            if i > 0:
                if empty:
                    code_lines.append('\n')
            code_lines.extend(self._get_routine_opening(routine))
            if i == 0:
                if routine.name != prefix:
                    raise ValueError('Octave function name should match prefix')
                if header:
                    code_lines.append('%' + prefix.upper() +
                                      '  Autogenerated by diofant\n')
                    code_lines.append(''.join(self._get_header()))
            code_lines.extend(self._declare_arguments(routine))
            code_lines.extend(self._declare_globals(routine))
            code_lines.extend(self._declare_locals(routine))
            if empty:
                code_lines.append('\n')
            code_lines.extend(self._call_printer(routine))
            if empty:
                code_lines.append('\n')
            code_lines.extend(self._get_routine_ending(routine))

        code_lines = self._indent_code(''.join(code_lines))

        if code_lines:
            f.write(code_lines)

    dump_m.extension = code_extension  # type: ignore[attr-defined]
    dump_m.__doc__ = CodeGen.dump_code.__doc__

    # This list of dump functions is used by CodeGen.write to know which dump
    # functions it has to call.
    dump_fns = [dump_m]


def get_code_generator(language, project):
    CodeGenClass = {'C': CCodeGen, 'F95': FCodeGen,
                    'OCTAVE': OctaveCodeGen}.get(language.upper())
    if CodeGenClass is None:
        raise ValueError(f"Language '{language}' is not supported.")
    return CodeGenClass(project)


#
# Friendly functions
#


def codegen(name_expr, language, prefix=None, project='project',
            to_files=False, header=True, empty=True, argument_sequence=None,
            global_vars=None):
    """Generate source code for expressions in a given language.

    Parameters
    ==========

    name_expr : tuple, or list of tuples
        A single (name, expression) tuple or a list of (name, expression)
        tuples.  Each tuple corresponds to a routine.  If the expression is
        an equality (an instance of class Equality) the left hand side is
        considered an output argument.  If expression is an iterable, then
        the routine will have multiple outputs.

    language : string
        A string that indicates the source code language.  This is case
        insensitive.  Currently, 'C', 'F95' and 'Octave' are supported.
        'Octave' generates code compatible with both Octave and Matlab.

    prefix : string, optional
        A prefix for the names of the files that contain the source code.
        Language-dependent suffixes will be appended.  If omitted, the name
        of the first name_expr tuple is used.

    project : string, optional
        A project name, used for making unique preprocessor instructions.
        [default: "project"]

    to_files : bool, optional
        When True, the code will be written to one or more files with the
        given prefix, otherwise strings with the names and contents of
        these files are returned. [default: False]

    header : bool, optional
        When True, a header is written on top of each source file.
        [default: True]

    empty : bool, optional
        When True, empty lines are used to structure the code.
        [default: True]

    argument_sequence : iterable, optional
        Sequence of arguments for the routine in a preferred order.  A
        CodeGenError is raised if required arguments are missing.
        Redundant arguments are used without warning.  If omitted,
        arguments will be ordered alphabetically, but with all input
        aguments first, and then output or in-out arguments.

    global_vars : iterable, optional
        Sequence of global variables used by the routine.  Variables
        listed here will not show up as function arguments.

    Examples
    ========

    >>> [(c_name, c_code), (h_name, c_header)] = codegen(
    ...     ('f', x+y*z), 'C', 'test', header=False, empty=False)
    >>> print(c_name)
    test.c
    >>> print(c_code)
    #include "test.h"
    #include <math.h>
    double f(double x, double y, double z) {
      double f_result;
      f_result = x + y*z;
      return f_result;
    }
    >>> print(h_name)
    test.h
    >>> print(c_header)
    #ifndef PROJECT__TEST__H
    #define PROJECT__TEST__H
    double f(double x, double y, double z);
    #endif

    Another example using Equality objects to give named outputs.  Here the
    filename (prefix) is taken from the first (name, expr) pair.

    >>> from diofant.abc import f, g
    >>> [(c_name, c_code),
    ...  (h_name, c_header)] = codegen([('myfcn', x + y),
    ...                                 ('fcn2', [Eq(f, 2*x), Eq(g, y)])],
    ...                                'C', header=False, empty=False)
    >>> print(c_name)
    myfcn.c
    >>> print(c_code)
    #include "myfcn.h"
    #include <math.h>
    double myfcn(double x, double y) {
       double myfcn_result;
       myfcn_result = x + y;
       return myfcn_result;
    }
    void fcn2(double x, double y, double *f, double *g) {
       (*f) = 2*x;
       (*g) = y;
    }

    If the generated function(s) will be part of a larger project where various
    global variables have been defined, the 'global_vars' option can be used
    to remove the specified variables from the function signature

    >>> [(f_name, f_code), header] = codegen(
    ...     ('f', x+y*z), 'F95', header=False, empty=False,
    ...     argument_sequence=(x, y), global_vars=(z,))
    >>> print(f_code)
    REAL*8 function f(x, y)
    implicit none
    REAL*8, intent(in) :: x
    REAL*8, intent(in) :: y
    f = x + y*z
    end function

    """
    # Initialize the code generator.
    code_gen = get_code_generator(language, project)

    if isinstance(name_expr[0], str):
        # single tuple is given, turn it into a singleton list with a tuple.
        name_expr = [name_expr]

    if prefix is None:
        prefix = name_expr[0][0]

    # Construct Routines appropriate for this code_gen from (name, expr) pairs.
    routines = []
    for name, expr in name_expr:
        routines.append(code_gen.routine(name, expr, argument_sequence,
                                         global_vars))

    # Write the code.
    return code_gen.write(routines, prefix, to_files, header, empty)


def make_routine(name, expr, argument_sequence=None,
                 global_vars=None, language='F95'):
    """A factory that makes an appropriate Routine from an expression.

    Parameters
    ==========

    name : string
        The name of this routine in the generated code.

    expr : expression or list/tuple of expressions
        A Diofant expression that the Routine instance will represent.  If
        given a list or tuple of expressions, the routine will be
        considered to have multiple return values and/or output arguments.

    argument_sequence : list or tuple, optional
        List arguments for the routine in a preferred order.  If omitted,
        the results are language dependent, for example, alphabetical order
        or in the same order as the given expressions.

    global_vars : iterable, optional
        Sequence of global variables used by the routine.  Variables
        listed here will not show up as function arguments.

    language : string, optional
        Specify a target language.  The Routine itself should be
        language-agnostic but the precise way one is created, error
        checking, etc depend on the language.  [default: "F95"].

    A decision about whether to use output arguments or return values is made
    depending on both the language and the particular mathematical expressions.
    For an expression of type Equality, the left hand side is typically made
    into an OutputArgument (or perhaps an InOutArgument if appropriate).
    Otherwise, typically, the calculated expression is made a return values of
    the routine.

    Examples
    ========

    >>> from diofant.abc import f, g
    >>> r = make_routine('test', [Eq(f, 2*x), Eq(g, x + y)])
    >>> [arg.result_var for arg in r.results]
    []
    >>> [arg.name for arg in r.arguments]
    [x, y, f, g]
    >>> [arg.name for arg in r.result_variables]
    [f, g]
    >>> r.local_vars
    set()

    Another more complicated example with a mixture of specified and
    automatically-assigned names.  Also has Matrix output.

    >>> r = make_routine('fcn', [x*y, Eq(f, 1), Eq(g, x + g), Matrix([[x, 2]])])
    >>> [arg.result_var for arg in r.results]
    [result_...]
    >>> [arg.expr for arg in r.results]
    [x*y]
    >>> [arg.name for arg in r.arguments]
    [x, y, f, g, out_...]

    We can examine the various arguments more closely:

    >>> [a.name for a in r.arguments if isinstance(a, InputArgument)]
    [x, y]

    >>> [a.name for a in r.arguments if isinstance(a, OutputArgument)]
    [f, out_...]
    >>> [a.expr for a in r.arguments if isinstance(a, OutputArgument)]
    [1, Matrix([[x, 2]])]

    >>> [a.name for a in r.arguments if isinstance(a, InOutArgument)]
    [g]
    >>> [a.expr for a in r.arguments if isinstance(a, InOutArgument)]
    [g + x]

    """
    # initialize a new code generator
    code_gen = get_code_generator(language, 'nothingElseMatters')

    return code_gen.routine(name, expr, argument_sequence, global_vars)
