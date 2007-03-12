#!/usr/bin/env python
"""
Python shell for sympy. Imports sympy and defines the symbols x, y, z. 

command line options: 
  -c : permits to specify a python interactive interpreter, currently only 
       python or ipython. Example usage: 
           ./shell.py -c python
       default is set to ipython
       
  -h : prints this help message
"""

welcome_msg = """
Python console for sympy.

Modules already imported: sympy
Symbols defined: x, y, z
"""

welcome_msg_ipython = """
Python console for sympy.

Modules already imported: sympy
Symbols defined: x, y, z

?       -> Introduction to IPython's features.
%magic  -> Information about IPython's 'magic' % functions.
help    -> Python's own help system.
object? -> Details about 'object'. ?object also works, ?? prints more.
"""

import sys
try:
    import sympy
    from sympy import *
    # leave both so that the user can do help(sympy) and the like
    
except ImportError:
    print "Could not find sympy\n...exiting"
    sys.exit()


import getopt

def run_ipython_interpreter():
	
	from IPython.Shell import IPShellEmbed
	
	ipshell = IPShellEmbed()

	x = Symbol('x')
	y = Symbol('y')
	z = Symbol('z')

	# Now start an embedded ipython.
	ipshell(welcome_msg_ipython)
	sys.exit("Exiting ...")

def run_python_interpreter():
	print """
Couldn't locate IPython. Having IPython installed is greatly recomended.	
See http://ipython.scipy.org for more details. If you use Debian, just install
the "ipython" package and start isym.py again.\n"""
	
	import code
	import readline
	import atexit
	import os
	
	
	class HistoryConsole(code.InteractiveConsole):
	    def __init__(self, locals=None, filename="<console>",
	                 histfile=os.path.expanduser("~/.sympy-history")):
	        code.InteractiveConsole.__init__(self)
	        self.init_history(histfile)
	
	    def init_history(self, histfile):
	        readline.parse_and_bind("tab: complete")
	        if hasattr(readline, "read_history_file"):
	            try:
	                readline.read_history_file(histfile)
	            except IOError:
	                pass
	            atexit.register(self.save_history, histfile)
	
	    def save_history(self, histfile):
	        readline.write_history_file(histfile)
		
	sh = HistoryConsole()
	sh.runcode("from sympy import *")
	sh.runcode("x = Symbol('x')")
	sh.runcode("y = Symbol('y')")
	sh.runcode("z = Symbol('z')")
	sh.interact(welcome_msg)
	sys.exit("Exiting ...")



argv = sys.argv[1:]

optlist, rest = getopt.getopt(argv, "c:h")

for opt in optlist:
	if opt[0] == '-c':
		if opt[1] == 'python':
			run_python_interpreter()
		elif opt[1] == 'ipython':
			run_ipython_interpreter()
	elif opt[0] == '-h':
		print __doc__
		sys.exit()

try:
	run_ipython_interpreter()
	
except ImportError:
	run_python_interpreter()

