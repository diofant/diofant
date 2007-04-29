
from sympy.core.basic import Basic


	
def pretty(expr):
	"""Returns the pretty representation for expr (as a string)
	"""
	if hasattr(expr, '__pretty__'):
		return str( expr.__pretty__() )
	else:
		return str(expr)

def pretty_print(expr):
	"""
	Prints expr in pretty form. 
	
	pprint is just a shortcut for this function
	"""
	print pretty(expr)
	
pprint = pretty_print