"""
I was missing a decent pretty printer, so I wrote one myself.
Although it is a separate file, it is supposed to be sprinkled over
the code, but maybe you want to do it another way.
Anyway, the code is structured so that different operators all need a
little code to support pretty printing.
You may want to use this code however you like, or throw it away if
you have something better.

Jurjen N.E. Bos
Deukelven 48
1963 SR Heemskerk
0251-252426
"""
import sys, operator

test = True
class stringPict:
	"""A ASCII picture.
	The pictures are represented as a list of equal length strings.
	"""
	#special value for stringPict.below
	LINE = 'line'

	def __init__(self, s, baseline=0):
		"""Initialize from string.
		Multiline strings are centered.
		"""
		#picture is a string that just can be printed
		self.picture = stringPict.equalLengths(s.splitlines())
		#baseline is the line number of the "base line"
		self.baseline = baseline

	@staticmethod
	def equalLengths(lines):
		width = max(len(line) for line in lines)
		return [line.center(width) for line in lines]

	def height(self):
		return len(self.picture)

	def width(self):
		return len(self.picture[0])

	@staticmethod
	def next(*args):
		"""Put a string of stringPicts next to each other.
		"""
		#convert everything to stringPicts
		objects = []
		for arg in args:
			if isinstance(arg, str): arg = stringPict(arg)
			objects.append(arg)

		#make a list of pictures, with equal height and baseline
		newBaseline = max(obj.baseline for obj in objects)
		newHeightBelowBaseline = max(
			obj.height()-obj.baseline
			for obj in objects)
		newHeight = newBaseline + newHeightBelowBaseline

		pictures = []
		for obj in objects:
			oneEmptyLine = [' '*obj.width()]
			basePadding = newBaseline-obj.baseline
			totalPadding = newHeight-obj.height()
			pictures.append(
				oneEmptyLine * basePadding +
				obj.picture +
				oneEmptyLine * (totalPadding-basePadding))

		result = [''.join(lines) for lines in zip(*pictures)]
		return stringPict('\n'.join(result),newBaseline)

	def right(self, *args):
		"""Put pictures next to this one.
		(Multiline) strings are allowed, and are given a baseline of 0.
		>>> stringPicture("10").right("+",stringPict("1\r-\r2",1))
		   1
		10+-
		   2
		"""
		return stringPict.next(self, *args)

	def left(self, *args):
		"""Put pictures (left to right) at left.
		"""
		return stringPict.next(*(args+(self,)))

	@staticmethod
	def stack(*args):
		"""Put pictures on top of each other,
		from top to bottom.
		The baseline is the baseline of the second picture.
		Everything is centered.
		Baseline is the baseline of the second picture.
		Strings are allowed.
		The special value stringPict.LINE is a row of '-' extended to the width.
		"""
		#convert everything to stringPicts; keep LINE
		objects = []
		for arg in args:
			if arg is not stringPict.LINE and isinstance(arg, str):
				arg = stringPict(arg)
			objects.append(arg)

		#compute new width
		newWidth = max(
			obj.width()
			for obj in objects
			if obj is not stringPict.LINE)

		lineObj = stringPict('-'*newWidth)

		#replace LINE with proper lines
		for i, obj in enumerate(objects):
			if obj is stringPict.LINE:
				objects[i] = lineObj

		#stack the pictures, and center the result
		newPicture = []
		for obj in objects:
			newPicture.extend(obj.picture)
		newPicture = [line.center(newWidth) for line in newPicture]
		newBaseline = objects[0].height()+objects[1].baseline
		return stringPict('\n'.join(newPicture), newBaseline)

	def below(self, *args):
		"""Put pictures under this picture.
		>>> stringPict("x+3").below(stringPict.LINE, '3')
		x+3
		---
		 3
		"""
		return stringPict.stack(self, *args)

	def top(self, *args):
		"""Put pictures (top to bottom) at top."""
		return stringPict.stack(*(args+(self,)))

	def parens(self):
		"""Put parentheses around self.
		"""
		height = self.height()
		if height==1:
			return stringPict('(').right(self, ')')
		else:
			verticalBar = '\n'+'|\n'*(self.height()-2)
			lparen = stringPict('/'+verticalBar+'\\',self.baseline)
			rparen = stringPict('\\'+verticalBar+'/',self.baseline)
			return lparen.right(self, rparen)

	def leftslash(self):
		"""Precede object by a slash of the proper size.
		"""
		height = max(
			self.baseline,
			self.height()-1-self.baseline)*2+1
		slash = '\n'.join(
			' '*(height-i-1)+'/'+' '*i
			for i in range(height)
			)
		return self.left(stringPict(slash, height//2))

	def root(self, n=None):
		"""Produce a nice root symbol.
		Produces ugly results for big n inserts.
		"""
		#put line over expression
		result = self.top(stringPict.LINE)
		#construct right half of root symbol
		height = self.height()
		slash = '\n'.join(
			' '*(height-i-1)+'/'+' '*i
			for i in range(height)
			)
		slash = stringPict(slash, height-1)
		#left half of root symbol
		if height>2:
			downline = stringPict('\\ \n \\',1)
		else:
			downline = stringPict('\\')
		#put n on top, as low as possible
		if n is not None and n.width()>downline.width():
			downline = downline.left(' '*(n.width()-downline.width()))
			downline = downline.top(n)
		#build root symbol
		root = downline.right(slash)
		#glue it on at the proper height
		#normally, the root symbel is as high as self
		#which is one less than result
		#this moves the root symbol one down
		#if the root became higher, the baseline has to grow too
		root.baseline = result.baseline-result.height()+root.height()
		return result.left(root)

	def __str__(self):
		return '\n'.join(self.picture)

	def __repr__(self):
		return "stringPict(%r,%d)"%('\n'.join(self.picture), self.baseline)

if test:
	print stringPict("10").right("+",stringPict("1\r-\r2",1))
	print stringPict("x+3").below(stringPict.LINE, '3')
	print stringPict("1\r-\r2",1).parens().right('+','10')

class expression:
	"""A simple class to demonstrate the features of stringPict."""
	def __init__(self, *args):
		"""Apply function to args."""
		self.function = self.__class__.__name__
		self.args = args

	def __str__(self):
		return str(self.pretty())

	def pretty(self):
		"Override this for your class: should return a stringPict."

	def __repr__(self):
		return "%s(%r)"%(self.function, ','.join(self.args))

class atom(expression):
	def __init__(self, value):
		self.function = 'atom'
		self.args = [value]
	def pretty(self):
		return stringPict(self.args[0])

class neg(expression):
	def __init__(self, value):
		self.function = 'neg'
		self.args = [value]
	def pretty(self):
		return self.args[0].pretty().left('-')

class add(expression):
	def pretty(self):
		"""Make a pretty addition.
		Parentheses are not needed here.
		Addition of negative numbers is simplified.
		"""
		result = []
		for arg in self.args:
			if isinstance(arg, neg):
				result.append('-')
				arg = arg.args[0]
			elif result:
				result.append('+')
			result.append(arg.pretty())
		return stringPict.next(*result)

class sub(expression):
	def __init__(self, a, b):
		"""a - b"""
		self.function = 'sub'
		self.args = [a, b]
	def pretty(self):
		"""Make a pretty subtraction.
		Parentheses go around the second argument if it is + or -.
		"""
		[a, b] = self.args
		if isinstance(b, neg):
			"""Replace -- by +, saving parentheses."""
			return stringPict.next(a.pretty(), '+', b.args[0].pretty())
		bpretty = b.pretty()
		if isinstance(b, (add, sub)): bpretty = bpretty.parens()
		return stringPict.next(a.pretty(), '-', bpretty)

class inv(expression):
	def __init__(self, value):
		self.function = 'inv'
		self.args = [value]
	def pretty(self):
		return self.args[0].pretty().top('1', stringPict.LINE)

class mul(expression):
	def pretty(self):
		"""Make a pretty addition.
		Parentheses are needed around +, - and neg.
		"""
		result = []
		for arg in self.args:
			#todo: Multiplication by inverse simplified to a slash fraction.
			if isinstance(arg, inv):
				if not result: result.append(stringPict('1'))
				argpretty = arg.pretty()
				if isinstance(arg, (add, neg, sub, mul)):
					argpretty = argpretty.parens()
				result.append(argpretty.leftslash())
				continue
			if result:
				result.append('*')
			argpretty = arg.pretty()
			if isinstance(arg, (add, sub, neg)):
				argpretty = argpretty.parens()
			result.append(argpretty)
		return stringPict.next(*result)

class div(expression):
	def __init__(self, a, b):
		"""a / b"""
		self.function = 'div'
		self.args = [a, b]
	def pretty(self):
		"""Make a pretty division.
		Parentheses aren't needed.
		"""
		[a, b] = self.args
		if isinstance(b, inv):
			"""Replace /(1/) by *, saving space."""
			return stringPict.next(a.pretty(), '*', b.args[0].pretty())
		return stringPict.stack(a.pretty(), stringPict.LINE, b.pretty())

class power(expression):
	def __init__(self, a, b):
		"""a ** b"""
		self.function = 'power'
		self.args = [a, b]
	def pretty(self):
		"""Make a pretty power.
		Parentheses aren't needed.
		"""
		[a, b] = self.args
		apretty = a.pretty()
		if isinstance(a, (neg, add, sub, inv, mul, div, power)):
			apretty = apretty.parens()
		bpretty = b.pretty()
		exponent = bpretty.left(' '*apretty.width())
		apretty = apretty.right(' '*bpretty.width())
		return apretty.top(exponent)

class root(expression):
	def __init__(self, value):
		self.function = 'root'
		self.args = [value]
	def pretty(self):
		return self.args[0].pretty().root()

class xroot(expression):
	def __init__(self, a, b):
		"""b th root of a"""
		self.function = 'xroot'
		self.args = [a, b]
	def pretty(self):
		"""Make a pretty x-th root.
		Parentheses aren't needed, but it sometimes is ugly.
		"""
		[a, b] = self.args
		return a.pretty().root(b.pretty())

if test:
	x = atom('x')
	one = atom('1')
	two = atom('2')
	print x
	print neg(x)
	print add(one,x,neg(x))
	print sub(x,one)
	print sub(two, neg(x))
	print inv(add(x,one))
	print mul(inv(x),inv(x),two) #gaat mis
	print div(one, two)
	print div(two, inv(x))
	print power(add(x,one),power(two,x))
	print root(add(power(x, two), two))
	print xroot(one, power(x, two))
	print root(div(one, two))
