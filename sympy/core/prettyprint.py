"""Prettyprinter idea by Jurjen Bos 2007

(I hate spammers: mail me at pietjepuk314 at the reverse of ku.oc.oohay).

Adapted for sympy by Ondrej Certik 2007.

"""

class StringPict:
	"""A ASCII picture.
	The pictures are represented as a list of equal length strings.
	"""
	#special value for StringPict.below
	LINE = 'line'

	def __init__(self, s, baseline=0):
		"""Initialize from string.
		Multiline strings are centered.
		"""
		#picture is a string that just can be printed
		self.picture = StringPict.equalLengths(s.splitlines())
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
			if isinstance(arg, str): arg = StringPict(arg)
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
		return StringPict('\n'.join(result),newBaseline)

	def right(self, *args):
		"""Put pictures next to this one.
		(Multiline) strings are allowed, and are given a baseline of 0.
		
		>>> StringPict("10").right("+",StringPict("1\r-\r2",1))
		   1
		10+-
		   2
		"""
		return StringPict.next(self, *args)

	def left(self, *args):
		"""Put pictures (left to right) at left.
		"""
		return StringPict.next(*(args+(self,)))

	@staticmethod
	def stack(*args):
		"""Put pictures on top of each other,
		from top to bottom.
		The baseline is the baseline of the second picture.
		Everything is centered.
		Baseline is the baseline of the second picture.
		Strings are allowed.
		The special value StringPict.LINE is a row of '-' extended to the width.
		"""
		#convert everything to stringPicts; keep LINE
		objects = []
		for arg in args:
			if arg is not StringPict.LINE and isinstance(arg, str):
				arg = StringPict(arg)
			objects.append(arg)

		#compute new width
		newWidth = max(
			obj.width()
			for obj in objects
			if obj is not StringPict.LINE)

		lineObj = StringPict('-'*newWidth)

		#replace LINE with proper lines
		for i, obj in enumerate(objects):
			if obj is StringPict.LINE:
				objects[i] = lineObj

		#stack the pictures, and center the result
		newPicture = []
		for obj in objects:
			newPicture.extend(obj.picture)
		newPicture = [line.center(newWidth) for line in newPicture]
		newBaseline = objects[0].height()+objects[1].baseline
		return StringPict('\n'.join(newPicture), newBaseline)

	def below(self, *args):
		"""Put pictures under this picture.
		Baseline is baseline of top picture
		>>> StringPict("x+3").below(StringPict.LINE, '3')
		x+3
		---
		 3
		"""
		result = StringPict.stack(self, *args)
		result.baseline = self.baseline
		return result

	def top(self, *args):
		"""Put pictures (top to bottom) at top.
		Baseline is baseline of bottom picture."""
		result = StringPict.stack(*(args+(self,)))
		result.baseline = result.height()-self.height()+self.baseline
		return result

	def parens(self):
		"""Put parentheses around self.
		"""
		height = self.height()
		if height==1:
			return StringPict('(').right(self, ')')
		else:
			verticalBar = '\n' + '|\n' * (self.height()-2)
			lparen = StringPict('/'+verticalBar+'\\',self.baseline)
			rparen = StringPict('\\'+verticalBar+'/',self.baseline)
			return lparen.right(" ", self, " ", rparen)

	def leftslash(self):
		"""Precede object by a slash of the proper size.
		"""
		height = max(
			self.baseline,
			self.height()-1-self.baseline)*2 + 1
		slash = '\n'.join(
			' '*(height-i-1)+'/'+' '*i
			for i in range(height)
			)
		return self.left(StringPict(slash, height//2))

	def root(self, n=None):
		"""Produce a nice root symbol.
		Produces ugly results for big n inserts.
		"""
		#put line over expression
		result = self.top('_'*self.width())
		#construct right half of root symbol
		height = self.height()
		slash = '\n'.join(
			' ' * (height-i-1) + '/' + ' ' * i
			for i in range(height)
			)
		slash = StringPict(slash, height-1)
		#left half of root symbol
		if height > 2:
			downline = StringPict('\\ \n \\',1)
		else:
			downline = StringPict('\\')
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
		return "StringPict(%r,%d)"%('\n'.join(self.picture), self.baseline)
