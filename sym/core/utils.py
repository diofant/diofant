def sign(x):
    """Return the sign of x, that is, 
    1 if x is positive, 0 if x == 0 and -1 if x is negative
    """
    if x < 0: return -1
    elif x==0: return 0
    else: return 1
    
def isnumber(x):
    """Return True if x is a number. 
    """
    from numbers import Number
    return isinstance(x, (Number, int, float, long))
