def sign(x):
    if x < 0: return -1
    elif x==0: return 0
    else: return 1
    
def isnumber(x):
    """Return true if x is a number. 
    """
    from numbers import Number
    return isinstance(x, (Number, int, float, long))
