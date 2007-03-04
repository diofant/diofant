def sign(x):
    if x < 0: return -1
    elif x==0: return 0
    else: return 1
    
def isnumber(x):
    """Return true if x is a number. 
    """
    # we suppose that if a class has the __float__ method, then it contains
    # a numberic data type
    return hasattr(x, '__float__')
