def sign(x):
    """Return the sign of x, that is, 
    1 if x is positive, 0 if x == 0 and -1 if x is negative
    """
    if x < 0: return -1
    elif x==0: return 0
    else: return 1
    
def isnumber(x):
    """Return True if x is a number. False otherwise. 
    """
    
    from numbers import Number
    from basic import Basic
    from decimal import Decimal
    
    if isinstance(x, (Number, int, float, long, Decimal)):
        return True
    elif isinstance(x, Basic):
        try:
            x.evalf()
            # if x has symbols it will raise
            # an exception, which we will catch
            # and return False
            return True
        except ValueError:
            return False
    return False
