def frange(a_min, a_max, a_steps):
    """
    Helper function which returns an array containing a_step+1 elements
    ranging from a_min to a_max.
    """
    a_delta = (a_max-a_min)/float(a_steps)
    return [a_min+a_delta*i for i in range(a_steps+1)]

def interpolate(a_min, a_max, a_ratio):
    return a_min + a_ratio * (a_max - a_min)

def rinterpolate(a_min, a_max, a_value):
    a_range = a_max-a_min
    if a_range == 0:
        a_range = 1.0
    return (a_value - a_min) / float(a_range)

def interpolatecolor(color1, color2, ratio):
    return [interpolate(color1[i], color2[i], ratio) for i in range(3)]