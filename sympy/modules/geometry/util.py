def intersection(e1, e2):
    """
    Finds the intersection between two GeometryEntity instances. Returns a list
    of all the intersections, or None if there is none. Will raise a
    NotImplementedError exception if unable to calculate the intersection.
    """
    from entity import GeometryEntity
    try:
        return GeometryEntity.do_intersection(e1, e2)
    except NotImplementedError:
        n1 = e1.__class__.__name__
        n2 = e2.__class__.__name__
        raise NotImplementedError("Cannot find intersection between %s and %s" % (n1, n2))

def convex_hull(*args):
    """
    Returns a Polygon representing the convex hull of a set of 2D points.
    Can only be performed on a set of non-symbolic points.
    """
    from point import Point
    from line import Segment
    from polygon import Polygon

    p = args[0]
    if isinstance(args[0], Point):
        p = args

    def tarea(a, b, c):
        return (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])

    # Basic checks
    if len(p) == 1:
        return p[0]
    elif len(p) == 2:
        return Segment(p[0], p[1])

    # Find lowest+rightmost point
    m = 0
    for i in xrange(1, len(p)):
        if (p[i][1] < p[m][1]) or ((p[i][1] == p[m][1]) and (p[i][0] > p[m][0])):
            m = i
    p[0], p[m] = p[m], p[0]

    # Sort points
    destroy = {}
    p0 = p[0]
    def pcompare(p1, p2):
        a = tarea(p0, p1, p2)
        if a > 0:
            return -1
        elif a < 0:
            return 1
        else:
            x = abs(p1[0] - p0[0]) - abs(p2[0] - p0[0])
            y = abs(p1[1] - p0[1]) - abs(p2[1] - p0[1])
            if (x < 0) or (y < 0):
                destroy[p1] = True
                return -1
            elif (x > 0) or (y > 0):
                destroy[p2] = True
                return 1
            else:
                destroy[p1] = True
                return 0
    p = p[1:]
    p.sort(pcompare)
    p.insert(0, p0)

    # Destroy points as found by sorting
    for i in xrange(len(p)-1, -1, -1):
        if p[i] in destroy:
            del p[i]

    # Graham scan
    def isleft(a, b, c):
        return (tarea(a, b, c) > 0)

    top = [p[0], p[1]]
    i = 2
    while i < len(p):
        p1 = top[-2]
        p2 = top[-1]
        if isleft(p1, p2, p[i]):
            top.append(p[i])
            i += 1
        else:
            top.pop()
    return Polygon(top)