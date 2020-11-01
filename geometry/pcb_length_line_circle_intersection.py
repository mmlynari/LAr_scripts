from math import tan, radians, sqrt

# inner_radius is the radius at starting point of electrode, outer radius is the radius at end point of electrod, angle in degrees
def get_intersection_line_circle(inner_radius, outer_radius, inclination):
    # line equation y=ax+b
    a = tan(radians(inclination))
    b = -1 * tan(radians(inclination)) * inner_radius

    # resolution of circle to line intersection equation
    # from wolfram
    #x1 = (-sqrt(a**2 * outer_radius**2 - b**2 + outer_radius**2) - a * b)/(a**2 + 1)
    #y1 = (b - a * sqrt(a**2 * outer_radius**2 - b**2 + outer_radius**2))/(a**2 + 1)
    #x2 = (sqrt(a**2 * outer_radius**2 - b**2 + outer_radius**2) - a * b)/(a**2 + 1)
    #y2 = (a * sqrt(a**2 * outer_radius**2 - b**2 + outer_radius**2) + b)/(a**2 + 1)

    # done by hand
    rho = (2 * a * b) ** 2 - 4 * (a ** 2 + 1) * (b ** 2 - outer_radius ** 2)
    x1 = (-2 * a * b + sqrt(rho)) / (2 * (a ** 2 +1))
    x2 = (-2 * a * b - sqrt(rho)) / (2 * (a ** 2 +1))
    y1 = a * x1 + b
    y2 = a * x2 + b
    print "The circle and the line intersect in (%f, %f) and (%f, %f)..."%(x1, y1, x2, y2)
    p1x = inner_radius
    p1y = 0
    # take the positive solution 
    if (x1 > 0 and y1 > 0):
        p2x = x1
        p2y = y1
    else:
        p2x = x2
        p2y = y2
    length = sqrt((p2x - p1x)**2 + (p2y - p1y)**2)
    print "Length of the readout parallel to it is %f"%length

print "get_intersection_line_circle(216, 256, 50)"
get_intersection_line_circle(216, 256, 50)

# presampler 1.5 regular cells 3.5
print "get_intersection_line_circle(216, 217.5, 50)"
get_intersection_line_circle(216, 217.5, 50)

print "get_intersection_line_circle(216, 217.5, 50)"
get_intersection_line_circle(217.5, 221, 50)

