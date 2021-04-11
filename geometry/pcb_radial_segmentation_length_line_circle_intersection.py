from math import tan, radians, sqrt, sin, asin, degrees

# inner_radius is the radius at starting point of electrode, outer radius is the radius at end point of electrod, angle in degrees
def get_cell_length_from_intersection_line_circle(inner_radius, outer_radius, inclination):
    # line equation y=ax+b, assuming the segment starts with touching the x axis
    a = tan(radians(inclination))
    b = -1 * tan(radians(inclination)) * inner_radius

    # resolution of circle to line intersection equation
    # from wolfram
    x1 = (-sqrt(a**2 * outer_radius**2 - b**2 + outer_radius**2) - a * b)/(a**2 + 1)
    y1 = (b - a * sqrt(a**2 * outer_radius**2 - b**2 + outer_radius**2))/(a**2 + 1)
    x2 = (sqrt(a**2 * outer_radius**2 - b**2 + outer_radius**2) - a * b)/(a**2 + 1)
    y2 = (a * sqrt(a**2 * outer_radius**2 - b**2 + outer_radius**2) + b)/(a**2 + 1)

    # done by hand
    rho = (2 * a * b) ** 2 - 4 * (a ** 2 + 1) * (b ** 2 - outer_radius ** 2)
    #x1 = (-2 * a * b + sqrt(rho)) / (2 * (a ** 2 +1))
    #x2 = (-2 * a * b - sqrt(rho)) / (2 * (a ** 2 +1))
    #y1 = a * x1 + b
    #y2 = a * x2 + b

    #print "The circle and the line intersect in (%f, %f) and (%f, %f)..."%(x1, y1, x2, y2)
    # choose easy coordinate, out the segment start on the x axis
    p1x = inner_radius
    p1y = 0
    # take the correct solution 
    if (x1 > 0 and y1 > 0):
        p2x = x1
        p2y = y1
    else:
        p2x = x2
        p2y = y2
    #print "Chosen solution: %f %f"%(p2x, p2y)
    #print "Radius - x solution %f"%(p2x-inner_radius)
    #print "y0 - y solution %f"%(p2y-p1y)
    length = sqrt((p2x - p1x)**2 + (p2y - p1y)**2)
    #print "Length of the readout parallel to it is %f"%length
    return length

inner_radius = 216
inclination = 50 
first_cell_radial_length = 1.5
other_cells_radial_lengths = 3.5
n_cell_after_ps = 11
outer_radius = inner_radius + first_cell_radial_length + n_cell_after_ps*other_cells_radial_lengths
calo_depth = outer_radius - inner_radius


total_pcb_length_parallel = get_cell_length_from_intersection_line_circle(inner_radius, outer_radius, inclination)
print "Total length from %d to %d with %d inclination: "%(inner_radius, outer_radius, inclination), total_pcb_length_parallel
print "Pre sampler radial length %f, other cell radial length %f, number of other cells %d."%(first_cell_radial_length, other_cells_radial_lengths, n_cell_after_ps)
pre_sampler_length_parallel = get_cell_length_from_intersection_line_circle(inner_radius, inner_radius+first_cell_radial_length, inclination)

cell_lengths_along_pcb_direction = []
cell_lengths_along_pcb_direction.append(pre_sampler_length_parallel)

regular_cells_radial_segementation = [inner_radius + first_cell_radial_length + i * other_cells_radial_lengths for i in range(0, n_cell_after_ps+1)]
print "Radius of separation for the regular cells: ", regular_cells_radial_segementation

for index in range(len(regular_cells_radial_segementation)-1):
    cell_lengths_along_pcb_direction.append(get_cell_length_from_intersection_line_circle(regular_cells_radial_segementation[index], regular_cells_radial_segementation[index+1], inclination))

print "All cell lengths: ", cell_lengths_along_pcb_direction
print "Total length from %d to %d  with %d inclination, by summing all cells: "%(inner_radius, outer_radius, inclination), sum(cell_lengths_along_pcb_direction)
print "Looks like a numerical precision issue?"

print "Length of the pre-sampler cell: %f"%(pre_sampler_length_parallel)
print "Length of the cells after pre-shower to reach total radius: %f"%((total_pcb_length_parallel-pre_sampler_length_parallel)/float(n_cell_after_ps))
print "--------------------------------"

ratio_first_cell = first_cell_radial_length/calo_depth
ratio_other_cells = other_cells_radial_lengths/calo_depth
pre_sampler_length_parallel = total_pcb_length_parallel * ratio_first_cell
normal_cells_length_parallel = total_pcb_length_parallel * ratio_other_cells

print "Length of the pre-sampler cell parallel: %f"%(pre_sampler_length_parallel)
print "Length of the other cells parallel: %f"%(normal_cells_length_parallel)
print "Total lenght parallel: %f"%(pre_sampler_length_parallel + n_cell_after_ps * normal_cells_length_parallel)

# angle comprise by lines from  1) IP to inner right edge of a cell, 2) IP to outer left edge of a cell (useful to get the plate angle with radial direction that changes with increasing R)
# based on scalene triangle sine law A/sin(a) = B/sin(b) = C/sin(c) (outer left edge aligned on the Y axis)
delta_phi = degrees(asin(total_pcb_length_parallel/outer_radius * sin(radians(inclination))))
print "Angle sustained by the readout projection: %f"%delta_phi

radial_inclination_end = inclination - delta_phi


