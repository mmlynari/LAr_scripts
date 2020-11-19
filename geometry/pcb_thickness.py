from math import sqrt, log # log is neper logarithm

# everything in microns
h_hv = 65 # default in CDR is 100
h_m = 242.5 # default in CDR is 285
h_s = 170
trace_thickness = 35
other_conductor_thickness = 35
w_t = 127
dielectric_constant = 4.6
total_thickness_no_trace_no_conductor = (h_hv + h_m + h_s) * 2
total_thickness_no_trace = total_thickness_no_trace_no_conductor + 2 * other_conductor_thickness * 3 # there are three conductor layers per half PCB
total_thickness = total_thickness_no_trace + trace_thickness # there are three conductor layers per half PCB
print "Total PCB thickness: %f mm\n"%(total_thickness/1000.0)

# check that all dimensions comply
# avoid spark between HV layer and signal pads (H_hv)
electrical_rigidity_FR4 = 20 #kV/mm    electric field that it cain sustain without creating shortcircuit/
wished_electric_field_in_gap = 1 # kV/mm
LAr_gap_thickness = 1.24 # mm
potential_hv_needed = wished_electric_field_in_gap * LAr_gap_thickness # kV
minimum_space_hv_sigpad = potential_hv_needed / electrical_rigidity_FR4
print "Mimimum spacing between signal pads and HV plate: %d microns"%(minimum_space_hv_sigpad*1000)
print "Foreseen space: %d microns\n"%(h_hv)

# check the transmission line impedence
def get_impedence(space_between_inner_side_ground, trace_width, trace_thickness, dielectric_constant):
    return 60 * log(1.9*space_between_inner_side_ground / (0.8 * trace_width + trace_thickness))  / sqrt(dielectric_constant) # from martin's talk
space_between_inner_side_ground = h_s * 2 + trace_thickness 
impedence = get_impedence(space_between_inner_side_ground, w_t, trace_thickness, dielectric_constant)
print "Transmission line impedence for h_s=%d, w_t=%d, trace_thickness=%d and dielectric constant=%f: %f Ohm"%(h_s, w_t, trace_thickness, dielectric_constant, impedence)
