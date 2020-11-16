from math import tan, log, atan, exp, degrees, radians
from numpy import arange

#everything in degrees
def get_theta(eta):
    theta = 2 * atan(exp(-eta))
    return degrees(theta)
def get_eta(theta):
    eta = -log(tan(radians(theta/2)))
    return eta

calo_half_dz = 220
calo_rMin = 216
calo_rMax = 256

eta_seg = arange(0, 0.89, 0.01)
#print "Eta: ", eta_seg
#print "\n"
theta_seg = [get_theta(eta) for eta in eta_seg]
#print "Theta: ", [round(a, 2) for a in theta_seg]
#print "\n"
theta_delta = []
cell_width_inner = []
cell_width_outer = []
total_width_sofar_inner = 0
total_width_sofar_outer = 0
total_angle_sofar = 0
for i in range(0, eta_seg.size - 1):
    delta_theta = theta_seg[i] - theta_seg[i+1]
    theta_delta.append(delta_theta)
    total_angle_sofar += delta_theta
    li = tan(radians(total_angle_sofar)) * calo_rMin - total_width_sofar_inner
    total_width_sofar_inner += li
    lo = tan(radians(total_angle_sofar)) * calo_rMax - total_width_sofar_outer
    total_width_sofar_outer += lo

    cell_width_inner.append(li)
    cell_width_outer.append(lo)

##print "Delta Theta: ", [round(a, 3) for a in theta_delta]
#print "Delta Theta: ", [a for a in theta_delta]
#print "\n"
#print "Inner radius widths: ", [round(a, 2) for a in cell_width_inner]
#print "\n"
#print "Outer radius widths: ", [round(a, 2) for a in cell_width_outer]
#print "\n"
##Derive angle covering based on calo length
#covering_angle = degrees(atan(calo_half_dz/calo_rMin))
#print "For the calo dimension given, covering angle is up to eta ", get_eta(covering_angle)
#print "Covering angle from both sides of radial direction: ", covering_angle
#print "We cover from %f to %f theta degrees"%(90-covering_angle, 90+covering_angle)
#print "Number of eta cells: %d * 2"%len(theta_delta)
#print "--------------------------------------"

calo_half_dz = 220
calo_rMin = 216
calo_rMax = 256
covering_angle = 45 #degrees(atan(calo_half_dz/calo_rMin))
delta_theta = 0.5625
print "Delta theta = %f degrees or %f miliradians"%(delta_theta, radians(delta_theta)*1000)
print "Corresponds to delta eta of %f at eta = 0 (fixed delta eta mean smaller delta theta as we go toward the beam)"%(get_eta(90-delta_theta) - get_eta(90))

theta_seg = arange(0, covering_angle + 2*delta_theta, delta_theta)
print theta_seg

theta_delta = []
cell_width_inner = []
cell_width_outer = []
total_width_sofar_inner = 0
total_width_sofar_outer = 0
total_angle_sofar = 0
for i in range(0, theta_seg.size - 1):
    theta_delta.append(delta_theta)
    total_angle_sofar += delta_theta
    li = tan(radians(total_angle_sofar)) * calo_rMin - total_width_sofar_inner
    total_width_sofar_inner += li
    lo = tan(radians(total_angle_sofar)) * calo_rMax - total_width_sofar_outer
    total_width_sofar_outer += lo

    cell_width_inner.append(li)
    cell_width_outer.append(lo)

#print "Delta Theta: ", [round(a, 3) for a in theta_delta]
print "Delta Theta: ", [a for a in theta_delta]
print "\n"
print "Inner radius widths: ", [round(a, 2) for a in cell_width_inner]
print "\n"
print "Outer radius widths: ", [round(a, 2) for a in cell_width_outer]
print "\n"
#Derive angle covering based on calo length
covering_angle = degrees(atan(calo_half_dz/calo_rMin))
print "For the calo dimension given, covering angle is up to eta ", get_eta(covering_angle)
print "Covering angle from both sides of radial direction: ", covering_angle
print "We cover from %f to %f theta degrees"%(90-covering_angle, 90+covering_angle)
print "Number of eta cells: %d * 2"%len(theta_delta)
print "Calo total width: %f"%(total_width_sofar_inner)
print "--------------------------------------"

# check when we have to cut the lines 
n_longitudinal_layer = 12
delta_theta = 0.5625/4
theta_seg = arange(0, covering_angle + 2*delta_theta, delta_theta)
cell_width_inner = []
cell_width_outer = []
total_width_sofar_outer = 0
total_angle_sofar = 0
found_breaking_point = False
for i in range(0, theta_seg.size - 1):
    total_angle_sofar += delta_theta
    li = tan(radians(total_angle_sofar)) * calo_rMin - total_width_sofar_inner
    lo = tan(radians(total_angle_sofar)) * calo_rMax - total_width_sofar_outer
    total_width_sofar_outer += lo
    if (total_width_sofar_outer > total_width_sofar_inner) and not found_breaking_point:
        found_breaking_point = True
        print "For the strip layer, we have to break lines starting from cell %d \n"%i
        print "For the regular layers, we have to break lines starting from cell %f \n"%(i/4.0)
        #print "In terms of absolute number of line (1 per radial segment per eta segment): %f \n"%(i*n_longitudinal_layer/4.0)
print "Total number fo diagonal lines for strip layer: %d"%(len(theta_delta)*4)
print "Total number fo diagonal lines for regular cells: %d"%(n_longitudinal_layer*len(theta_delta))



