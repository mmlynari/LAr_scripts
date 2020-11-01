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
print "Eta: ", eta_seg
print "\n"
theta_seg = [get_theta(eta) for eta in eta_seg]
print "Theta: ", [round(a, 2) for a in theta_seg]
print "\n"
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

print "Delta Theta: ", [round(a, 3) for a in theta_delta]
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
print "PCB length: %f" % (calo_rMax - calo_rMin)
