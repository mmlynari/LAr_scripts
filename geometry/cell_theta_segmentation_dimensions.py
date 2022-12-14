from math import tan, log, atan, exp, degrees, radians
from numpy import arange

#everything in degrees
def get_theta(eta):
    theta = 2 * atan(exp(-eta))
    return degrees(theta)
def get_eta(theta):
    eta = -log(tan(radians(theta/2)))
    return eta

calo_rMin = 216
cryo_start = 210
calo_rMax = 256
#calo_half_dz = 220
dr_ps = 1.5
dr_other = 3.5
cryo_stat_thickness_side = 10 #cm
covering_angle_sensitive = 45 #degrees(atan(calo_half_dz/calo_rMin)), this is meant to not include crystat
delta_theta = 0.5625
calo_sensitive_half_dz = calo_rMin * tan(radians(covering_angle_sensitive)) # including cryostat
print("Calo sensitive half width for the covering angle %f and radius %f : %f"%(covering_angle_sensitive, calo_rMin, calo_sensitive_half_dz))
print("Sensitive covering angle in eta: %f"%(get_eta(covering_angle_sensitive)))
calo_half_dz = calo_sensitive_half_dz + cryo_stat_thickness_side
covering_angle = degrees(atan(calo_half_dz/calo_rMin))
print("Total calo half width including cryostat for the sensitive covering angle %f and radius %f : %f which has covering angle of %f corresponding to eta = %f"%(covering_angle_sensitive, calo_rMin, calo_half_dz, covering_angle, get_eta(covering_angle)))
print("We cover from %f to %f theta degrees including cryostat"%(90-covering_angle, 90+covering_angle))
print("Top put in the detector description is dimension including cryostat (it is removed later): %f"%calo_half_dz)
print("Top put in the detector description for the theta offset in radians: %f"%(-1*atan(calo_half_dz/float(cryo_start))))
print("Top put in the latex drawing: %f with delta_theta of %f"%(calo_sensitive_half_dz, delta_theta))
print("Delta theta = %f degrees or %f miliradians"%(delta_theta, radians(delta_theta)*1000))
print("Corresponds to delta eta of %f at eta = 0 (fixed delta eta mean smaller delta theta as we go toward the beam)"%(get_eta(90-delta_theta) - get_eta(90)))


theta_seg = arange(0, covering_angle_sensitive+delta_theta, delta_theta)
print(theta_seg)

theta_delta = []
cell_width_inner = []
cell_width_outer = []
total_width_sofar_inner = 0
total_width_sofar_outer = 0
total_angle_sofar = 0
for i in range(0, theta_seg.size+1):
    theta_delta.append(delta_theta)
    total_angle_sofar += delta_theta
    li = tan(radians(total_angle_sofar)) * calo_rMin - total_width_sofar_inner
    total_width_sofar_inner += li
    lo = tan(radians(total_angle_sofar)) * calo_rMax - total_width_sofar_outer
    total_width_sofar_outer += lo

    cell_width_inner.append(li)
    cell_width_outer.append(lo)

#print "Delta Theta: ", [round(a, 3) for a in theta_delta]
print("Delta Theta: ", [a for a in theta_delta])
print("\n")
print("Inner radius widths: ", [round(a, 2) for a in cell_width_inner])
print("\n")
print("Outer radius widths: ", [round(a, 2) for a in cell_width_outer])
print("\n")
#Derive angle covering based on calo length
#covering_angle = degrees(atan(calo_half_dz/calo_rMin))
print("Number of eta cells: %d * 2"%(len(theta_delta)-2))#list contain 0 and 45
print("--------------------------------------")

# check when we have to cut the lines 
n_longitudinal_layer = 12
cell_width_inner = []
cell_width_outer = []
total_width_sofar_outer = 0
total_angle_sofar = 0
found_breaking_point = False
# regular cells
for i in range(0, theta_seg.size):
    total_angle_sofar += delta_theta
    lo = tan(radians(total_angle_sofar)) * calo_rMax - total_width_sofar_outer
    total_width_sofar_outer += lo
    if (total_width_sofar_outer >= calo_sensitive_half_dz) and not found_breaking_point:
        found_breaking_point = True
        print("For the regular layers, we have to break lines starting from cell %f"%(round(i)))
print("Total number of diagonal lines for regular cells: %d"%(len(theta_seg)-1))
# strip layer
delta_theta = delta_theta/4.0
theta_seg = arange(0, covering_angle_sensitive + delta_theta, delta_theta)
#print "Strip theta seg: "
#print theta_seg
cell_width_inner = []
cell_width_outer = []
total_width_sofar_outer = 0
total_angle_sofar = 0
found_breaking_point = False
for i in range(0, theta_seg.size * 2):
    total_angle_sofar += delta_theta
    lo = tan(radians(total_angle_sofar)) * (calo_rMin + dr_ps) - total_width_sofar_outer
    total_width_sofar_outer += lo
    if (total_width_sofar_outer >= calo_sensitive_half_dz) and not found_breaking_point:
        found_breaking_point = True
        print("For the strip layer, we have to break lines starting from cell %d"%(i-3))
        break
print("Total number of diagonal lines for strip layer: %d"%(len(theta_seg)-3)) # because you have 0 and 45 in this list



## centimeters
#calo_half_dz = 220
#calo_rMin = 216
#calo_rMax = 256
#dr_ps = 1.5
#dr_other = 3.5

#eta_seg = arange(0, 0.89, 0.01)
##print "Eta: ", eta_seg
##print "\n"
#theta_seg = [get_theta(eta) for eta in eta_seg]
##print "Theta: ", [round(a, 2) for a in theta_seg]
##print "\n"
#theta_delta = []
#cell_width_inner = []
#cell_width_outer = []
#total_width_sofar_inner = 0
#total_width_sofar_outer = 0
#total_angle_sofar = 0
#for i in range(0, eta_seg.size - 1):
#    delta_theta = theta_seg[i] - theta_seg[i+1]
#    theta_delta.append(delta_theta)
#    total_angle_sofar += delta_theta
#    li = tan(radians(total_angle_sofar)) * calo_rMin - total_width_sofar_inner
#    total_width_sofar_inner += li
#    lo = tan(radians(total_angle_sofar)) * calo_rMax - total_width_sofar_outer
#    total_width_sofar_outer += lo
#
#    cell_width_inner.append(li)
#    cell_width_outer.append(lo)

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

