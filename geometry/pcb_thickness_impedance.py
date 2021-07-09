from math import sqrt, log # log is neper logarithm

# everything in microns
h_hv = 65 # default in CDR is 100
h_m = 242.5 # default in CDR is 285
h_s = 170
trace_thickness = 35
other_conductor_thickness = 35
w_t = 127
dielectric_constant = 4.0
total_thickness_no_trace_no_conductor = (h_hv + h_m + h_s) * 2
total_thickness_no_trace = total_thickness_no_trace_no_conductor + 2 * other_conductor_thickness * 3 # there are three conductor layers per half PCB (signal trace will be added after since only one is present)
total_thickness = total_thickness_no_trace + trace_thickness
print("Total PCB thickness: %f mm\n"%(total_thickness/1000.0))

# check that all dimensions comply
# avoid spark between HV layer and signal pads (H_hv)
electrical_rigidity_FR4 = 20 #kV/mm    electric field that it can sustain without creating shortcircuit
wished_electric_field_in_gap = 1 # kV/mm
LAr_gap_thickness = 1.24 # mm
potential_hv_needed = wished_electric_field_in_gap * LAr_gap_thickness # kV
minimum_space_hv_sigpad = potential_hv_needed / electrical_rigidity_FR4
print("For a LAr gap of %f mm and a voltage on HV plate of %f kV:"%(LAr_gap_thickness, potential_hv_needed))
print("\t Mimimum spacing between signal pads and HV plate to avoid electrical breackdown: %d microns"%(minimum_space_hv_sigpad*1000))
print("\t Foreseen space: %d microns\n"%(h_hv))

# check the transmission line impedence
def get_impedence(height, trace_width, trace_thickness, dielectric_constant):
    return 60 * log(1.9 * (2 * height + trace_thickness) / (0.8 * trace_width + trace_thickness))  / sqrt(dielectric_constant)
space_between_inner_side_ground = h_s * 2 + trace_thickness 
impedence = get_impedence(h_s, w_t, trace_thickness, dielectric_constant)
print("Transmission line impedence for h_s=%d, w_t=%d, trace_thickness=%d and dielectric constant=%f: %f Ohm"%(h_s, w_t, trace_thickness, dielectric_constant, impedence))

import os, sys, math
import matplotlib.pyplot as plt
import numpy as np

plot_dir = "plots_impedance"
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

impedance_vs_thickness = []
impedance_vs_width = []
impedance_vs_height = []
impedance_vs_diel = []

thicknesses = range(10, 70, 1)
print(thicknesses)
for thickness in thicknesses:
    impedance_vs_thickness.append(get_impedence(h_s, w_t, thickness, dielectric_constant))
impedance_vs_thickness_array = np.array(impedance_vs_thickness)

plt.figure("Stripline impedance")
plt.plot(np.array(thicknesses), impedance_vs_thickness_array)
axes = plt.gca()
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('Trace thickness [$\mu$m]')
plt.ylabel('Stripline impedance [$\Omega$]')
plt.grid()
plt.title('Impedance vs trace thickness')
plt.gcf().savefig(os.path.join(plot_dir, "impedance_vs_thickness.png"))
plt.clf()

widthes = range(50, 250, 1)
print(widthes)
for width in widthes:
    impedance_vs_width.append(get_impedence(h_s, width, trace_thickness, dielectric_constant))
impedance_vs_width_array = np.array(impedance_vs_width)

plt.figure("Stripline impedance")
plt.plot(np.array(widthes), impedance_vs_width_array)
axes = plt.gca()
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('Trace width [$\mu$m]')
plt.ylabel('Stripline impedance [$\Omega$]')
plt.grid()
plt.title('Impedance vs trace width')
plt.gcf().savefig(os.path.join(plot_dir, "impedance_vs_width.png"))
plt.clf()

heightes = range(75, 455, 1)
print(heightes)
for height in heightes:
    impedance_vs_height.append(get_impedence(height, w_t, trace_thickness, dielectric_constant))
impedance_vs_height_array = np.array(impedance_vs_height)

plt.figure("Stripline impedance")
plt.plot(np.array(heightes), impedance_vs_height_array)
axes = plt.gca()
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('H [$\mu$m]')
plt.ylabel('Stripline impedance [$\Omega$]')
plt.grid()
plt.title('Impedance vs H')
plt.gcf().savefig(os.path.join(plot_dir, "impedance_vs_height.png"))
plt.clf()


dieles = np.arange(3.0, 5.0, 0.01)
print(dieles)
for diel in dieles:
    impedance_vs_diel.append(get_impedence(h_s, w_t, trace_thickness, diel))
impedance_vs_diel_array = np.array(impedance_vs_diel)

plt.figure("Stripline impedance")
plt.plot(np.array(dieles), impedance_vs_diel_array)
axes = plt.gca()
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('Dielectric constant')
plt.ylabel('Stripline impedance [$\Omega$]')
plt.grid()
plt.title('Impedance vs dielectric constant')
plt.gcf().savefig(os.path.join(plot_dir, "impedance_vs_diel.png"))
plt.clf()


# do the same for microstrip


