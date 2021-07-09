#from sympy.parsing.mathematica import mathematica
#from sympy import var
from cmath import exp
import cmath
import os
#from numpy import complex

#s, tdr, i0, = var('s tdr i0')
#current_ps_mat = "(E-s tdr (1-Es tdr+Es tdr s tdr) i0)/(s2 tdr)"
#print(mathematica(current_ps_mat))

# cross-talk
frequencies = []
s_parameter = []
counter = 0
freq_str = "{"
s_parameter_str = "{"
with open("cross_talk_table.txt") as table_data:
    for line in table_data.readlines():
        if counter == 0:
            counter += 1
            continue
        entries = line.split()
        frequencies.append(float(entries[0]))
        s_parameter.append(float(entries[13]))
        freq_str += "%s, "%(entries[0].replace('E', '*^'))
        #print entries[0], entries[13]
        counter += 1
freq_str += "}"




def get_signal_current_laplace(s, i0 = 0.021, tdr = 295 * 10**-9):
    return i0 * (-1 + exp(-s * tdr) + s * tdr) / (s**2 * tdr)

def get_signal_current_ps(w, i0 = 0.021, tdr = 295 * 10**-9):
    return abs(i0 * (-1 + exp(-complex(0, w) * tdr) + complex(0, w) * tdr) / (complex(0, w)**2 * tdr))**2

def get_signal_current_fs(w, i0 = 0.021, tdr = 295 * 10**-9):
    return abs(i0 * (-1 + exp(-complex(0, w) * tdr) + complex(0, w) * tdr) / (complex(0, w)**2 * tdr))
import numpy as np

plot_dir = "plots_signalSpectrum"
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

s_range = range(10**6, 10**9, 1000)
signal_current_laplace = []
signal_current_ps = []

signal_current_laplace_weighted = []
signal_current_ps_weighted = []
signal_current_fs_weighted = []

for s in s_range:
    absolute_difference_function = lambda list_value : abs(list_value - s)
    closest_value = min(frequencies, key=absolute_difference_function)
    sparam = s_parameter[frequencies.index(closest_value)]
    signal_current_laplace.append(get_signal_current_laplace(s))
    signal_current_ps.append(get_signal_current_ps(s))

    signal_current_laplace_weighted.append(get_signal_current_laplace(s) * sparam)
    signal_current_ps_weighted.append(get_signal_current_ps(s) * sparam)
    signal_current_fs_weighted.append(get_signal_current_fs(s) * sparam)

import matplotlib.pyplot as plt
import numpy as np

s_parameter_array = np.array(s_parameter)
frequencies_array = np.array(frequencies)
signal_current_laplace_weighted_array = np.array(signal_current_laplace_weighted)
signal_current_ps_weighted_array = np.array(signal_current_ps_weighted)
signal_current_fs_weighted_array = np.array(signal_current_fs_weighted)

plt.figure("S parameters")
plt.plot(frequencies_array, s_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('w')
plt.ylabel('S parameter')
plt.grid()
plt.title('S parameters')
plt.gcf().savefig(os.path.join(plot_dir, "s_parameters.png"))
plt.clf()

plt.figure("signal_current_laplace_weighted")
plt.plot(np.array(s_range), signal_current_laplace_weighted_array)
axes = plt.gca()
axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('w')
plt.ylabel('$\mathcal{L}(i(t))_{s=w} * S_{21}$')
plt.grid()
plt.title('Cross-talk current laplace spectrum')
plt.gcf().savefig(os.path.join(plot_dir, "signal_current_laplace_weighted.png"))
plt.clf()

plt.figure("signal_current_fs_weighted")
plt.plot(np.array(s_range), signal_current_fs_weighted_array)
axes = plt.gca()
axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('w')
plt.ylabel('$|\mathcal{L}(i(t))_{s=iw}| * S_{21}$')
plt.grid()
plt.title('Cross-talk current frequency spectrum')
plt.gcf().savefig(os.path.join(plot_dir, "signal_current_fs_weighted.png"))
plt.clf()

plt.figure("signal_current_ps_weighted")
plt.plot(np.array(s_range), signal_current_ps_weighted_array)
axes = plt.gca()
axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('w')
plt.ylabel('$|\mathcal{L}(i(t))_{s=iw}|^{2} * S_{21}$')
plt.grid()
plt.title('Cross-talk current power spectrum')
plt.gcf().savefig(os.path.join(plot_dir, "signal_current_ps_weighted.png"))
plt.clf()

signal_current_laplace_array = np.array(signal_current_laplace)

signal_current_laplace_array = np.array(signal_current_laplace)

plt.figure("Signal current laplace")
plt.plot(np.array(s_range), signal_current_laplace_array)
axes = plt.gca()
axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('w')
plt.ylabel('$\mathcal{L}(i(t))_{s=w}$')
plt.grid()
plt.title('Signal current laplace spectrum')
plt.gcf().savefig(os.path.join(plot_dir, "original_signal_laplaceSpectrum.png"))
plt.clf()

signal_current_ps_array = np.array(signal_current_ps)

plt.figure("Signal current power density")
plt.plot(np.array(s_range), signal_current_ps_array)
axes = plt.gca()
axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('w')
plt.ylabel('$|\mathcal{L}(i(t))_{s=iw}|^{2}$')
plt.grid()
plt.title('Signal current power density')
plt.gcf().savefig(os.path.join(plot_dir, "original_signal_spectrum.png"))
plt.clf()


#print(frequencies)
#print(s_parameter)
#print(freq_str)
#    #string = table_data.read()
#    #print string
