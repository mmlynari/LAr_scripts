#from sympy.parsing.mathematica import mathematica
#from sympy import var
from cmath import exp
import cmath
import os
import numpy as np
from numpy import degrees, angle, radians
import matplotlib.pyplot as plt
#from numpy import complex

inputFile_path = "/eos/user/b/brfranco/EMDesignFiles/cross_talk_results/baseline_portWithGND_noAbsorber_Ymatrrix.tab"

#s, tdr, i0, = var('s tdr i0')
#current_ps_mat = "(E-s tdr (1-Es tdr+Es tdr s tdr) i0)/(s2 tdr)"
#print(mathematica(current_ps_mat))
#! Terminal data exported
#! Port[1] = cell6_I_TP19
#! Port[2] = cell6_O_TP6
#! Port[3] = cell8_I_TP21
#! Port[4] = cell8_O_TP8
#! Port[5] = cell15_I_TP30
#! Port[6] = cell15_O_TP28

# cross-talk
frequencies = []
#y_parameters = {"cell6_I-cell6_O" : [], "cell6_I-cell8_O" : []}
y_parameters = []
#treating current induced on cell15_O_TP28 from cell6_I_TP19
y11s_entry = 1
y11s = []
y22s_entry = 71
y22s = []
y12s_entry = 11 # or 61 for y21
y12s = []
y21s_entry = 61
y21s = []
counter = 0
freq_str = "{"
y_parameter_str = "{"
with open(inputFile_path) as table_data:
    for line in table_data.readlines():
        if counter == 0:
            counter += 1
            continue
        entries = line.split()
        #print(float(entries[0]))
        frequencies.append(float(entries[0]))
        #print(entries[y11s_entry], " ", entries[y11s_entry+1])
        #y11s.append(complex(float(entries[y11s_entry]), float(entries[y11s_entry + 1])))

        y11 = float(entries[y11s_entry]) * exp(1j * radians(float(entries[y11s_entry + 1])))
        y22 = float(entries[y22s_entry]) * exp(1j * radians(float(entries[y22s_entry + 1])))
        y12 = float(entries[y12s_entry]) * exp(1j * radians(float(entries[y12s_entry + 1])))
        y21 = float(entries[y21s_entry]) * exp(1j * radians(float(entries[y21s_entry + 1])))

        y11s.append(y11)
        y22s.append(y22)
        y12s.append(y12)
        y21s.append(y21)

        #print(abs(y11), " ", degrees(angle(y11)))
        #y22s.append(complex(float(entries[y22s_entry]), float(entries[y22s_entry + 1])))
        #y12s.append(complex(float(entries[y12s_entry]), float(entries[y12s_entry + 1])))
        #print(complex(float(entries[y11s_entry]), float(entries[y11s_entry + 1])))
        #print(entries[y12s_entry], " ", entries[61])
        #print(entries[y22s_entry])
        #y_parameters.append(float(entries[13]))
        freq_str += "%s, "%(entries[0].replace('E', '*^'))
        #print entries[0], entries[13]
        counter += 1
        #if counter == 10:
        #    break
freq_str += "}"

def get_output_current_fs(input_current, y11, y12, y22, z = 50):
    return input_current * 1/(y11 * (1 - y22 * z) / y12 + y12 * z)
    #return input_current * 1/(y11 * (1 - y22 * z) / y21 + y12 * z)

def get_signal_current_laplace(s, i0 = 21, tdr = 300 * 10**-9):
    return i0 * (-1 + exp(-s * tdr) + s * tdr) / (s**2 * tdr)

def get_signal_current_ps(w, i0 = 21, tdr = 300 * 10**-9):
    return abs(i0 * (-1 + exp(-complex(0, w) * tdr) + complex(0, w) * tdr) / (complex(0, w)**2 * tdr))**2

def get_signal_current_fs(w, i0 = 21, tdr = 300 * 10**-9):
    return abs(i0 * (-1 + exp(-complex(0, w) * tdr) + complex(0, w) * tdr) / (complex(0, w)**2 * tdr))

def get_signal_current_fs_complex(w, i0 = 21, tdr = 300 * 10**-9):
    return i0 * (-1 + exp(-complex(0, w) * tdr) + complex(0, w) * tdr) / (complex(0, w)**2 * tdr)


def get_signal_current(t, i0 = 21, tdr = 300, tstart = 10):
    if t > tdr + tstart or t < tstart:
        return 0
    else:
        return i0*(1-(t-tstart)/tdr)


step = 0.5
times = np.arange(0, 500, step)
#print(times)
input_current_time_domain = []
for time in times:
    input_current_time_domain.append(get_signal_current(time))


plot_dir = "plots_signalSpectrum_yparams"
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

# Input signal time domain
plt.figure("Input current in the time domain")
plt.plot(np.array(times), np.array(input_current_time_domain))
axes = plt.gca()
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('t [ns]')
plt.ylabel('i [nA]')
plt.grid()
plt.title('Input current in the time domain')
plt.gcf().savefig(os.path.join(plot_dir, "input_current_time_domain.png"))
plt.clf()


from numpy.fft import fft, ifft, fftfreq

#input signal frequency domain from python FFT
n_sample = len(input_current_time_domain)
fftfreq = fftfreq(n_sample, step)[:n_sample//2]
#fftfreq = fftfreq(n_sample, step)[:n_sample//2]
input_current_frequency_domain = fft(input_current_time_domain)

plt.figure("Input signal Frequency domain from FFT")
#print(fftfreq)
#print(input_current_frequency_domain)
#plt.plot(fftfreq, 2.0/float(n_sample) * np.abs(input_current_frequency_domain[0:n_sample//2]))
plt.plot(fftfreq, np.abs(input_current_frequency_domain[0:n_sample//2]))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|FFT|')
plt.grid()
plt.title('Input signal Frequency domain from FFT')
plt.gcf().savefig(os.path.join(plot_dir, "input_sig_freq_domain_fromFFT.png"))
plt.clf()

sig_from_fft = ifft(input_current_frequency_domain)
plt.figure("Input current in the time domain from FFT")
plt.plot(np.array(times), np.abs(sig_from_fft))
axes = plt.gca()
plt.xlabel('t [ns]')
plt.ylabel('i [nA]')
plt.grid()
plt.title('Input current in the time domain from FFT')
plt.gcf().savefig(os.path.join(plot_dir, "input_current_time_domain_fromFFT.png"))
plt.clf()

# Y parameters plots
plt.figure("Y11 mag")
plt.plot(frequencies, np.abs(y11s))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|Y_{11}|')
plt.grid()
plt.title('|Y_{11}|')
plt.gcf().savefig(os.path.join(plot_dir, "y11_mag.png"))
plt.clf()

plt.figure("Y11 phase")
plt.plot(frequencies, np.angle(y11s, deg=True))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|Y_{11}|')
plt.grid()
plt.title('|Y_{11}|')
plt.gcf().savefig(os.path.join(plot_dir, "y11_phase.png"))
plt.clf()

plt.figure("Y12 mag")
plt.plot(frequencies, np.abs(y12s))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|Y_{12}|')
plt.grid()
plt.title('|Y_{12}|')
plt.gcf().savefig(os.path.join(plot_dir, "y12_mag.png"))
plt.clf()

plt.figure("Y12 phase")
plt.plot(frequencies, np.angle(y12s, deg=True))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|Y_{12}|')
plt.grid()
plt.title('|Y_{12}|')
plt.gcf().savefig(os.path.join(plot_dir, "y12_phase.png"))
plt.clf()

plt.figure("Y21 mag")
plt.plot(frequencies, np.abs(y21s))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|Y_{21}|')
plt.grid()
plt.title('|Y_{21}|')
plt.gcf().savefig(os.path.join(plot_dir, "y21_mag.png"))
plt.clf()

plt.figure("Y21 phase")
plt.plot(frequencies, np.angle(y21s, deg=True))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|Y_{21}|')
plt.grid()
plt.title('|Y_{21}|')
plt.gcf().savefig(os.path.join(plot_dir, "y21_phase.png"))
plt.clf()

plt.figure("Y22 mag")
plt.plot(frequencies, np.abs(y22s))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|Y_{22}|')
plt.grid()
plt.title('|Y_{22}|')
plt.gcf().savefig(os.path.join(plot_dir, "y22_mag.png"))
plt.clf()

plt.figure("Y22 phase")
plt.plot(frequencies, np.angle(y22s, deg=True))
#plt.plot(frequencies_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('w')
plt.ylabel('|Y_{22}|')
plt.grid()
plt.title('|Y_{22}|')
plt.gcf().savefig(os.path.join(plot_dir, "y22_phase.png"))
plt.clf()


#frequency_range = range(10**6, 10**9, 1000)
#frequency_range = range(0, 10**9, 500)
frequency_range = range(0, 10**9, 500)
signal_current_laplace = []
signal_current_ps = []
signal_current_fs_fromMathematica = []

signal_current_laplace_weighted = []
signal_current_ps_weighted = []
signal_current_fs_weighted = []

count = 0
crosstalk_current_freq_dom = []
frequency_range = np.arange(0, 10**9, 500)
#frequency_range = np.linspace(0, 10**9, 1000)
#frequency_range = range(0, 10**9, 500)
print(len(frequency_range))
#for w in frequencies:
for w in frequency_range:
    if w == 0:#needed to compute the limit when w --> 0 from Mathematica, otherwise it is 0/0 
        print("Treated 0 frequency")
        phase_0 = 0
        abs_0 = 3.0975*10**-9  
        signal_current_fs_fromMathematica.append(abs_0)
        crosstalk_current_freq_dom.append(get_output_current_fs(3.0975*10**-9 + 1j*0, y11s[count], y12s[count], y22s[count]))
    else:
        absolute_difference_function = lambda list_value : abs(list_value - w)
        closest_value = min(frequencies, key=absolute_difference_function)
        signal_current_fs_fromMathematica.append(get_signal_current_fs(w))
        #crosstalk_current_freq_dom.append(get_output_current_fs(get_signal_current_fs_complex(w), y11s[count], y12s[count], y22s[count]))
        crosstalk_current_freq_dom.append(get_output_current_fs(get_signal_current_fs_complex(w), y11s[frequencies.index(closest_value)], y12s[frequencies.index(closest_value)], y22s[frequencies.index(closest_value)]))
    count += 1

plt.figure("signal_current_fs_fromMathematica")
plt.plot(np.array(frequency_range), np.array(signal_current_fs_fromMathematica))
axes = plt.gca()
axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('w')
plt.ylabel('$|\mathcal{L}(i(t))_{s=iw}|$')
plt.grid()
plt.title('Signal current frequency spectrum from mathematica')
plt.gcf().savefig(os.path.join(plot_dir, "signal_current_fs_fromMathematica.png"))
plt.clf()

plt.figure("Crosstalk signal frequency domain")
plt.plot(np.array(frequency_range), np.abs(crosstalk_current_freq_dom))
#plt.plot(frequency_range_array, y_parameter_array)
axes = plt.gca()
axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('w')
plt.ylabel('|FFT|')
plt.grid()
plt.title('Crosstalk signal frequency domain')
plt.gcf().savefig(os.path.join(plot_dir, "crosstalk_sig_freq_domain.png"))
plt.clf()

print("Start to compute the inverse FFT")
#Cross-talk signal time domain
crosstalk_current_freq_dom_forFFT = crosstalk_current_freq_dom + crosstalk_current_freq_dom[::-1][:-1] # need to append the negative frequency (without the DC 0 Hz point) starting from the most negative frequence (::-1 reverse the list order) 
#crosstalk_current_freq_dom_forFFT = crosstalk_current_freq_dom + crosstalk_current_freq_dom[1:]
crosstalk_current_time_domain_fromiFFT = ifft(crosstalk_current_freq_dom_forFFT)

plt.figure("Cross-talk current in the time domain")
plt.plot(np.abs(crosstalk_current_time_domain_fromiFFT[:len(crosstalk_current_time_domain_fromiFFT)//2]))
axes = plt.gca()
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('t [ns]')
plt.ylabel('i [nA]')
plt.grid()
plt.title('Cross talk current in the time domain')
plt.gcf().savefig(os.path.join(plot_dir, "crosstalk_current_time_domain.png"))
plt.clf()

#print(len(crosstalk_current_freq_dom_forFFT), " ", len(crosstalk_current_freq_dom), " ", len(crosstalk_current_time_domain_fromiFFT), " ", len(times))
print(crosstalk_current_time_domain_fromiFFT)


#for s in frequency_range:
#    absolute_difference_function = lambda list_value : abs(list_value - s)
#    closest_value = min(frequencies, key=absolute_difference_function)
#    y11 = y11s[frequencies.index(closest_value)]
#    y22 = y22s[frequencies.index(closest_value)]
#    y21 = y21s[frequencies.index(closest_value)]
#    y12 = y12s[frequencies.index(closest_value)]

#    sparam = y_parameters[frequencies.index(closest_value)]
#    signal_current_laplace.append(get_signal_current_laplace(s))
#    signal_current_ps.append(get_signal_current_ps(s))
#
#    signal_current_laplace_weighted.append(get_signal_current_laplace(s) * sparam)
#    signal_current_ps_weighted.append(get_signal_current_ps(s) * sparam)
#    signal_current_fs_weighted.append(get_signal_current_fs(s) * sparam)
#
#
#y_parameter_array = np.array(y_parameters)
#frequencies_array = np.array(frequencies)
#signal_current_laplace_weighted_array = np.array(signal_current_laplace_weighted)
#signal_current_ps_weighted_array = np.array(signal_current_ps_weighted)
#signal_current_fs_weighted_array = np.array(signal_current_fs_weighted)
#
#plt.figure("S parameters")
#plt.plot(frequencies_array, y_parameter_array)
#axes = plt.gca()
#axes.set_xscale('log')
##axes.set_xlim([0-step, rise_time + drift_time + step])
##axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
##axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
#plt.xlabel('w')
#plt.ylabel('S parameter')
#plt.grid()
#plt.title('S parameters')
#plt.gcf().savefig(os.path.join(plot_dir, "y_parameters.png"))
#plt.clf()
#
#plt.figure("signal_current_laplace_weighted")
#plt.plot(np.array(frequency_range), signal_current_laplace_weighted_array)
#axes = plt.gca()
#axes.set_xscale('log')
##axes.set_xlim([0-step, rise_time + drift_time + step])
##axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
##axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
#plt.xlabel('w')
#plt.ylabel('$\mathcal{L}(i(t))_{s=w} * S_{21}$')
#plt.grid()
#plt.title('Cross-talk current laplace spectrum')
#plt.gcf().savefig(os.path.join(plot_dir, "signal_current_laplace_weighted.png"))
#plt.clf()
#
#plt.figure("signal_current_fs_weighted")
#plt.plot(np.array(frequency_range), signal_current_fs_weighted_array)
#axes = plt.gca()
#axes.set_xscale('log')
##axes.set_xlim([0-step, rise_time + drift_time + step])
##axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
##axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
#plt.xlabel('w')
#plt.ylabel('$|\mathcal{L}(i(t))_{s=iw}| * S_{21}$')
#plt.grid()
#plt.title('Cross-talk current frequency spectrum')
#plt.gcf().savefig(os.path.join(plot_dir, "signal_current_fs_weighted.png"))
#plt.clf()
#
#plt.figure("signal_current_ps_weighted")
#plt.plot(np.array(frequency_range), signal_current_ps_weighted_array)
#axes = plt.gca()
#axes.set_xscale('log')
##axes.set_xlim([0-step, rise_time + drift_time + step])
##axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
##axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
#plt.xlabel('w')
#plt.ylabel('$|\mathcal{L}(i(t))_{s=iw}|^{2} * S_{21}$')
#plt.grid()
#plt.title('Cross-talk current power spectrum')
#plt.gcf().savefig(os.path.join(plot_dir, "signal_current_ps_weighted.png"))
#plt.clf()
#
#signal_current_laplace_array = np.array(signal_current_laplace)
#
#signal_current_laplace_array = np.array(signal_current_laplace)
#
#plt.figure("Signal current laplace")
#plt.plot(np.array(frequency_range), signal_current_laplace_array)
#axes = plt.gca()
#axes.set_xscale('log')
##axes.set_xlim([0-step, rise_time + drift_time + step])
##axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
##axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
#plt.xlabel('w')
#plt.ylabel('$\mathcal{L}(i(t))_{s=w}$')
#plt.grid()
#plt.title('Signal current laplace spectrum')
#plt.gcf().savefig(os.path.join(plot_dir, "original_signal_laplaceSpectrum.png"))
#plt.clf()
#
#signal_current_ps_array = np.array(signal_current_ps)
#
#plt.figure("Signal current power density")
#plt.plot(np.array(frequency_range), signal_current_ps_array)
#axes = plt.gca()
#axes.set_xscale('log')
##axes.set_xlim([0-step, rise_time + drift_time + step])
##axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
##axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
#plt.xlabel('w')
#plt.ylabel('$|\mathcal{L}(i(t))_{s=iw}|^{2}$')
#plt.grid()
#plt.title('Signal current power density')
#plt.gcf().savefig(os.path.join(plot_dir, "original_signal_spectrum.png"))
#plt.clf()
#
#
##print(frequencies)
##print(y_parameter)
##print(freq_str)
##    #string = table_data.read()
##    #print string
