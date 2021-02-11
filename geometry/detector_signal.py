import os
import numpy as np
from numpy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt
from matplotlib import rc

lar_gap = 1.24 # mm
drift_velocity = 4.75 # mm/microsecond = micron/ns (from fig 6 of http://cds.cern.ch/record/683899, 10 kV/cm and 87 Kelvin)
drift_time = 1000 * lar_gap / drift_velocity # ns
rise_time = 0.2 # ns  --> frequencies will be in GHz

step = 0.2 # ns granularity with which defining the functions

plot_dir = "signal_plots_risetime0dot2nsstep0dot2ns"
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

def pol1(x, a, b): # y = a * x + b
    return a * x + b

def pol1_from_two_points(x, x1, y1, x2, y2):
    a = (y2 -y1) / float(x2 - x1)
    b = y1 - (a * x1)
    return pol1(x, a, b)

# Allow using latex in plots
rc('text', usetex=True)

# signal definition
# signal height
average_energyDeposit_per_interaction = 23.3 # eV
average_energyDeposit_per_sensitiveCell = 5 # MeV derived from simulation without any sampling fraction applied
average_number_of_interaction = average_energyDeposit_per_sensitiveCell * 10**6 / average_energyDeposit_per_interaction # each interaction create one electron and one ion, but only the electrons are considered for signal height
print "Avergae number of interaction in one cell: %f"%average_number_of_interaction
electron_charge = 1.6 * 10**-19 # Amper * seconds (Coulomb)

# weighting potential from Riegler lecture 2 slide 4 and slide 7: for parralel plate it is equal to the potential V_w and does not depend on the position inside the sensitive gap where charge is deposited
# this lecture actually gives the induce current directly
i_induced_on_electrodes = (average_number_of_interaction * electron_charge * 10**6 * drift_velocity / lar_gap) * 10**6 # keep drift velocity in mm/microsec, take the electron charge in Amper * microsec --> the microsec cancel, keep LAr gap in mm --> mm from drift velocity cancel --> we have something in Amps --> multiply by 10**6 to get in microAmps
print "i_induced_on_electrodes: ", i_induced_on_electrodes, " microAmps"
signal_height = i_induced_on_electrodes # micro amps

#f_i_e = electron_charge / (drift_velocity * drift_time)  # nA/MeV
#FIXME E is not the electric field! it is the weighting field
#i_induced_on_electrodes = average_number_of_interaction * electron_charge * 10**6 * drift_velocity * electric_field # https://en.wikipedia.org/wiki/Shockley%E2%80%93Ramo_theorem   i = Eqv for one particle, put electric field at 1, set drift velocity in mm/microsec, take the electron charge in Amper * microsec --> the microsec cancel --> we get something in Amp --> multiply by 10**6 to get in microAmps
#i_induced_on_electrodes = average_number_of_interaction * electron_charge * 10**6 * drift_velocity * electric_field * 10**3  # https://en.wikipedia.org/wiki/Shockley%E2%80%93Ramo_theorem   i = Eqv for one particle, put electric field in volt/mm, keep drift velocity in mm/microsec --> the mm cancel and take the electron charge in Amper * microsec --> the microsec cancel --> we get something in Amp --> multiply by 10**6 to get in microAmps



# rising edge
p1x_rising = 0
p1y_rising = 0
p2x_rising = rise_time
p2y_rising = signal_height
x_rising = np.arange(0, rise_time, step)
y_rising = np.array([pol1_from_two_points(x, p1x_rising, p1y_rising, p2x_rising, p2y_rising) for x in x_rising])

# falling edge
p1x_falling = rise_time
p1y_falling = signal_height
p2x_falling = rise_time + drift_time 
p2y_falling = 0
x_falling = np.arange(rise_time, rise_time + drift_time, step)
y_falling = np.array([pol1_from_two_points(x, p1x_falling, p1y_falling, p2x_falling, p2y_falling) for x in x_falling])

# combine
x_signal = np.concatenate([x_rising, x_falling])
y_signal = np.concatenate([y_rising, y_falling])

# number of sampled points
n_sample = len(x_signal)
fftfreq = fftfreq(n_sample, step)[:n_sample//2]

print "Number of sampled points %d"%n_sample

sig_ft = fft(y_signal) # sig_fft gives an array with same length than y_signal (called N) whose entry k is the coefficient c_k such that f(t) = Sum c_k exp(iwkt). Entry 0 is the 'integral' of y_signal, entry 1 to N/2 - 1 corresponds to the positive frequen y terms, entry N/2 to N-1 are the negative frequency terms.
sig_from_fft = ifft(sig_ft)
#xf = fftfreq(N, T)[:N//2]

plt.figure("Signal shape")
plt.plot(x_signal, y_signal)
axes = plt.gca()
axes.set_xlim([0-step, rise_time + drift_time + step])
axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('time (ns)')
plt.ylabel(r'signal height ($\mu$A)')
plt.grid()
plt.title('Detector signal shape')
plt.gcf().savefig(os.path.join(plot_dir, "detector_signal.png"))

plt.figure("Signal FFT")
#xft = np.linspace(0.0, 1.0/(2.0*step), int((rise_time + drift_time)/float(step*2)))
plt.title('Fourier Transform')
#plt.plot(xft, np.abs(sig_ft))
#print  np.abs(sig_ft)
plt.plot(fftfreq, 2.0/float(n_sample) * np.abs(sig_ft[0:n_sample//2]))
plt.xlabel('Frequency (GHz)')
plt.gcf().savefig(os.path.join(plot_dir, "detector_signal_fft.png"))

#from scipy.signal import blackman
#w = blackman(n_sample)
#ywf = fft(y_signal*w)
#plt.figure("Signal FFT w. window")
#plt.plot(fftfreq, 2.0/float(n_sample) * np.abs(ywf[0:n_sample//2]))
#plt.gcf().savefig(os.path.join(plot_dir, "detector_signal_fft_withWindow.png"))

plt.figure("Signal FFT (logy)")
#plt.semilogy(np.abs(sig_ft))
#xft = np.linspace(0.0, 1.0/(2.0*step), int((rise_time + drift_time)/float(step*2)))
plt.title('Fourier Transform')
#plt.plot(xft, np.abs(sig_ft))
#print  np.abs(sig_ft)
plt.semilogy(fftfreq, 2.0/float(n_sample) * np.abs(sig_ft[0:n_sample//2]))
plt.xlabel('Frequency (GHz)')
plt.gcf().savefig(os.path.join(plot_dir, "detector_signal_fft_logy.png"))

#
#plt.figure("Signal FFT (logx)")
#plt.gcf().savefig(os.path.join(plot_dir, "detector_signal_fft_logx.png"))
#
#plt.figure("Signal FFT (log)")
#plt.figure()
#plt.yscale('log')
#plt.xscale('log')
#plt.plot(np.abs(sig_ft))
#plt.gcf().savefig(os.path.join(plot_dir, "detector_signal_fft_log.png"))

plt.figure("Signal from fft")
plt.figure()
plt.plot(np.abs(sig_from_fft))
plt.gcf().savefig(os.path.join(plot_dir, "detector_signal_fromfft.png"))


