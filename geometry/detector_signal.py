import os
import numpy as np
from scipy import fft
import matplotlib.pyplot as plt
from matplotlib import rc

lar_gap = 1.24 # mm
drift_time = 275 # ns
rise_time = 1 # ns
signal_height = 1 # micro amps

step = 0.01 # granularity with which defining the functions

plot_dir = "signal_plots"
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
x_falling = np.arange(rise_time, drift_time, step)
y_falling = np.array([pol1_from_two_points(x, p1x_falling, p1y_falling, p2x_falling, p2y_falling) for x in x_falling])

# combine
x_signal = np.concatenate([x_rising, x_falling])
y_signal = np.concatenate([y_rising, y_falling])

plt.subplot(511)
plt.plot(x_signal, y_signal)
axes = plt.gca()
axes.set_xlim([0-step, rise_time + drift_time + step])
axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('time (ns)')
plt.ylabel(r'signal height ($\mu$A)')
plt.grid()
plt.title('Detector signal shape')

plt.subplot(512)
sig_ft = fft(y_signal)
xft = np.linspace(0.0, 1.0/(2.0*step), int((rise_time + drift_time)/float(step*2)))
plt.title('Fourrier Transform')
#plt.plot(xft, np.abs(sig_ft))
print  np.abs(sig_ft)
plt.plot(np.abs(sig_ft))

plt.subplot(513)
sig_ft = fft(y_signal)
plt.semilogy(np.abs(sig_ft))

plt.subplot(514)
sig_ft = fft(y_signal)
plt.semilogx(np.abs(sig_ft))

plt.subplot(515)
sig_ft = fft(y_signal)
plt.yscale('log')
plt.xscale('log')
plt.plot(np.abs(sig_ft))
plt.show()
#y = fft(x)
#print y
plt.savefig(os.path.join(plot_dir, "detector_signal.png"))
