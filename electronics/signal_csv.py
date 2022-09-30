# write the input current with time_stamp current_value
import matplotlib.pyplot as plt
import numpy as np


# cross-talk
times = np.arange(0, 500, 0.1)

def get_signal_current(t, i0 = 0.021, tdr = 300, tstart = 10, trise = 1):
    if t > tdr + tstart or t < tstart:
        return 0
    elif tstart <= t <= tstart + trise:
        return i0*(t-tstart)/trise
    else:
        return i0*(1-(t-(tstart+trise))/tdr)

current = []
with open("input_signal.txt", 'w') as current_data:
    current_data.write("# ns uA\n")
    for time in times:
        current_data.write("%f %f\n"%(time, get_signal_current(time)))
        current.append(get_signal_current(time))
    current_data.close()


plt.figure("Input current")
plt.plot(np.array(times), np.array(current))
axes = plt.gca()
#axes.set_xscale('log')
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('t')
plt.ylabel('Input current [nA]')
plt.grid()
plt.title('Signal input current')
plt.gcf().savefig("input_current_forSimu.png")
plt.clf()
