import os, sys, math
import matplotlib.pyplot as plt
import numpy as np

plot_dir = "plots_noise"
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

noise_value_per_cell = 10 #MeV   0.1 cold pre amp and 10 pF,  0.5 for warm pre-amp and 10 pF cells
number_fired_cells = np.arange(1, 800) # average number of fired cell for: 1GeV elec all theta new geo no zero supp --> 78 (61 if above 1 MeV, 20 if above 10 MeV), 10GeV elec all theta new geo no zero supp --> 370
total_noise_list = []
for ncell in number_fired_cells:
    total_noise = math.sqrt(sum([noise_value_per_cell**2 for i in range(ncell)]))
    total_noise_list.append(total_noise)
total_noise_array = np.array(total_noise_list)

plt.figure("Cluster noise VS #fired cells")
plt.plot(number_fired_cells, total_noise_array)
axes = plt.gca()
#axes.set_xlim([0-step, rise_time + drift_time + step])
#axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
#axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
plt.xlabel('#fired cells')
plt.ylabel('Cluster noise [MeV]')
plt.grid()
plt.title('Cluster noise VS #fired cells')
plt.gcf().savefig(os.path.join(plot_dir, "noise_vs_fired_cells_noise_%f.png"%noise_value_per_cell))
