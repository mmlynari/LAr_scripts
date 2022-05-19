#!/usr/bin/env python
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
import os

# Make sure you have ticked the option "uniform step" in the data export from ANSYS

#f_sig = "signal_csvs/cell7_input_hv_0_shields_allCells.csv"
#f_sig = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/210927/LAr_scripts/electronics/signal_csvs/cell7_input_hv_0_shields_allCells_83Ohms_longXrange.csv"
#f_sig = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/210927/LAr_scripts/electronics/signal_csvs/cell7_input_hv_0_shields_allCells_83Ohms_longXrange_zeroPaddedUntil200ns.csv"
#f_sig = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/210927/LAr_scripts/electronics/signal_csvs/cell7_input_hv_0_shields_allCells_83Ohms_100psStepUniform_2usecrange.csv"
f_sig = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/210927/LAr_scripts/electronics/signal_csvs/cell7_input_hv_1_shields_allCells_83Ohm_2usecXrange.csv"
#plotfolder = "xtalk_plots_cell7input_0shields_longXrange_zeroPaddedTill200ns"
plotfolder = "xtalk_plots_cell7input_1shields_83Ohms_100psStepUniform"
if not os.path.isdir(plotfolder):
    os.mkdir(plotfolder)

x_axis_string_from_csv = "Time"
string_aggressor_cell = "Output"

df_sig = pd.read_csv(f_sig, sep=',', skiprows=0)
n_columns = df_sig.shape[1]
headers = df_sig.columns.tolist()
headers = headers[-1:] + headers[:-1]
df_sig = df_sig[headers]
for header in headers:
    if x_axis_string_from_csv in header:
        x_axis_string_from_csv = header
        continue
    df_sig = df_sig.rename(columns={header:header.replace(" [nA]", "").replace("I(", "").replace(")", "").replace("Inegative(", "").replace("_O", "").replace('cell', 'Cell ')})


# plot all the curves from the csv
my_timeLimit = 400
if 'us' in x_axis_string_from_csv:
    print("Warning: not sure everything is still ok for time axis in micro secs! Please check!")
    my_timeLimit = my_timeLimit / 1000.0

plt.clf()
df_sig.plot(x = x_axis_string_from_csv, ylabel = "Current [nA]", subplots = True, xlim = (0, my_timeLimit), layout = (int(n_columns/2.0) + 1, 2), figsize = (12, 12), grid = True)
plt.savefig(os.path.join(plotfolder, 'all_currents.png'))

#do not place this before to do the plot above 
if 'us' in x_axis_string_from_csv: # since tau is defined in nanosec, everything must be in nano sec
    df_sig[x_axis_string_from_csv] *= 1000
max_time_from_ansys = df_sig[x_axis_string_from_csv].max()

dict_shapingTime_cell_peakShaper = {}
sampling_frequency = 1 / df_sig[x_axis_string_from_csv].diff()[1]
print("Sampling frequency: ", sampling_frequency)
taus = [0, 20, 50, 70, 100, 150, 200, 300]
for tau in taus:
    print("Shaping time: ", tau)
    # define CR-RC^2 filter
    # order of the filter, critical (cut-off) frequency of the filter, filter type, analog or digital, type of output (pole zero, num/dem or second order section), sampling frequency
    if tau == 0:# means no shaper
        ylabel = "Current at PCB output port [nA]"
        extra_title = ", no shaper"
        CR = signal.butter(0, 1, btype='highpass', analog=False, output='sos', fs=sampling_frequency)
        RC = signal.butter(0, 1, btype='lowpass', analog=False, output='sos', fs=sampling_frequency)
        maxXrange = my_timeLimit
    else:
        ylabel = "Shaper output (arbitrary unit)"
        extra_title = ", shaping time " + str(tau) + " ns"
        CR = signal.butter(1, 1/(2*np.pi*tau), btype='highpass', analog=False, output='sos', fs=sampling_frequency)
        RC = signal.butter(2, 1/(2*np.pi*tau), btype='lowpass', analog=False, output='sos', fs=sampling_frequency)
        maxXrange = min(25 * tau, max_time_from_ansys)
    dict_shapingTime_cell_peakShaper[str(tau)] = {}
    cell_header = ""
    column_string = ""
    cells = []
    for label, currents in df_sig.items():
        if label == x_axis_string_from_csv or label == "Input":
            continue
        if label == string_aggressor_cell:
            plot_label = "Cell 7 (aggressor)" + extra_title
        else:
            cell_header += " " + label + " &"
            plot_label = label + " (victim)" + extra_title
        # shape signal 
        sig_shaped = signal.sosfilt(np.asarray([CR[0], RC[0]]), currents)
        plt.clf()
        plt.rcParams.update({'font.size': 18})
        plt.xlim(0, maxXrange)
        plt.plot(df_sig[x_axis_string_from_csv], sig_shaped)
        plt.xlabel('Time [ns]')
        plt.ylabel(ylabel)
        plt.title(plot_label)
        #plt.xticks(np.arange(0, maxXrange + 1, 10))
        plt.grid(True)
        print(label, " ", sig_shaped.max())
        plot_path = os.path.join(plotfolder, 'sig_shaped_' + label + '_tau_' + str(tau) + '.png')
        plt.savefig(plot_path.replace(" ", ""), bbox_inches='tight')
        dict_shapingTime_cell_peakShaper[str(tau)][label] = sig_shaped.max()
        column_string += "c"
        cells.append(label)
print("Plot written in %s/"%plotfolder)
cell_header = cell_header[:-1]

template_string_xtalk_table = """
\\begin{table}[h!]
\centering
\\begin{tabular}{c|COLUMN}
\stackunder{Cross-talk (\%)}{Shaping time (ns) $\downarrow$} & CELLS \\\\
\hline
BODY
\end{tabular}
\caption{Peak-to-peak ratio between victim and aggressor cells (cell 7) output.}
\label{xtalk}
\end{table}
""".replace("COLUMN", column_string).replace('CELLS', cell_header)

body_string = ""
for tau in taus:
    if tau == 0:
        body_string += "No shaper &"
    else:
        body_string += str(tau) + " &"
    for cell in cells:
        if cell == string_aggressor_cell:
            continue
        x_talk = round(dict_shapingTime_cell_peakShaper[str(tau)][cell] * 100 / dict_shapingTime_cell_peakShaper[str(tau)][string_aggressor_cell], 2)
        body_string += " " + str(x_talk) + " &"
    body_string = body_string[:-1] + " \\\\\n"
body_string = body_string[:-1]
#print(body_string)
template_string_xtalk_table = template_string_xtalk_table.replace('BODY', body_string)
print(template_string_xtalk_table)










