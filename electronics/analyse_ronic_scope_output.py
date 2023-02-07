#!/usr/bin/env python
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import signal, integrate
import os, glob, sys

input_data_folder = "/Users/brieucfrancois/Document/Fellowship/ElectrodesDesign/Prototype/V0/Measurements/paris/only_v3/" #RefCurve_input_tower02Cell07_output_tower02Cell06_50OhmInj_burst_goodOne.csv
data_file_pattern = "RefCurve_input_tower04cell07*_v3.Wfm.csv"
#data_file_pattern = "*input-tower08cell*bypass*.txt"
#data_file_pattern = "C1--input-tower08cell15-output-tower08cell15-1megaohm-better--00000.txt"
#data_file_pattern = "C1--input-tower08cell08-output-tower08cell08-1megaohm--00000.txt"
#output_plot_folder = os.path.join(input_data_folder, "plots_trialWithStrips")
output_plot_folder = os.path.join(input_data_folder, "plots_paris_50OhmInjector_2Shields_injectOnCell7")
if not os.path.isdir(output_plot_folder):
    os.mkdir(output_plot_folder)

data_files = glob.glob(os.path.join(input_data_folder, data_file_pattern))
if len(data_files) == 0:
    print("Error: no data file found in %s"%os.path.join(input_data_folder, data_file_pattern))
    exit(1)
file_reference_signal = ""
n_ref = 0
for data_file in data_files:
    words = os.path.basename(data_file).split("_")
    if len(words) < 11:
        continue
    if(words[2] == words [4]):
        file_reference_signal = data_file
        ref_tower = str(int(words[2].split('tower')[1][0:2]))
        ref_cell = str(int(words[2].split('cell')[1][0:2]))
        n_ref += 1
    print(words)

if file_reference_signal == "":
    print("Error: no reference file found to define the time window of a single pulse!")
    sys.exit(1)
#if n_ref > 1:
#    print("Error: two reference files found to define the time window of a single pulse! Please change the criteria to get only one ref file.")
#    sys.exit(1)

print("Reference file used: ", file_reference_signal)

ref_data = pd.read_csv(file_reference_signal, sep=',', skiprows=0, header=None)
time_step = 0.1  # ns
n_points = 5000
time_series = [x * time_step for x in range(0, n_points)]
time_data = pd.DataFrame(time_series)
# switch to nanoseconds
#ref_data[3] = ref_data[3] * 1e9
ref_data[1] = time_data
sampling_frequency = round(1 / time_step)
# plot the raw signals
plt.clf()
ref_data.plot(x = 1, y = 0, xlabel = "Time (ns)", ylabel = "Voltage [V]", grid = True)
plot_path = os.path.join(output_plot_folder, 'raw_main_signal.png')
plt.savefig(plot_path)
print(plot_path, " written")
print("Sampling frequency per nano second: ", sampling_frequency)
min_idx = 0
min_voltage = ref_data[0][min_idx]
min_time = time_data[0][min_idx]
print("Start time = ", min_time)
max_time = time_step * (n_points - 1)
max_idx = n_points - 1
print("Stop time = ", time_data[0][max_idx])
print("Time duration: ", time_data[0][max_idx] - time_data[0][min_idx])
peak_index = ref_data[0][min_idx:max_idx].idxmax()
peak_voltage = ref_data[0][peak_index]
peak_time = time_data[0][peak_index]
signal_heigt = peak_voltage - min_voltage
print("Signal height = ", signal_heigt)
print("Peaking time = ", peak_time - min_time)
input_tower_type = ""
if ref_tower == str(2):
    input_tower_type = " baseline "
if ref_tower == str(8):
    input_tower_type = " no shield "
if ref_tower == str(7):
    input_tower_type = " one shield "
if ref_tower == str(4):
    input_tower_type = " shield cell7/8 + common GND "
if ref_tower == str(5):
    input_tower_type = " GND plate "
if ref_tower == str(6):
    input_tower_type = " no hv via "
if ref_tower == str(9):
    input_tower_type = " doubled width shields "
if ref_tower == str(10):
    input_tower_type = " halved width shields "
if ref_tower == str(11):
    input_tower_type = " one shield doubled width "
if ref_tower == str(12):
    input_tower_type = " 70 Ohms "
if ref_tower == str(13):
    input_tower_type = " outer radius cells 6 and 7 "
if ref_tower == str(14):
    input_tower_type = " GND btw strips "
if ref_tower == str(15):
    input_tower_type = " baseline "
if ref_tower == str(16):
    input_tower_type = " baseline "

# from the main signal, find the start/end minmax, extract those time from the dataframes and pad the rest with zeros

x_axis_string_from_csv = "Time"
string_aggressor_cell = "Output"
input_string = ""
output_string = ""

cell_header = ""
cell_headers = []
column_string = ""

dict_shapingTime_towercellstring_peakShaper = {}
dict_shapingTime_ref_peakShaper = {}
cells = []
for data_file in data_files:
    main_signal = False
    print("Treating ", data_file)
    if os.path.getsize(data_file) == 0:
        print("\t Empty file, skipping!")
        continue
    words = os.path.basename(data_file).split("_")
    if len(words) >= 11:
        input_string = words[2]
        output_string = words [4]
        if(input_string == output_string):
            main_signal = True
    else:
        input_string = "input"
        output_string = "noOutput"
    input_tower = str(int(input_string.split('tower')[1][0:2]))
    input_cell = str(int(input_string.split('cell')[1][0:2]))
    output_tower = str(int(output_string.split('tower')[1][0:2]))
    output_cell = str(int(output_string.split('cell')[1][0:2]))
    cells.append(output_tower + "-" + output_cell)
    df_sig = pd.read_csv(data_file, sep=',', header=None)
    df_sig[1] = time_data
    print(df_sig)
    # switch to nanoseconds, and shift to 0
    #df_sig[3] = df_sig[3] * 1e9 - min_time
    start_value = df_sig[0][0]
    #start_value = df_sig[4][0:10].mean()

    df_sig[0] = df_sig[0] - start_value

    #df_sig = df_sig.append(zeros) do this only if you want to see what enters the filter, otherwise the plots shows the zero padded curve
    n_columns = df_sig.shape[1]

    plt.clf()
    df_sig.plot(x = 1, y = 0, xlabel = "Time (ns)", ylabel = "Voltage [V]", grid = True)
    #df_sig.plot(x = x_axis_string_from_csv, ylabel = "Current [nA]", subplots = True, xlim = (0, my_timeLimit), layout = (int(n_columns/2.0) + 1, 2), figsize = (12, 12), grid = True)
    #integral = integrate.trapz(df_sig[4], x = df_sig[3])
    #print("Integral: ", integral)
    label = os.path.basename(data_file)
    plot_path = os.path.join(output_plot_folder, '%s.png'%(label))
    plt.savefig(plot_path)
    print(plot_path, " written")

    if not (input_tower == output_tower and input_cell == output_cell):
        column_string += "c"
        if not input_tower == output_tower:
            cell_headers.append("T " + output_tower + " C " + output_cell + " (nbr) &")
        else:
            cell_headers.append("T " + output_tower + " C " + output_cell + " &")

    plot_label = "Input signal on cell " + input_cell + " tower " + input_tower + "(" +  input_tower_type + ")" + "\n" + "Output signal on cell " + output_cell + " tower " + output_tower

    #taus = [0, 20, 50, 100, 150, 200, 300]
    taus = [0, 20, 50, 100, 150]
    for tau in taus:
        extra_plot_label = "\n Shaping time: " + str(tau) + " ns"
        #print("Shaping time: ", tau)
        # define CR-RC^2 filter
        # order of the filter, critical (cut-off) frequency of the filter, filter type, analog or digital, type of output (pole zero, num/dem or second order section), sampling frequency
        if tau == 0:# means no shaper
            ylabel = "Voltage at PCB output port [V]"
            extra_title = ", no shaper"
            CR = signal.butter(0, 1, btype='highpass', analog=False, output='sos', fs=sampling_frequency)
            RC = signal.butter(0, 1, btype='lowpass', analog=False, output='sos', fs=sampling_frequency)
            maxXrange = max_time
        else:
            ylabel = "Shaper output (arbitrary unit)"
            extra_title = ", shaping time " + str(tau) + " ns"
            CR = signal.butter(1, 1/(2*np.pi*tau), btype='highpass', analog=False, output='sos', fs=sampling_frequency)
            RC = signal.butter(2, 1/(2*np.pi*tau), btype='lowpass', analog=False, output='sos', fs=sampling_frequency)
            maxXrange = min(25 * tau, max_time)
        # add zeros after the signal so that the filter behave well, 0th component will be voltage set to min_value, 1th component will be time
        nrows_zero_padding = round((6 * tau)/df_sig[1].diff()[1])
        times_zero_padding = [max_time + (x * time_step) for x in range(1, nrows_zero_padding + 1)]
        values_zero_padding = [df_sig[0][0]] * nrows_zero_padding 
        zeros = pd.DataFrame(values_zero_padding)
        zeros[1] = times_zero_padding
        print(zeros)
        # add zeros before the signal (looks like it is not needed)
        #zeros_before = zeros.copy()
        #zeros_before[3] = zeros_before[3] - ((max_idx - min_idx) * df_sig[3].diff()[1] + round(6 * tau))

        #save the signal feeding the shaper
        plt.clf()
        plt.plot(df_sig.append(zeros)[1], df_sig.append(zeros)[0])
        #plt.plot(zeros_before.append(df_sig.append(zeros))[3], zeros_before.append(df_sig.append(zeros))[4])
        plt.xlabel('Time [ns]')
        plt.ylabel(ylabel)
        plt.title(plot_label + extra_plot_label)
        #plt.xticks(np.arange(0, maxXrange + 1, 10))
        plt.grid(True)
        #print(label, " ", sig_shaped.max())
        plot_path = os.path.join(output_plot_folder, 'sig_forShaper_' + label + '_tau_' + str(tau) + '.png')
        plt.savefig(plot_path.replace(" ", ""), bbox_inches='tight')

        #sig_shaped = signal.sosfilt(np.asarray([CR[0], RC[0]]), zeros_before.append(df_sig.append(zeros))[4])
        sig_shaped = signal.sosfilt(np.asarray([CR[0], RC[0]]), df_sig.append(zeros)[0])

        plt.clf()
        #plt.rcParams.update({'font.size': 18})
        #plt.xlim(0, maxXrange - min_time)
        plt.plot(df_sig.append(zeros)[1], sig_shaped)
        #plt.plot(zeros_before.append(df_sig.append(zeros))[3], sig_shaped)
        plt.xlabel('Time [ns]')
        plt.ylabel(ylabel)
        plt.title(plot_label + extra_plot_label)
        #plt.xticks(np.arange(0, maxXrange + 1, 10))
        plt.grid(True)
        #print(label, " ", sig_shaped.max())
        plot_path = os.path.join(output_plot_folder, 'sig_afterShaper_' + label.replace(".txt", "") + '_tau_' + str(tau) + '.png')
        plt.savefig(plot_path.replace(" ", ""), bbox_inches='tight')
        # the division by four in the signal range to check maximum is to make sure we do not select a second maximum
        try:
            #dict_shapingTime_towercellstring_peakShaper[str(tau)][output_tower+output_cell] = sig_shaped.max()
            dict_shapingTime_towercellstring_peakShaper[str(tau)][output_tower+output_cell] = sig_shaped[0:round(len(sig_shaped)/4.0)].max()
        except:
            dict_shapingTime_towercellstring_peakShaper[str(tau)] = {}
            #dict_shapingTime_towercellstring_peakShaper[str(tau)][output_tower+output_cell] = sig_shaped.max()
            dict_shapingTime_towercellstring_peakShaper[str(tau)][output_tower+output_cell] = sig_shaped[0:round(len(sig_shaped)/4.0)].max()
        if(main_signal):
            dict_shapingTime_ref_peakShaper[str(tau)] = sig_shaped.max()

for tmp_cell_header in sorted(cell_headers):# careful that this does not put the wrong header on top of the wrong column !! --> make sure cells and cell_headers will be sorted the same way
    cell_header += tmp_cell_header

cell_header = cell_header[:-1]
template_string_xtalk_table = """
\\begin{table}[h!]
\centering
\\begin{tabular}{c|COLUMN}
\stackunder{Cross-talk (\%)}{Shaping time (ns) $\downarrow$} & CELLS \\\\
\hline
BODY
\end{tabular}
\caption{Peak-to-peak ratio between victim and aggressor cells (cell NUM) output, \\textbf{TOWERTYPE tower}.}
\label{xtalk}
\end{table}
""".replace("COLUMN", column_string).replace('CELLS', cell_header).replace('TOWERTYPE', input_tower_type).replace("NUM", ref_cell)

body_string = ""
for tau in taus:
    if tau == 0:
        body_string += "No shaper &"
    else:
        body_string += str(tau) + " &"
    for cell in sorted(cells):
        current_tower = cell.split("-")[0]
        current_cell = cell.split("-")[1]
        if current_tower == ref_tower and current_cell == ref_cell:
            continue
        #x_talk = round(dict_shapingTime_towercellstring_peakShaper[str(tau)][cell] * 100 / dict_shapingTime_towercellstring_peakShaper[str(tau)][string_aggressor_cell], 2)
        #print(cell, " ", tau, " ", dict_shapingTime_towercellstring_peakShaper[str(tau)][cell.replace("-", "")])
        x_talk = round(dict_shapingTime_towercellstring_peakShaper[str(tau)][cell.replace("-", "")] * 100 / dict_shapingTime_ref_peakShaper[str(tau)], 2)
        body_string += " " + str(x_talk) + " &"
    body_string = body_string[:-1] + " \\\\\n"
body_string = body_string[:-1]
template_string_xtalk_table = template_string_xtalk_table.replace('BODY', body_string)

with open(os.path.join(output_plot_folder, "table.txt"), "a") as texFile:
    texFile.write(template_string_xtalk_table)
print(template_string_xtalk_table)
print(os.path.join(output_plot_folder, "table.txt"), " written.")


