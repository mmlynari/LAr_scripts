#!/usr/bin/env python
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import signal, integrate
import os, glob, sys

input_file_dir = "/Users/brieucfrancois/Document/Fellowship/ElectrodesDesign/Prototype/V0/Measurements/Input_signals/"
input_file_regex = "C1--input-tower07cell07-output-tower07cell02-input*--00000.txt"
input_files = glob.glob(os.path.join(input_file_dir, input_file_regex))
for input_file in input_files:
    print(input_file)
    basename = os.path.basename(input_file)
    output_file = os.path.join(os.path.join(input_file_dir, basename.replace(".txt", "_forANSYS.csv")))
    ref_data = pd.read_csv(input_file, sep=',', skiprows=6, header=None)
    # switch to nanoseconds
    ref_data[3] = ref_data[3] * 1e9
    sampling_frequency = round(1 / ref_data[3].diff()[1])
    min_idx = ref_data[4][100:round((100+301)/ref_data[3].diff()[1])].idxmin()
    min_voltage = ref_data[4][min_idx]
    min_time = ref_data[3][min_idx]
    print("Start time = ", min_time)
    max_time = min_time + 301
    max_idx = min_idx + round(301/ref_data[3].diff()[1])
    ref_data[4] = ref_data[4] - min_voltage
    y_axis_values = []
    x_axis_values = []
    with open(output_file, 'w') as output:
        output.write("# ns V\n")
        current_idx = 0
        #first put zeros for 10 ns
        for idx in range(round(10/ref_data[3].diff()[1])):
            time = current_idx * ref_data[3].diff()[1]
            signal_height = 0
            x_axis_values.append(time)
            y_axis_values.append(signal_height)
            output.write("%f %f\n"%(time, signal_height))
            current_idx += 1
        #then put the real signal
        for idx in range(min_idx, max_idx):
            time = current_idx * ref_data[3].diff()[1]
            signal_height = ref_data[4][idx]
            x_axis_values.append(time)
            y_axis_values.append(signal_height)
            output.write("%f %f\n"%(time, signal_height))
            current_idx += 1
        #then put again 0's for 200 ns
        for idx in range(round(200/ref_data[3].diff()[1])):
            time = current_idx * ref_data[3].diff()[1]
            signal_height = 0
            x_axis_values.append(time)
            y_axis_values.append(signal_height)
            output.write("%f %f\n"%(current_idx * ref_data[3].diff()[1], 0))
            current_idx += 1
    print(output_file, " written.")

    plt.figure("Input signal")
    plt.plot(np.array(x_axis_values), np.array(y_axis_values))
    axes = plt.gca()
    #axes.set_xscale('log')
    #axes.set_xlim([0-step, rise_time + drift_time + step])
    #axes.set_xlim([p1x_rising - 5, drift_time + rise_time + 5])
    #axes.set_ylim([p1y_rising, p2y_rising + p2y_rising * 0.1])
    plt.xlabel('t [ns]')
    plt.ylabel('Input signal [V]')
    plt.grid()
    plt.title('Input signal from oscilloscope')
    plotFileName = os.path.join(input_file_dir, basename.replace(".txt", "_forANSYS.png"))
    plt.gcf().savefig(plotFileName)
    plt.clf()
    print(plotFileName, " written")
