import ROOT
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

ROOT.gROOT.SetBatch(ROOT.kTRUE)

plotfolder = "pileup_studies"
if not os.path.isdir(plotfolder):
    os.mkdir(plotfolder)


generator = ROOT.TRandom3()
collision_rate = 100e3 # Hz
bunch_crossing_rate = 50e6 # Hz
bunch_crossing_timeStep = 1 / bunch_crossing_rate# seconds
collision_probability = collision_rate / bunch_crossing_rate

print(collision_probability)

signal_durations = np.linspace(0, 2000e-9, 400) #s
n_overlaps = []

# Out of time pile up from uniform distribution of collisions over one second
#random_numbers = []
#for i in range(int(collision_rate)):
#    time = generator.Uniform(0, 1e9)
#    print(time)
#    random_numbers.append(time)
#
#print(random_numbers)
#times_array = np.array(random_numbers)

#from scipy.spatial import distance
#overlap_time = 1500 # ns
#for time in random_numbers:
#    print(distance.cdist(np.array([time]), times_array).min(axis=1))
 
# draw a TH1 with number of occurence where they are closer than X

# Out of time pile up from "collision or not collision" per bx
#for second in range(0, 10):
last_collision_bx = 0
n_collision = 0
max_times_between_collisions = 10000
th1_times_between_two_collisions = ROOT.TH1F("Time between two collisions", "Time between two collisions", 100, 0, max_times_between_collisions)
th1_times_between_two_collisions.GetYaxis().SetTitle("Entries")
th1_times_between_two_collisions.GetXaxis().SetTitle("Time between consecutive collisions [ns]")
n_overlaps = [0 for idx in signal_durations]
for bx_id in range(int(bunch_crossing_rate)):
#for bx_id in range(int(100000)):
    # emulate the collision probability
    rndm_number = generator.Uniform(0, 1)
    if rndm_number < collision_probability:
        idx_signal_duration = 0
        time_between_two_collisions = (bx_id - last_collision_bx) * bunch_crossing_timeStep
        th1_times_between_two_collisions.Fill(time_between_two_collisions * 1e9)
        for signal_duration in signal_durations:
            if time_between_two_collisions < signal_duration:
                n_overlaps[idx_signal_duration] += 1
                #print('overlap')
                #print((bx_id - last_collision_bx) * bunch_crossing_timeStep * 1e9)
            idx_signal_duration += 1
        last_collision_bx = bx_id
        n_collision += 1
print("N_bx: ", len(range(int(bunch_crossing_rate))))
print("n_collision: ", n_collision)
print("n_overlap: ", n_overlaps)
print("n_overlap length: ", len(n_overlaps))
print("n_signal duration: ", len(signal_durations))

print(th1_times_between_two_collisions.GetMean())

print("Trasform n_overlap in array")
overlap_percentage = []
for n_overlap in n_overlaps:
    overlap_percentage.append(n_overlap * 100 / n_collision)

print("Trasform overlap_percentage in array")
array_overlap_percentage = np.array(overlap_percentage)
print("Trasform signal_duration in array")
array_signal_durations = np.array([signal_duration * 1e9 for signal_duration in signal_durations])

print("Plot things")
plt.clf()
plt.rcParams.update({'font.size': 13})
#plt.xlim(0, maxXrange)
#print(array_signal_durations)
#print(array_overlap_percentage)
m, b = np.polyfit(array_signal_durations, array_overlap_percentage, 1)
plot = plt.plot(array_signal_durations, array_overlap_percentage, 'o')
plt.plot(array_signal_durations, m * array_signal_durations + b)
#plt.figure().set_yscale('log')
plt.plot()

plt.xlabel('Signal duration [ns]')
plt.ylabel("Percentage of events with overlap (%)")
plt.text(10, 16.5, "Bunch crossing rate: %d MHz\nPhysics collision rate: %d kHz"%(int(bunch_crossing_rate/1e6), int(collision_rate/1e3)))
plt.text(1050, 7.5, "%f * x + %f"%(m, b))
#plt.title(plot_label)
#plt.xticks(np.arange(0, maxXrange + 1, 10))
plt.grid(True)
plot_path = os.path.join(plotfolder, 'events_with_overlap_vs_signal_duration.png')
plt.savefig(plot_path)

plt.clf()
fig = plt.figure()
ax = fig.add_subplot()
ax.set_yscale('log')
plt.rcParams.update({'font.size': 13})
#plt.xlim(0, maxXrange)
#print(array_signal_durations)
#print(array_overlap_percentage)
plot = plt.plot(array_signal_durations, array_overlap_percentage, 'o')
plt.plot(array_signal_durations, m * array_signal_durations + b)
#plt.figure().set_yscale('log')
plt.plot()

plt.xlabel('Signal duration [ns]')
plt.ylabel("Percentage of events with overlap (%)")
plt.text(750, 1, "Bunch crossing rate: %d MHz\nPhysics collision rate: %d kHz"%(int(bunch_crossing_rate/1e6), int(collision_rate/1e3)))
plt.text(1050, 4, "%f * x + %f"%(m, b))
#plt.title(plot_label)
#plt.xticks(np.arange(0, maxXrange + 1, 10))
plt.grid(True)
plot_path = os.path.join(plotfolder, 'events_with_overlap_vs_signal_duration_logy.png')
plt.savefig(plot_path)

canvas = ROOT.TCanvas("Time between two collisions", "Time between two collisions")
#f1 = ROOT.TF1("f1", "[0]*TMath::Poisson(x,[1])", 0, max_times_between_collisions) 
#f1 = ROOT.TF1("f1","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, 4000) #"xmin" = 0, "xmax" = 10
#f1.SetParameters(1, 4000) #you MUST set non-zero initial values for parameters
#th1_times_between_two_collisions.Fit("f1", "R") #"R" = fit between "xmin" and "xmax" of the "f1"
th1_times_between_two_collisions.Draw()
canvas.Print(os.path.join(plotfolder, "time_between_two_collisions.png"))












