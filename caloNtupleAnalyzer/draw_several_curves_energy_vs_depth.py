import os, sys, argparse
import ROOT
from datetime import date

import gStyle


parser = argparse.ArgumentParser()
parser.add_argument("-inputFiles", default = ["plots_performances_210201/relResol_vs_energy.root"], help = "List of input files.", nargs="+")
parser.add_argument("-labels", default = [], help = "List of labels for the legend, same order as input files.", nargs="+")
parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d"), help = "Postfix to append to the output folder.", type = str)
args = parser.parse_args()

plot_dir_name = 'plots_combined_energy_vs_depth_'+ args.outputPostfix
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

legend = gStyle.topRight_legend
colors = gStyle.colors
if len(colors) < len(args.inputFiles):
    print("Not enough color to draw all the asked curved")
    sys.exit(1)

dict_legend_tprof = {}
canvas = ROOT.TCanvas('energy_vs_depth', 'energy_vs_depth')

count = 0
energy_legends = []
max_y = 0 
for inputFile in args.inputFiles:
    print("Treating ", inputFile)
    energy_str = str(inputFile.split('_')[-1].split('.')[0])
    energy = float(energy_str.replace('dot', '.'))
    unit = " GeV"
    if energy < 1:
        energy = energy * 1000
        unit = " MeV"
    energy_legend = str(int(energy)) + unit
    if len(args.labels) != 0:
        energy_legend = args.labels[count]
    energy_legends.append(energy_legend)

    inputRootFile = ROOT.TFile(inputFile, "r")
    #inputRootFile.cd()
    name = "tprof_energy_vs_depth_" + energy_str
    #graph = ROOT.TProfile(inputRootFile.Get(name))
    #graph = inputRootFile.Get(name)
    #graph = ROOT.TProfile(inputRootFile.GetObject(name.c_str()))
    graph = ROOT.TProfile(inputRootFile.Get(name))
    #graph.SetMaximum(1.2*graph.GetMaximum())
    #graph.SetDirectory(0)
    print("tprof_energy_vs_dept_"+energy_str)
    dict_legend_tprof[energy_legend] = graph
    dict_legend_tprof[energy_legend].SetDirectory(0)
    dict_legend_tprof[energy_legend].SetLineColor(colors[count])
    dict_legend_tprof[energy_legend].SetTitle("Avergage energy deposit vs depth")
    #legend.AddEntry(graph, energy_legend)
    if count == 0:
        dict_legend_tprof[energy_legend].Draw()
    else:
        dict_legend_tprof[energy_legend].Draw('same')
    count += 1

for energy_legend in energy_legends:
    print(dict_legend_tprof[energy_legend])
    legend.AddEntry(dict_legend_tprof[energy_legend], energy_legend)

legend.Draw()
canvas.Print(os.path.join(plot_dir_name, 'tprof_energy_vs_dept.png'))


