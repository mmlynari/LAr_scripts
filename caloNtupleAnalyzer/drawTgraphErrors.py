import os, sys
import ROOT
from datetime import date

import gStyle


parser = argparse.ArgumentParser()
parser.add_argument("-inputFiles", default = ["plots_performances_210201/relResol_vs_energy.root"], help = "List of input files.", nargs="+")
parser.add_argument("-labels", default = ["12 layers LAr"], help = "List of labels for the legend, same order as input files.", nargs="+")
parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d"), help = "Postfix to append to the output folder.", type = str)
args = parser.parse_args()

plot_dir_name = 'plots_combined_performances_'+ args.outputPostfix
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

legend = gStyle.topRight_legend


for inputFile in inputFiles:
    inputRootFile = ROOT.TFile(inputFile, "r")
    graph = inputRootFile.Get("Graph")
    graph.Draw()

