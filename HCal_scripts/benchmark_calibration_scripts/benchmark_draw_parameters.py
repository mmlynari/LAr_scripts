import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
import os, sys, glob
import numpy as np
import argparse
from math import sqrt
from datetime import date
#from shutil import copy
from copy import copy

import gStyle

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")



parser = argparse.ArgumentParser()
## benchmark_21Nov_addUpstream_baseline
## benchmark_calib_addUpstreamParam_local 
parser.add_argument("-inputFiles", default = "/afs/cern.ch/user/m/mmlynari/workspace/FCC_MR_benchmark_16Apr24/LAr_scripts/FCCSW_ecal/benchmark_parameters_24Apr24/benchmark_*.root", help = "Regex for input files.", type = str)
parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d") + "_benchmark_parameters_baseline", help = "Postfix to append to the output folder.", type = str)
parser.add_argument("-color", default = 46, help = "Color of the graph", type = int)
parser.add_argument("-markerStyle", default = 21, help = "Style of the graph markers", type = int)
parser.add_argument("-cells", default = False, help = "Also produce plots with total energy from cells (no clustering) -- needs some refreshing", type = bool)
parser.add_argument("-clusterCollection", default = "CaloClusters", help = "Name of the cluster collection to use", type = str)
parser.add_argument("-lowEne", default = False, help = "Focus on low energy", type = bool)
parser.add_argument("-highEne", default = False, help = "Focus on high energy", type = bool)
parser.add_argument("-allEne", default = False, help = "Process all energies", type = bool)
parser.add_argument("-addHCal", default = False, help = "Free parameter 1 for HCal energy scaling", type = bool)
parser.add_argument("-below5", default = False, help = "Process all energies", type = bool)


args = parser.parse_args()

plot_dir_name = 'plots_performances_'+ args.outputPostfix
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

individual_plots = os.path.join(plot_dir_name, "resolution") 
if not os.path.isdir(individual_plots):
    os.mkdir(individual_plots)

inputFiles = glob.glob(args.inputFiles)
if not inputFiles:
    print("No input file found")
    sys.exit(1)

max_evt = 100000

dict_p0_error = {}
dict_p1_error = {}
dict_p2_error = {}
dict_p3_error = {}
dict_p4_error = {}
dict_p5_error = {}
energies_gev_float = []
energies_gev_float_err = [] 
p0 = []

for inputFile in inputFiles:
    print("Treating %s..."%inputFile)
    rootfile_path = inputFile
    energy = inputFile.split('_pMin_')[1].split("_")[0]
    energy_gev_float = float(energy)/1000.0
    energies_gev_float.append(energy_gev_float)
    energies_gev_float_err.append(0)
    energy = str(energy_gev_float).replace(".","dot")
    print('energy is ', energy)
    dict_p0_error[energy] = []
    dict_p1_error[energy] = []
    dict_p2_error[energy] = []
    dict_p3_error[energy] = []
    dict_p4_error[energy] = []
    dict_p5_error[energy] = []
    prefix = energy + "GeV_"

    f = ROOT.TFile(rootfile_path)
    #f.ls()
    histo = f.Get("parameters")
    if not histo:
        print("Failed to get data histogram")
        sys.exit(1)
    histo.SetDirectory(0)

    ## fill dictionaries with values stored in each histogram bin 
    dict_p0_error[energy].append(histo.GetBinContent(1))
    dict_p0_error[energy].append(histo.GetBinError(1))
    dict_p1_error[energy].append(histo.GetBinContent(2))
    dict_p1_error[energy].append(histo.GetBinError(2))
    dict_p2_error[energy].append(histo.GetBinContent(3))
    dict_p2_error[energy].append(histo.GetBinError(3))
    dict_p3_error[energy].append(histo.GetBinContent(4))
    dict_p3_error[energy].append(histo.GetBinError(4))
    dict_p4_error[energy].append(histo.GetBinContent(5))
    dict_p4_error[energy].append(histo.GetBinError(5))
    ## dict_p5_error[energy].append(histo.GetBinContent(6))
    ## dict_p5_error[energy].append(histo.GetBinError(6))    
    
    print("energy ", energy_gev_float, "parameter 1", histo.GetBinContent(1)) 

energies_gev_float.sort()
idx = 0
p0_vs_energy_graph = ROOT.TGraphErrors("p0_vs_energy_graph")
p1_vs_energy_graph = ROOT.TGraphErrors("p1_vs_energy_graph")
p2_vs_energy_graph = ROOT.TGraphErrors("p2_vs_energy_graph")
p3_vs_energy_graph = ROOT.TGraphErrors("p3_vs_energy_graph")
p4_vs_energy_graph = ROOT.TGraphErrors("p4_vs_energy_graph")
## p5_vs_energy_graph = ROOT.TGraphErrors("p5_vs_energy_graph")

for energy_float in energies_gev_float:
    energy_str = str(energy_float).replace(".", "dot")
    # fill graphs using dictionaries as inputs 
    p0_vs_energy_graph.SetPoint(idx, energy_float, dict_p0_error[energy_str][0])
    p0_vs_energy_graph.SetPointError(idx, 0, dict_p0_error[energy_str][1])

    p1_vs_energy_graph.SetPoint(idx, energy_float, dict_p1_error[energy_str][0])
    p1_vs_energy_graph.SetPointError(idx, 0, dict_p1_error[energy_str][1])

    p2_vs_energy_graph.SetPoint(idx, energy_float, dict_p2_error[energy_str][0])
    p2_vs_energy_graph.SetPointError(idx, 0, dict_p2_error[energy_str][1])

    p3_vs_energy_graph.SetPoint(idx, energy_float, dict_p3_error[energy_str][0])
    p3_vs_energy_graph.SetPointError(idx, 0, dict_p3_error[energy_str][1])

    p4_vs_energy_graph.SetPoint(idx, energy_float, dict_p4_error[energy_str][0])
    p4_vs_energy_graph.SetPointError(idx, 0, dict_p4_error[energy_str][1])

    ## p5_vs_energy_graph.SetPoint(idx, energy_float, dict_p5_error[energy_str][0])
    ## p5_vs_energy_graph.SetPointError(idx, 0, dict_p5_error[energy_str][1])

    idx +=1

ROOT.gStyle.SetOptFit(000)

## for each parameter we use different fitting functions - still can be improved by finding better functions 
def plot_parameter_and_fit(parameter_name, param_vs_energy_graph, write_formula = False):

    #tf1 = ROOT.TF1("tf1", "tf1", energies_gev_float[0], energies_gev_float[-1])

    if args.lowEne:
        if parameter_name == 'p0':
            y_axis_label = "parameter 0"
            tf1 = ROOT.TF1("tf1", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomLeft_legend
        if parameter_name == 'p1':
            y_axis_label = "parameter 1"
            tf1 = ROOT.TF1("tf1", "[0]+[1]/(x*x+[2])", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomRight_legend
        if parameter_name == 'p2':
            y_axis_label = "parameter 2"
            if args.addHCal:
                tf1 = ROOT.TF1("tf1", "[0]+[1]*(x+[2])**2", energies_gev_float[0], energies_gev_float[-1])
            else:
                tf1 = ROOT.TF1("tf1", "[0]+[1]/(x+[2])**2", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomRight_legend
        if parameter_name == 'p3':
            y_axis_label = "parameter 3"
            tf1 = ROOT.TF1("tf1", "[0]+[1]/(x)", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomRight_legend
        if parameter_name == 'p4':
            y_axis_label = "parameter 4"
            # [0]+[1]/(x*x+[2])
            tf1 = ROOT.TF1("tf1", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.topLeft_legend
        ## constant parameter
        if parameter_name == 'p5':
            y_axis_label = "parameter 5" 
            tf1 = ROOT.TF1("tf1", "[0]+[1]*x+[2]*x*x", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomLeft_legend

    if args.highEne or args.allEne:
        if parameter_name == 'p0':
            y_axis_label = "parameter 0"
            if args.below5: 
                tf1 = ROOT.TF1("tf1", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", energies_gev_float[0], energies_gev_float[-1])
                legend_style = gStyle.bottomLeft_legend
            else:
                tf1 = ROOT.TF1("tf1", "[0]+[1]/sqrt(x)", energies_gev_float[0], energies_gev_float[-1])
                legend_style = gStyle.topRight_legend
        if parameter_name == 'p1':
            y_axis_label = "parameter 1"
            tf1 = ROOT.TF1("tf1", "[0]+[1]/x+[2]*x", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomRight_legend
        if parameter_name == 'p2':
            y_axis_label = "parameter 2"
            if args.addHCal and args.below5:
                tf1 = ROOT.TF1("tf1", "[0]+[1]*x+[2]*x*x", energies_gev_float[0], energies_gev_float[-1])
                legend_style = gStyle.topRight_legend
            else: 
                tf1 = ROOT.TF1("tf1", "[0]+[1]/(x+[2])", energies_gev_float[0], energies_gev_float[-1])
                legend_style = gStyle.bottomRight_legend
        if parameter_name == 'p3':
            y_axis_label = "parameter 3"
            tf1 = ROOT.TF1("tf1", "[0]+[1]/(x)", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomRight_legend
        if parameter_name == 'p4':
            y_axis_label = "parameter 4"
            tf1 = ROOT.TF1("tf1", "[0]+[1]/(x+[2])", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomRight_legend
        ## constant parameter
        if parameter_name == 'p5':
            y_axis_label = "parameter 5" 
            tf1 = ROOT.TF1("tf1", "[0]+[1]*x+[2]*x*x", energies_gev_float[0], energies_gev_float[-1])
            legend_style = gStyle.bottomLeft_legend

    tf1.SetLineColor(args.color)
    canvas_param_vs_energy = ROOT.TCanvas(parameter_name + "_vs_energy", parameter_name + "_vs_energy")
    param_vs_energy_legend = copy(legend_style)
    param_vs_energy_legend.AddEntry(ROOT.nullptr,"#bf{FCC-ee simulation}","")
    param_vs_energy_graph.SetMarkerStyle(21);
    param_vs_energy_graph.GetYaxis().SetTitle(y_axis_label)
    param_vs_energy_graph.GetXaxis().SetTitle("Energy [GeV]")
    param_vs_energy_fit = param_vs_energy_graph.Fit(tf1, "SQ")
    a = str(round(param_vs_energy_fit.Get().Parameter(0), 4))
    b = str(round(param_vs_energy_fit.Get().Parameter(1), 4))
    c = str(round(param_vs_energy_fit.Get().Parameter(2), 4))
    d = str(round(param_vs_energy_fit.Get().Parameter(3), 4))
    ''' 
    if a>1.:
        a = str(round(param_vs_energy_fit.Get().Parameter(0), 2))
    if b>1.:
        b = str(round(param_vs_energy_fit.Get().Parameter(1), 2))
    if c>1.:
        c = str(round(param_vs_energy_fit.Get().Parameter(2), 2))
    if d>1.:
        d = str(round(param_vs_energy_fit.Get().Parameter(3), 2))
    ''' 


    if args.lowEne: 
        if parameter_name == 'p0':
            param_vs_e_formula = "#color[%d]{%s #plus %s*E #plus %s*E^2 #plus %s*E^3}"%(args.color, a, b, c, d)
        if parameter_name == 'p1':
            param_vs_e_formula = "#color[%d]{%s #plus #frac{%s}{E^2+%s}}"%(args.color, a, b, c)
        if parameter_name == 'p2':
            if args.addHCal:
                param_vs_e_formula = "#color[%d]{%s #plus %s*(E+%s)^2}"%(args.color, a, b, c)
            else: 
                param_vs_e_formula = "#color[%d]{%s #plus %s/(E+%s)^2}"%(args.color, a, b, c)
        if parameter_name == 'p3':
            param_vs_e_formula = "#color[%d]{%s #plus #frac{%s}{E}}"%(args.color, a, b)
        if parameter_name == 'p4':
            param_vs_e_formula = "#color[%d]{%s #plus %s*E #plus %s*E^2 #plus %s*E^3}"%(args.color, a, b, c, d)

    if args.highEne or args.allEne: 
        if parameter_name == 'p0':
            if args.below5: 
                param_vs_e_formula = "#color[%d]{%s #plus %s*E #plus %s*E^2 #plus %s*E^3}"%(args.color, a, b, c, d)
            else:
                param_vs_e_formula = "#color[%d]{%s #plus #frac{%s}{#sqrt{E}}}"%(args.color, a, b)
        if parameter_name == 'p1': 
            param_vs_e_formula = "#color[%d]{%s #plus #frac{%s}{E} #plus %s*E}"%(args.color, a, b, c)
        if parameter_name == 'p2':
            if args.addHCal and args.below5:
                param_vs_e_formula = "#color[%d]{%s #plus %s*E #plus %s*E^2}"%(args.color, a, b, c)
            else: 
                param_vs_e_formula = "#color[%d]{%s #plus #frac{%s}{E+%s}}"%(args.color, a, b, c) 
        if parameter_name == 'p3':
            param_vs_e_formula = "#color[%d]{%s #plus #frac{%s}{E}}"%(args.color, a, b)
        if parameter_name == 'p4':
            param_vs_e_formula = "#color[%d]{%s #plus #frac{%s}{E+%s}}"%(args.color, a, b, c)

        ''' 
    if parameter_name == 'p5':
        b = str(round(param_vs_energy_fit.Get().Parameter(1), 4))
        c = str(round(param_vs_energy_fit.Get().Parameter(2), 4))
        param_vs_e_formula = "#color[%d]{%s #plus %s*E #plus %s*E^2}"%(args.color, a, b, c) 
        '''
    if write_formula: 
        param_vs_energy_legend.AddEntry(ROOT.nullptr, param_vs_e_formula, "")
    param_vs_energy_graph.Draw("ap")
    param_vs_energy_legend.Draw()
    param_vs_energy_graph.Write()
    canvas_param_vs_energy.SetLogx()
    canvas_param_vs_energy.Print(os.path.join(plot_dir_name,parameter_name + "_vs_energy.png"))
    #canvas_param_vs_energy.Print(os.path.join("p0_vs_energy.png"))
    canvas_param_vs_energy.Write()


output_rootfile_graph_name = os.path.join(plot_dir_name, "parameters_vs_energy.root")
print("Writing %s"%output_rootfile_graph_name)
output_rootfile_graph = ROOT.TFile(output_rootfile_graph_name, "recreate")

p0_vs_energy_fit = plot_parameter_and_fit('p0', p0_vs_energy_graph, write_formula = True)
#p1_vs_energy_fit = plot_parameter_and_fit('p1', p1_vs_energy_graph, write_formula = True)
p2_vs_energy_fit = plot_parameter_and_fit('p2', p2_vs_energy_graph, write_formula = True)
p3_vs_energy_fit = plot_parameter_and_fit('p3', p3_vs_energy_graph, write_formula = True)
p4_vs_energy_fit = plot_parameter_and_fit('p4', p4_vs_energy_graph, write_formula = True)
## p5_vs_energy_fit = plot_parameter_and_fit('p5', p5_vs_energy_graph, write_formula = True)

## mess below can be ignored 
'''
canvas_p0_vs_energy = ROOT.TCanvas("p0_vs_energy")
p0_vs_energy_legend = copy(gStyle.topRight_legend_relEresol)
p0_vs_energy_legend.AddEntry(ROOT.nullptr,"#bf{FCC-ee simulation}","")
p0_vs_energy_graph.SetMarkerStyle(21);
p0_vs_energy_graph.GetYaxis().SetTitle("parameter 0")
p0_vs_energy_graph.GetXaxis().SetTitle("Energy [GeV]")
p0_vs_energy_fit = p0_vs_energy_graph.Fit(tf1_p0, "SQ")
print('p0: fit parameter 0 ', p0_vs_energy_fit.Get().Parameter(0))
print('p0: fit parameter 1 ', p0_vs_energy_fit.Get().Parameter(1))
a_0 = str(round(p0_vs_energy_fit.Get().Parameter(0), 2))
b_0 = str(round(p0_vs_energy_fit.Get().Parameter(1), 2))
p0_vs_e_formula = "#color[%d]{#frac{%s}{#sqrt{E}} #plus %s}"%(args.color, a_0, b_0) 
p0_vs_energy_legend.AddEntry(ROOT.nullptr, p0_vs_e_formula, "")
p0_vs_energy_graph.Draw("ap")
p0_vs_energy_legend.Draw()
p0_vs_energy_graph.Write()
canvas_p0_vs_energy.Print(os.path.join("p0_vs_energy.png"))
canvas_p0_vs_energy.Write()

canvas_p2_vs_energy = ROOT.TCanvas("p2_vs_energy")
p2_vs_energy_legend = copy(gStyle.bottomRight_legend)
p2_vs_energy_graph.SetMarkerStyle(21);
p2_vs_energy_graph.GetYaxis().SetTitle("parameter 1")
p2_vs_energy_graph.GetXaxis().SetTitle("Energy [GeV]")
p2_vs_energy_fit = p2_vs_energy_graph.Fit(tf1_p3, "SQ")
print('p2: fit parameter 0 ', p2_vs_energy_fit.Get().Parameter(0))
print('p2: fit parameter 1 ', p2_vs_energy_fit.Get().Parameter(1))
a_2 = str(round(p2_vs_energy_fit.Get().Parameter(0), 2))
b_2 = str(round(p2_vs_energy_fit.Get().Parameter(1), 2))
p2_vs_e_formula = "#color[%d]{#frac{%s}{log(E)} #plus %s}"%(args.color, a_2, b_2) 
p2_vs_energy_legend.AddEntry(ROOT.nullptr, p2_vs_e_formula, "")
p2_vs_energy_graph.Draw("ap")
p2_vs_energy_legend.Draw()
p2_vs_energy_graph.Write()
canvas_p2_vs_energy.Print(os.path.join("p2_vs_energy.png"))
canvas_p2_vs_energy.Write()

canvas_p3_vs_energy = ROOT.TCanvas("p3_vs_energy")
p3_vs_energy_legend = copy(gStyle.bottomRight_legend)
p3_vs_energy_graph.SetMarkerStyle(21);
p3_vs_energy_graph.GetYaxis().SetTitle("parameter 2")
p3_vs_energy_graph.GetXaxis().SetTitle("Energy [GeV]")
p3_vs_energy_fit = p3_vs_energy_graph.Fit(tf1_p3, "SQ")
print('p3: fit parameter 0 ', p3_vs_energy_fit.Get().Parameter(0))
print('p3: fit parameter 1 ', p3_vs_energy_fit.Get().Parameter(1))
a_3 = str(round(p3_vs_energy_fit.Get().Parameter(0), 2))
b_3 = str(round(p3_vs_energy_fit.Get().Parameter(1), 4))
p3_vs_e_formula = "#color[%d]{#frac{%s}{E} #plus %s}"%(args.color, a_3, b_3) 
p3_vs_energy_legend.AddEntry(ROOT.nullptr, p3_vs_e_formula, "")
p3_vs_energy_graph.Draw("ap")
p3_vs_energy_legend.Draw()
p3_vs_energy_graph.Write()
canvas_p3_vs_energy.Print(os.path.join("p3_vs_energy.png"))
canvas_p3_vs_energy.Write()

''' 

output_rootfile_graph.Close()


