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

print("Launch without sourcing k4 environment!!")

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")
ROOT.gROOT.ProcessLine(".L /afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/211210/FCCAnalyses/install/lib/libFCCAnalyses.C+")

parser = argparse.ArgumentParser()
parser.add_argument("-inputFiles", default = "/afs/cern.ch/user/m/mmlynari/workspace/ALLEGRO_PandoraPFA/HCal_standalone/outputs/241208/ALLEGRO_reco_pMin_100000_e.root", help = "Regex for input files.", type = str)
#parser.add_argument("-inputFiles", default = "/eos/user/m/mmlynari/FCC_rootfile_storage/HCal_v26May23/fcc_analysis_ouput/231220_energies_10kevt_cellsAndSW_noNoise_HCal_Sci_Steel_2x10_4x150_4x250_lowEne_SF_1/fccsw_output_pdgID_211_pMin_*_pMax_*_thetaMin_69.805_thetaMax_69.805.root", help = "Regex for input files.", type = str)
parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d") + "_v28Aug24_energies_10kevt_cells_noNoise_Barrel", help = "Postfix to append to the output folder.", type = str)
parser.add_argument("-color", default = 46, help = "Color of the graph", type = int)
parser.add_argument("-markerStyle", default = 21, help = "Style of the graph markers", type = int)
parser.add_argument("-cells", default = True, help = "Also produce plots with total energy from cells (no clustering) -- needs some refreshing", type = bool)
parser.add_argument("--specialLabel", help="Additional label to be plotted", type=str)
parser.add_argument("--materialLabel", help="Additional label to be plotted", type=str)
args = parser.parse_args()

plot_dir_name = 'plots_performances_'+ args.outputPostfix
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

individual_plots = os.path.join(plot_dir_name, "resolution") 
if not os.path.isdir(individual_plots):
    os.mkdir(individual_plots)

inputFiles = glob.glob(args.inputFiles)
if not inputFiles:
    print("No file found")
    sys.exit(1)

def get_response_and_resol(h, mean_guess=0, resol_guess=1):
    res = h.Fit("gaus", "LS", "", mean_guess-2*resol_guess, mean_guess+2*resol_guess)
    # Resolution should be corrected for the response (i.e as if response was brutally adjusted per energy)
    return (res.Parameter(1), res.Parameter(2)/(1+res.Parameter(1)))


# Single resolution histogram 

def draw_resol_canvas(th1, prefix, variable): 
    canvas_resol = ROOT.TCanvas(prefix + variable + "_resolution", prefix + variable + "_resolution")
    #fit_range_min = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin()) - 2 * th1.GetRMS()
    #fit_range_max = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin()) + 2 * th1.GetRMS()
    fit_range_min = th1.GetMean() - 2 * th1.GetRMS()
    fit_range_max = th1.GetMean() + 2 * th1.GetRMS()
    fit_result = th1.Fit("gaus", "SQ", "", fit_range_min, fit_range_max)
    if not fit_result.IsValid():
            raise ValueError("Fit result is invalid")
    th1.Draw()
    if 'rel' in variable:
        name = variable.replace("rel", "")
        th1.GetXaxis().SetTitle("%s_{Reco} - %s_{Gen}/%s_{Reco}"%(name, name, name))
        th1.GetXaxis().SetTitleOffset(1.2)
        if 'eGenDivided' in variable:
            th1.GetXaxis().SetTitle("E_{Reco} - E_{Gen}/E_{Gen}")
    elif 'esponse' in variable:
        th1.GetXaxis().SetTitle("#frac{E_{Reco}}{E_{Gen}}")
        th1.GetXaxis().SetTitleOffset(1.2)

    else:
        #th1.GetXaxis().SetTitle("#{0}_{Reco} - #{0}_{Gen}".format(variable))
        th1.GetXaxis().SetTitle("#%s_{Reco} - #%s_{Gen}"%(variable, variable))
    canvas_resol.Print(os.path.join(individual_plots, prefix + variable + "resol.png"))
    th1.Write()
    return fit_result

max_evt = 10000
invSF=1.

print('invSF ', invSF)
dict_energy_relEresol_error = {}
dict_energy_relEresol_eGenDivided_error = {}
dict_energy_energyResponseResol_error = {}
dict_energy_Ereco_error = {} 
dict_energy_linearity_error = {}
energies_gev_float = []
for inputFile in inputFiles:
    print("Treating %s..."%inputFile)
    rootfile_path = inputFile
    # lines below require a specific format of the input file name
    energy = inputFile.split('_pMin_')[1].split("_")[0]
    energy_gev_float = int(energy)/1000.0
    energies_gev_float.append(energy_gev_float)
    energy = str(energy_gev_float).replace(".","dot")
    print(energy)
    dict_energy_relEresol_error[energy] = []
    dict_energy_relEresol_eGenDivided_error[energy] = []
    dict_energy_energyResponseResol_error[energy] = []
    dict_energy_Ereco_error[energy] = []
    dict_energy_linearity_error[energy] = []

    prefix = energy + "GeV_"

    th1_Eresol = ROOT.TH1F(prefix + "energy_resolution", prefix + "energy_resolution", 500, -50, 50)
    th1_relEresol = ROOT.TH1F(prefix + "relative_energy_resolution", prefix + "relative_energy_resolution", 100, -0.6, 0.6)
    th1_relEresol_eGenDivided = ROOT.TH1F(prefix + "relative_energy_resolution_toEGen", prefix + "relative_energy_resolution_toEgen", 100, -0.6, 0.6)
    th1_energy_response = ROOT.TH1F(prefix + "energy_response", prefix + "energy_response", 100, 0, 2)
    #th1_Ereco = ROOT.TH1F(prefix + "reconstructed_energy", prefix + "reconstructed_energy", 500, energy_gev_float-0.5*energy_gev_float,energy_gev_float+0.5*energy_gev_float)
    th1_Ereco = ROOT.TH1F(prefix + "reconstructed_energy", prefix + "reconstructed_energy", 200, energy_gev_float-0.8*energy_gev_float,energy_gev_float+0.8*energy_gev_float)
    #th1_Ereco = ROOT.TH1F(prefix + "reconstructed_energy", prefix + "reconstructed_energy", 2000, 0, 2000)

    f = ROOT.TFile(rootfile_path)
    events = f.Get("events")
    events.SetBranchStatus("*", 0) # //disable all branches
    events.SetBranchStatus("genParticle_*", 1)
    events.SetBranchStatus("HCalBarrelReadoutPositioned_energy", 1)
    events.SetBranchStatus("HCalBarrelReadoutPositioned_theta", 1)
    events.SetBranchStatus("HCalBarrelReadoutPositioned_phi", 1)
    events.SetBranchStatus("HCalBarrelReadoutPositioned_layer", 1)

    evt = 0
    n_gen_particles = 0
    calo_label='HCal'

    for event in events:
        if evt >= max_evt and not max_evt == -1:
            break
        evt += 1
        total_energy = 0
        gen_particle_energy = 0
        for cell in range(len(getattr(event, "HCalBarrelReadoutPositioned_energy"))):
            cell_energy = event.HCalBarrelReadoutPositioned_energy[cell]*invSF
            total_energy += cell_energy

        #print('total_energy ', total_energy)
        for genParticle_idx in range(len(event.genParticle_status)):
            if event.genParticle_status[genParticle_idx] != 1:
                continue
            n_gen_particles += 1

            d_E = total_energy - event.genParticle_energy[genParticle_idx]
            d_relE = d_E/total_energy
            d_relE_eGenDivided = d_E/event.genParticle_energy[genParticle_idx]
            energyResponse =  total_energy / event.genParticle_energy[genParticle_idx]
            gen_particle_energy = event.genParticle_energy[genParticle_idx]

            th1_Ereco.Fill(total_energy)    
            th1_Eresol.Fill(d_E)
            th1_relEresol.Fill(d_relE_eGenDivided)
            th1_relEresol_eGenDivided.Fill(d_relE_eGenDivided)
            th1_energy_response.Fill(energyResponse)

        if evt % 1000 == 0:
            print("Event processed: %d"%evt)

    output_rootfile_path = os.path.join(individual_plots, prefix + "perfHistograms.root")
    output_rootfile = ROOT.TFile(output_rootfile_path, "recreate")
    output_rootfile.cd()

    ErecoFit = draw_resol_canvas(th1_Ereco, prefix, 'recoE')
    EresolFit = draw_resol_canvas(th1_Eresol, prefix, 'E')
    relEresolFit = draw_resol_canvas(th1_relEresol, prefix, 'relE')
    eResponseFit = draw_resol_canvas(th1_energy_response, prefix, 'energyResponse')
    

    relEresolFit_eGenDivided = draw_resol_canvas(th1_relEresol_eGenDivided, prefix, 'relE_eGenDivided')
    dict_energy_relEresol_eGenDivided_error[energy].append(relEresolFit_eGenDivided.Get().Parameter(2))
    dict_energy_relEresol_eGenDivided_error[energy].append(relEresolFit_eGenDivided.Get().ParError(2))
    
    ## alternative way to calculate resolution - same as what is done in BRT energy 
    #resp_e_v, resol_e_v = get_response_and_resol(th1_relEresol_eGenDivided, th1_relEresol_eGenDivided.GetMean(), th1_relEresol_eGenDivided.GetStdDev())
    #dict_energy_relEresol_eGenDivided_error[energy].append(resol_e_v)
    #dict_energy_linearity_error[energy].append(resp_e_v)
    #print('resp_e_v ', resp_e_v)

    dict_energy_Ereco_error[energy].append(ErecoFit.Get().Parameter(2)/ErecoFit.Get().Parameter(1))
    dict_energy_Ereco_error[energy].append(ErecoFit.Get().ParError(2))

    dict_energy_relEresol_error[energy].append(relEresolFit.Get().Parameter(2))
    dict_energy_relEresol_error[energy].append(relEresolFit.Get().ParError(2))

    dict_energy_energyResponseResol_error[energy].append(eResponseFit.Get().Parameter(2))
    dict_energy_energyResponseResol_error[energy].append(eResponseFit.Get().ParError(2))

    dict_energy_linearity_error[energy].append(ErecoFit.Get().Parameter(1)/gen_particle_energy)

    
    output_rootfile.Close()

postfix = ""

def plot_resolution_vs_energy_graph(variable_name, postfix, relEresol_vs_energy_graph, write_formula = False):
    setGridx = False
    setGridy = False
    if 'relEresol' in variable_name:
        plot_title = variable_name
        y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco} - E_{Gen}}{E_{Reco}}#right)"
        if 'eGenDivided' in variable_name:
            y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco} - E_{Gen}}{E_{Gen}}#right)"
    elif 'esponse' in variable_name:
        plot_title = variable_name
        y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco}}{E_{Gen}}#right)"
    elif 'esponse' in variable_name:
        plot_title = 'Energy response '
        y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco}}{E_{Gen}}#right)"
    elif 'Ereco' in variable_name:
        plot_title = variable_name
        y_axis_label = "#scale[1.4]{#sigma_{E}/#LTE#GT}"
        setGrid = False
    elif 'linearity' in variable_name:
        plot_title = variable_name
        #y_axis_label = "linearity"
        y_axis_label = "#scale[1.4]{#LTE_{rec}#GT/E_{true}}"
        setGridx = False
        setGridy = True
    print(y_axis_label)

    x_axis_label = "E_{Gen} [GeV]"
    #relEresol_vs_energy_graph.SetMarkerSize(1.5)
    relEresol_vs_energy_graph.SetMarkerStyle(args.markerStyle)
    relEresol_vs_energy_graph.SetMarkerColor(args.color)
    relEresol_vs_energy_graph.SetTitle(plot_title)
    relEresol_vs_energy_graph.SetName(plot_title)
    relEresol_vs_energy_graph.GetXaxis().SetTitle(x_axis_label)
    relEresol_vs_energy_graph.GetXaxis().SetTitleOffset(1.2)
    relEresol_vs_energy_graph.GetXaxis().SetLimits(2, 300)
    relEresol_vs_energy_graph.GetYaxis().SetTitle(y_axis_label)
    relEresol_vs_energy_graph.GetYaxis().SetTitleOffset(1.45)
    tf1_relEresol_vs_e = ROOT.TF1("tf1_" + variable_name, "sqrt(pow([0]/x, 2) + pow([1]/sqrt(x), 2) + pow([2], 2))", energies_gev_float[0], energies_gev_float[-1])
    #tf1_relEresol_vs_e = ROOT.TF1("tf1_" + variable_name, "sqrt([0]*[0] + pow([1]/sqrt(x),2))", energies_gev_float[0], energies_gev_float[-1])
    tf1_relEresol_vs_e.SetLineColor(args.color)
    if 'linearity' in variable_name:
        relEresol_vs_energy_graph.SetMinimum(0.5)
        relEresol_vs_energy_graph.SetMaximum(1.)
    else:
        relEresol_vs_energy_graph.SetMinimum(0)
    canvas_relEresol_vs_energy = ROOT.TCanvas(variable_name + "_vs_energy", variable_name + "_vs_energy")
    canvas_relEresol_vs_energy.SetLogx(1)
    canvas_relEresol_vs_energy.SetGrid(setGridx, setGridy)

    ## LEGEND
    if variable_name == "linearity":
        eResol_vs_e_legend = copy(gStyle.topRight_legend)
    else: 
        eResol_vs_e_legend = copy(gStyle.topRight_legend_relEresol)
    eResol_vs_e_legend.AddEntry(ROOT.nullptr,"#bf{FCC-ee simulation}","")
    eResol_vs_e_legend.AddEntry(ROOT.nullptr, "#pi^{#minus} @ #theta=69, B=0T", "")
    eResol_vs_e_legend.AddEntry(ROOT.nullptr,  "%s Barrel"%(calo_label), "")
    #eResol_vs_e_legend.AddEntry(ROOT.nullptr,  "%s Endcap, %s"%(calo_label,args.materialLabel), "")
    #eResol_vs_e_legend.AddEntry(ROOT.nullptr,  "%s"%(args.specialLabel), "")
    eResol_vs_e_legend.AddEntry(ROOT.nullptr, "cells, EM scale", "")

    if write_formula:
        relEresol_vs_e_fit = relEresol_vs_energy_graph.Fit(tf1_relEresol_vs_e, "SM")
        print(relEresol_vs_e_fit.Get().Parameter(0))
        a = str(round(relEresol_vs_e_fit.Get().Parameter(0)*100, 4))
        b = str(round(relEresol_vs_e_fit.Get().Parameter(1)*100, 2))
        c = str(round(relEresol_vs_e_fit.Get().Parameter(2)*100, 4))
        eResol_vs_e_formula = "#color[%d]{#frac{%s}{E} #oplus #frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, a, b, c) 
        #eResol_vs_e_formula = "#color[%d]{#frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, b, c) 
        #eResol_vs_e_legend = copy(gStyle.topRight_legend_relEresol)
        eResol_vs_e_legend.AddEntry(ROOT.nullptr, "", "")
        eResol_vs_e_legend.AddEntry(ROOT.nullptr, eResol_vs_e_formula, "")
    relEresol_vs_energy_graph.Draw("ap")
    eResol_vs_e_legend.Draw()
    #if 'Ereco'==variable_name:
    relEresol_vs_energy_graph.Write()
    canvas_relEresol_vs_energy.Print(os.path.join(plot_dir_name, variable_name + "_vs_energy.png"))
    canvas_relEresol_vs_energy.Write()


energies_gev_float.sort()
idx = 0
Ereco_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_Ereco_vs_energy_graph",prefix + postfix + "_Ereco_vs_energy_graph")
relEresol_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_relEresol_vs_energy_graph")
relEresol_eGenDivided_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_relEresol_eGenDivided_vs_energy_graph")
energyResponseResol_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_energyResponseResol_vs_energy_graph")
linearity_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_graph")

for energy_float in energies_gev_float:
    energy_str = str(energy_float).replace(".", "dot")

    # E resol
    Ereco_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_Ereco_error[energy_str][0])
    #Ereco_vs_energy_graph.SetPointError(idx, 0, dict_energy_Ereco_error[energy_str][1])
    relEresol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_relEresol_error[energy_str][0])
    relEresol_vs_energy_graph.SetPointError(idx, 0, dict_energy_relEresol_error[energy_str][1])
    relEresol_eGenDivided_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_relEresol_eGenDivided_error[energy_str][0])
    #relEresol_eGenDivided_vs_energy_graph.SetPointError(idx, energy_float, dict_energy_relEresol_eGenDivided_error[energy_str][1])
    energyResponseResol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_energyResponseResol_error[energy_str][0])
    energyResponseResol_vs_energy_graph.SetPointError(idx, 0, dict_energy_energyResponseResol_error[energy_str][1])
    linearity_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_linearity_error[energy_str][0])

    idx +=1

# Writing the VARIABLE_RESOL_VS_ENERGY
ROOT.gStyle.SetOptFit(000)

output_rootfile_graph_name = os.path.join(individual_plots, "relResol_vs_energy.root")
print("Writing %s"%output_rootfile_graph_name)
output_rootfile_graph = ROOT.TFile(output_rootfile_graph_name, "recreate")

Ereco_vs_energy_fit = plot_resolution_vs_energy_graph('Ereco', "", Ereco_vs_energy_graph, write_formula = True)
relEresol_vs_energy_fit = plot_resolution_vs_energy_graph('relEresol', "", relEresol_vs_energy_graph, write_formula = True)
relEresol_eGenDivided_vs_energy_fit = plot_resolution_vs_energy_graph('relEresol_eGenDivided', "", relEresol_eGenDivided_vs_energy_graph, write_formula = True)
energyResponseResol_vs_energy_fit = plot_resolution_vs_energy_graph('energyResponseResol', "", energyResponseResol_vs_energy_graph, write_formula = True)
linearity_vs_energy_fit = plot_resolution_vs_energy_graph('linearity', "", linearity_vs_energy_graph)



output_rootfile_graph.Close()
