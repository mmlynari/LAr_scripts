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
parser.add_argument("-inputFiles", default = "/eos/user/m/mmlynari/FCC_fellow/FCC_rootfile_storage/HCal_v28Aug24_FSR/fcc_analysis_ouput/240828_energies_10kevt_SW_noNoise_HCal/fccsw_output_pdgID_211_pMin_*.root", help = "Regex for input files.", type = str)
#parser.add_argument("-inputFiles", default = "/eos/user/b/brfranco/rootfile_storage/fcc_analysis/220308_energies_10kevt_topoAndSW_noise_splitLowEnergy30_caloReco/fccsw_output_pdgID_22_pMin_*.root", help = "Regex for input files.", type = str)
parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d") + "_CaloClusters", help = "Postfix to append to the output folder.", type = str)
parser.add_argument("-color", default = 46, help = "Color of the graph", type = int)
parser.add_argument("-markerStyle", default = 21, help = "Style of the graph markers", type = int)
parser.add_argument("-cells", default = True, help = "Also produce plots with total energy from cells (no clustering) -- needs some refreshing", type = bool)
parser.add_argument("-clusterCollection", default = "CaloClusters", help = "Name of the cluster collection to use", type = str)
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

## sfor the plot legend
scale = "EM scale"
# scale = "HAD scale"
calo_label='HCal'


# Single resolution histogram 

def draw_resol_canvas(th1, prefix, variable): 
    canvas_resol = ROOT.TCanvas(prefix + variable + "_resolution", prefix + variable + "_resolution")
    fit_range_min = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin()) - 2 * th1.GetRMS()
    fit_range_max = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin()) + 2 * th1.GetRMS()
    fit_result = th1.Fit("gaus", "SQ", "", fit_range_min, fit_range_max)
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
#cutoff_dR = 0.015
#cutoff_relE = 0.2
cutoff_dR = 0.03
cutoff_dPhi = 0.1
cutoff_dTheta = 0.1
cutoff_relE = 0.5
dict_energy_relEresol_error = {}
dict_energy_relEresol_eGenDivided_error = {}
dict_energy_energyResponseResol_error = {}
dict_energy_efficiency_error = {}
dict_energy_phiResol_error = {}
dict_energy_thetaResol_error = {}
dict_energy_cells_relEresol_error = {}
dict_energy_cells_Ereco_error = {} 
dict_energy_Ereco_error = {} 
dict_energy_linearity_error = {}
dict_energy_cells_linearity_error = {}


energies_gev_float = []
for inputFile in inputFiles:
    print("Treating %s..."%inputFile)
    rootfile_path = inputFile
    energy = inputFile.split('_pMin_')[1].split("_")[0]
    energy_gev_float = int(energy)/1000.0
    energies_gev_float.append(energy_gev_float)
    energy = str(energy_gev_float).replace(".","dot")
    print(energy)
    dict_energy_relEresol_error[energy] = []
    dict_energy_relEresol_eGenDivided_error[energy] = []
    dict_energy_energyResponseResol_error[energy] = []
    dict_energy_phiResol_error[energy] = []
    dict_energy_thetaResol_error[energy] = []
    dict_energy_cells_relEresol_error[energy] = []
    dict_energy_efficiency_error[energy] = []
    dict_energy_Ereco_error[energy] = []
    dict_energy_cells_Ereco_error[energy] = []
    dict_energy_linearity_error[energy] = []
    dict_energy_cells_linearity_error[energy] = []
    prefix = energy + "GeV_"

    phiresol_range = 0.002
    if energy_gev_float < 6:
        phiresol_range = 0.01
    if energy_gev_float < 6:
        phiresol_range = 0.03

    th1_phiresol = ROOT.TH1F(prefix + "phi_resolution", prefix + "phi_resolution", 100, -1 * phiresol_range, phiresol_range)
    th1_thetaresol = ROOT.TH1F(prefix + "theta_resolution", prefix + "theta_resolution", 100, -0.005, 0.005)
    th1_angularresol = ROOT.TH1F(prefix + "angular_resolution", prefix + "angular_resolution", 100, 0, 0.03)
    th1_Eresol = ROOT.TH1F(prefix + "energy_resolution", prefix + "energy_resolution", 500, -25, 20)
    th1_Eresol_cells = ROOT.TH1F(prefix + "energy_resolution_cells", prefix + "energy_resolution_cells", 500, -25, 20)
    th1_relEresol = ROOT.TH1F(prefix + "relative_energy_resolution", prefix + "relative_energy_resolution", 100, -0.6, 0.6)
    th1_relEresol_eGenDivided = ROOT.TH1F(prefix + "relative_energy_resolution_toEGen", prefix + "relative_energy_resolution_toEgen", 150, -0.6, 0.9)
    th1_energy_response = ROOT.TH1F(prefix + "energy_response", prefix + "energy_response", 100, 0, 2)
    th1_relEresol_cells = ROOT.TH1F(prefix + "relative_energy_resolution_cells", prefix + "relative_energy_resolution_cells", 100, -0.6, 0.6)
    th1_Ereco = ROOT.TH1F(prefix + "reconstructed_energy", prefix + "reconstructed_energy", 200, energy_gev_float-0.8*energy_gev_float,energy_gev_float+0.8*energy_gev_float)
    th1_Ereco_cells = ROOT.TH1F(prefix + "reconstructed_energy_cells", prefix + "reconstructed_energy_cells", 200, energy_gev_float-0.8*energy_gev_float,energy_gev_float+0.8*energy_gev_float)

    f = ROOT.TFile(rootfile_path)
    events = f.Get("events")
    events.SetBranchStatus("*", 0) # //disable all branches
    events.SetBranchStatus(args.clusterCollection + "_*", 1)
    events.SetBranchStatus("genParticle_*", 1)
    if args.cells:
        events.SetBranchStatus("HCalBarrelPositionedCells2_energy", 1)


    theta_cells = []
    evt = 0
    n_gen_particles = 0
    evt_with_cluster_matching_genParticle = 0

    for event in events:
        if evt >= max_evt and not max_evt == -1:
            break
        evt += 1
        total_energy_cells = 0
        gen_particle_energy = 0
        if args.cells:
            for cell_energy in event.HCalBarrelPositionedCells2_energy:
                total_energy_cells += cell_energy

        for genParticle_idx in range(len(event.genParticle_status)):
            if event.genParticle_status[genParticle_idx] != 1:
                continue
            n_gen_particles += 1
            best_genParticle_match_idx = -1
            best_d_R = 100000
            best_d_theta = 100000
            best_d_phi = 100000
            best_d_E = 100000
            best_d_relE = 10000
            best_d_relE_eGenDivided = 10000
            best_energyResponse = 10000
            for caloCluster_idx in range(len(getattr(event, args.clusterCollection + "_energy"))):
                d_theta = getattr(event, args.clusterCollection + "_theta")[caloCluster_idx] - event.genParticle_theta[genParticle_idx]
                d_phi = getattr(event, args.clusterCollection + "_phi")[caloCluster_idx] - event.genParticle_phi[genParticle_idx]
                d_R = sqrt(d_theta * d_theta + d_phi * d_phi)
                d_E = getattr(event, args.clusterCollection + "_energy")[caloCluster_idx] - event.genParticle_energy[genParticle_idx]
                d_relE = d_E/getattr(event, args.clusterCollection + "_energy")[caloCluster_idx]
                d_relE_eGenDivided = d_E/event.genParticle_energy[genParticle_idx]
                energy_response =  getattr(event, args.clusterCollection + "_energy")[caloCluster_idx] / event.genParticle_energy[genParticle_idx]
                E_reco_cluster = getattr(event, args.clusterCollection + "_energy")[caloCluster_idx]
                #d_relE = (event.CorrectedCaloClusters_energy[caloCluster_idx] - event.genParticle_energy[genParticle_idx])/event.genParticle_energy[genParticle_idx]

                if d_R < best_d_R:
                    best_genParticle_match_idx = caloCluster_idx
                    best_d_R = d_R
                    best_d_theta = d_theta
                    best_d_phi = d_phi
                    best_d_E = d_E
                    best_d_relE = d_relE
                    best_d_relE_eGenDivided = d_relE_eGenDivided
                    best_energyResponse = energy_response
                    gen_particle_energy = event.genParticle_energy[genParticle_idx]

            # for energy resol, ask for phi and theta closeness and vice versa
            if best_genParticle_match_idx == -1:
                continue
            if abs(best_d_relE) < cutoff_relE:
                th1_phiresol.Fill(best_d_phi)
                th1_thetaresol.Fill(best_d_theta)
                th1_angularresol.Fill(best_d_R)
                
            if best_d_R < cutoff_dR:
                th1_Eresol.Fill(best_d_E)
                th1_relEresol.Fill(best_d_relE)
                th1_relEresol_eGenDivided.Fill(best_d_relE_eGenDivided)
                th1_energy_response.Fill(best_energyResponse)
                th1_Ereco.Fill(E_reco_cluster)    

                if args.cells:
                    th1_Eresol_cells.Fill((total_energy_cells - gen_particle_energy))
                    th1_relEresol_cells.Fill((total_energy_cells - gen_particle_energy)/total_energy_cells)
                    th1_Ereco_cells.Fill(total_energy_cells)  
            # Ask both for dR and dE matching for the efficiency
            if best_d_R < cutoff_dR and abs(best_d_relE) < cutoff_relE:
                evt_with_cluster_matching_genParticle += 1


        if evt % 1000 == 0:
            print("Event processed: %d"%evt)

    output_rootfile_path = os.path.join(individual_plots, prefix + "perfHistograms.root")
    output_rootfile = ROOT.TFile(output_rootfile_path, "recreate")
    output_rootfile.cd()

    phiResolFit = draw_resol_canvas(th1_phiresol, prefix, 'Phi')
    thetaResolFit = draw_resol_canvas(th1_thetaresol, prefix, 'Theta')
    EresolFit = draw_resol_canvas(th1_Eresol, prefix, 'E')
    relEresolFit = draw_resol_canvas(th1_relEresol, prefix, 'relE')
    relEresolFit_eGenDivided = draw_resol_canvas(th1_relEresol_eGenDivided, prefix, 'relE_eGenDivided')
    eResponseFit = draw_resol_canvas(th1_energy_response, prefix, 'energyResponse')
    ErecoFit = draw_resol_canvas(th1_Ereco, prefix, 'recoE')

    dict_energy_efficiency_error[energy].append(evt_with_cluster_matching_genParticle * 100 / float(n_gen_particles))
    dict_energy_efficiency_error[energy].append(sqrt(evt_with_cluster_matching_genParticle) * 100 / float(n_gen_particles))

    dict_energy_relEresol_error[energy].append(relEresolFit.Get().Parameter(2))
    dict_energy_relEresol_error[energy].append(relEresolFit.Get().ParError(2))

    dict_energy_relEresol_eGenDivided_error[energy].append(relEresolFit_eGenDivided.Get().Parameter(2))
    dict_energy_relEresol_eGenDivided_error[energy].append(relEresolFit_eGenDivided.Get().ParError(2))

    dict_energy_energyResponseResol_error[energy].append(relEresolFit.Get().Parameter(2))
    dict_energy_energyResponseResol_error[energy].append(relEresolFit.Get().ParError(2))

    dict_energy_phiResol_error[energy].append(phiResolFit.Get().Parameter(2))
    dict_energy_phiResol_error[energy].append(phiResolFit.Get().ParError(2))

    dict_energy_thetaResol_error[energy].append(thetaResolFit.Get().Parameter(2))
    dict_energy_thetaResol_error[energy].append(thetaResolFit.Get().ParError(2))

    dict_energy_Ereco_error[energy].append(ErecoFit.Get().Parameter(2)/ErecoFit.Get().Parameter(1))
    dict_energy_Ereco_error[energy].append(ErecoFit.Get().ParError(2))

    dict_energy_linearity_error[energy].append(ErecoFit.Get().Parameter(1)/gen_particle_energy)

    canvas_angularresol = ROOT.TCanvas(prefix + "angular_resolution", prefix + "angular_resolution")
    th1_angularresol.Draw()
    th1_angularresol.GetXaxis().SetTitle("#Delta R = #sqrt{#Delta #Phi^{2} + #Delta #Theta^{2}}")
    th1_angularresol.GetXaxis().SetTitleOffset(1.2)
    canvas_angularresol.Print(os.path.join(individual_plots, prefix + "angularresol.png"))
    th1_angularresol.Write()

    if args.cells:
        ErecoFit_cells = draw_resol_canvas(th1_Ereco_cells, prefix, 'recoEcell')
        dict_energy_cells_Ereco_error[energy].append(ErecoFit_cells.Get().Parameter(2)/ErecoFit_cells.Get().Parameter(1))
        dict_energy_cells_Ereco_error[energy].append(ErecoFit_cells.Get().ParError(2))
        dict_energy_cells_linearity_error[energy].append(ErecoFit_cells.Get().Parameter(1)/gen_particle_energy)

        relEresolFit_cells = draw_resol_canvas(th1_relEresol_cells, prefix, 'relEcell')
        dict_energy_cells_relEresol_error[energy].append(relEresolFit_cells.Get().Parameter(2))
        dict_energy_cells_relEresol_error[energy].append(relEresolFit_cells.Get().ParError(2))
        #canvas_relEresol_cells = ROOT.TCanvas(prefix + "relE_resolution_cells", prefix + "relE_resolution_cells")
        #fit_range_min_cells = th1_relEresol_cells.GetXaxis().GetBinCenter(th1_relEresol_cells.GetMaximumBin()) - 2 * th1_relEresol_cells.GetRMS()
        #fit_range_max_cells = th1_relEresol_cells.GetXaxis().GetBinCenter(th1_relEresol_cells.GetMaximumBin()) + 2 * th1_relEresol_cells.GetRMS()
        ##fit_range_max = th1_relEresol_cells.GetMean() + 1 * th1_relEresol_cells.GetRMS()
        #relEresolFit_cells = th1_relEresol_cells.Fit("gaus", "SQ", "", fit_range_min_cells, fit_range_max_cells)
        ##relEresolFit_cells = th1_relEresol_cells.Fit("gaus", "SQ", "")
        #th1_relEresol_cells.GetXaxis().SetTitle("(E_{Reco} - E_{Gen})/E_{Reco}")
        #th1_relEresol_cells.GetXaxis().SetTitleOffset(1.2)
        #th1_relEresol_cells.Draw()
        #canvas_relEresol_cells.Print(os.path.join(plot_dir_name, prefix + "relEresol_cells.png"))
        #th1_relEresol_cells.Write()
        #dict_energy_cells_relEresol_error[energy].append(relEresolFit_cells.Get().Parameter(2))
        #dict_energy_cells_relEresol_error[energy].append(relEresolFit_cells.Get().ParError(2))

    output_rootfile.Close()

postfix = ""

def plot_resolution_vs_energy_graph(variable_name, postfix, relEresol_vs_energy_graph, write_formula = False):
    setGrid = False
    if 'relEresol' in variable_name:
        plot_title = 'HCAL energy resolution '
        y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco} - E_{Gen}}{E_{Reco}}#right)"
        if 'eGenDivided' in variable_name:
            y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco} - E_{Gen}}{E_{Gen}}#right)"
    elif 'esponse' in variable_name:
        plot_title = 'HCAL energy response '
        y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco}}{E_{Gen}}#right)"
    elif variable_name == 'phiResol':
        plot_title = 'HCAL #Phi resolution '
        y_axis_label = "#scale[1.9]{#sigma}#left(#Phi_{Reco} - #Phi_{Gen}#right) [mrad]"
    elif variable_name == 'thetaResol':
        plot_title = 'HCAL #theta resolution '
        y_axis_label = "#scale[1.9]{#sigma}#left(#theta_{Reco} - #theta_{Gen}#right) [mrad]"
    elif variable_name == "efficiency":
        plot_title = 'HCAL efficiency '
        y_axis_label = "Efficiency [%]"
        setGrid = True
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
    relEresol_vs_energy_graph.GetXaxis().SetTitle(x_axis_label)
    relEresol_vs_energy_graph.GetXaxis().SetTitleOffset(1.2)
    relEresol_vs_energy_graph.GetXaxis().SetLimits(1., 300)
    #relEresol_vs_energy_graph.GetXaxis().SetLimits(0.7, 300)
    relEresol_vs_energy_graph.GetYaxis().SetTitle(y_axis_label)
    relEresol_vs_energy_graph.GetYaxis().SetTitleOffset(1.45)
    #relEresol_vs_energy_graph = ROOT.TGraphErrors(len(energies_gev_float), energies_gev_float, energy_resolutions, 0, energy_resolution_errors)
    #relEresol_vs_energy_graph.GetYaxis().SetLabelSize(0.06)
    #relEresol_vs_energy_graph.GetYaxis().SetTitleSize(0.06)
    #relEresol_vs_energy_graph.GetYaxis().SetTitleSize(0.05)
    tf1_relEresol_vs_e = ROOT.TF1("tf1_" + variable_name, "sqrt(pow([0]/x, 2) + pow([1]/sqrt(x), 2) + pow([2], 2))", energies_gev_float[0], energies_gev_float[-1])
    tf1_relEresol_vs_e.SetLineColor(args.color)
    if 'linearity' in variable_name:
        relEresol_vs_energy_graph.SetMinimum(0.5)
        relEresol_vs_energy_graph.SetMaximum(1.)
    else:
        relEresol_vs_energy_graph.SetMinimum(0)
    canvas_relEresol_vs_energy = ROOT.TCanvas(variable_name + "_vs_energy", variable_name + "_vs_energy")
    canvas_relEresol_vs_energy.SetLogx(1)
    canvas_relEresol_vs_energy.SetGrid(setGrid, setGrid)
    if write_formula:
        relEresol_vs_e_fit = relEresol_vs_energy_graph.Fit(tf1_relEresol_vs_e, "MS")
        #print(relEresol_vs_e_fit.Get().Parameter(0))
        #print(relEresol_vs_e_fit.Get().Parameter(1))
        #print(relEresol_vs_e_fit.Get().Parameter(2))
        a = str(round(relEresol_vs_e_fit.Get().Parameter(0), 4))
        b = str(round(relEresol_vs_e_fit.Get().Parameter(1), 2))
        c = str(round(relEresol_vs_e_fit.Get().Parameter(2), 4))
        eResol_vs_e_formula = "#color[%d]{#frac{%s}{E} #oplus #frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, a, b, c) 
        #eResol_vs_e_formula = "#color[%d]{#frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, b, c) 
        eResol_vs_e_legend = copy(gStyle.topRight_legend_relEresol)
        eResol_vs_e_legend.AddEntry(ROOT.nullptr, eResol_vs_e_formula, "")
    relEresol_vs_energy_graph.Draw("ap")
    if write_formula:
        eResol_vs_e_legend.Draw()
    relEresol_vs_energy_graph.Write()
    canvas_relEresol_vs_energy.Print(os.path.join(plot_dir_name, variable_name + "_vs_energy.png"))
    canvas_relEresol_vs_energy.Write()
    #return relEresol_vs_e_fit

    ## LEGEND
    if 'linearity' in variable_name:
        eResol_vs_e_legend = copy(gStyle.bottomRight_legend)
    else: 
        eResol_vs_e_legend = copy(gStyle.topRight_legend_relEresol)
    eResol_vs_e_legend.AddEntry(ROOT.nullptr,"#bf{FCC-ee simulation}","")
    eResol_vs_e_legend.AddEntry(ROOT.nullptr, "#pi^{#minus} @ #theta=69, B=0T", "")
    eResol_vs_e_legend.AddEntry(ROOT.nullptr,  "%s Barrel"%(calo_label), "")
    #eResol_vs_e_legend.AddEntry(ROOT.nullptr,  "%s Endcap, %s"%(calo_label,args.materialLabel), "")
    #eResol_vs_e_legend.AddEntry(ROOT.nullptr,  "%s"%(args.specialLabel), "")
    #eResol_vs_e_legend.AddEntry(ROOT.nullptr, "cells, EM scale", "")
    if 'cell' in variable_name:
        eResol_vs_e_legend.AddEntry(ROOT.nullptr, "cells, %s"%(scale), "")
    else:
        eResol_vs_e_legend.AddEntry(ROOT.nullptr, "SW clusters, %s"%(scale), "")

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
relEresol_vs_energy_cells_graph = ROOT.TGraphErrors(prefix + postfix)
Ereco_vs_energy_cells_graph = ROOT.TGraphErrors(prefix + postfix + "_Ereco_vs_energy_cells_graph",prefix + postfix + "_Ereco_vs_energy_cells_graph")

#relEresol_vs_energy_cells_graph, relEerol_vs_energy_fit = create_resolution_vs_energy_graph('relEresol', "cells")

efficiency_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_efficiency_vs_energy_graph")
linearity_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_graph")
linearity_vs_energy_cells_graph = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_cells_graph")


phiResol_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_phiResol_vs_energy_graph")
thetaResol_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_thetaResol_vs_energy_graph")

for energy_float in energies_gev_float:
    energy_str = str(energy_float).replace(".", "dot")

    # E resol
    Ereco_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_Ereco_error[energy_str][0])
    relEresol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_relEresol_error[energy_str][0])
    relEresol_vs_energy_graph.SetPointError(idx, 0, dict_energy_relEresol_error[energy_str][1])
    relEresol_eGenDivided_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_relEresol_eGenDivided_error[energy_str][0])
    relEresol_eGenDivided_vs_energy_graph.SetPointError(idx, 0, dict_energy_relEresol_eGenDivided_error[energy_str][1])
    if args.cells:
        Ereco_vs_energy_cells_graph.SetPoint(idx, energy_float, dict_energy_cells_Ereco_error[energy_str][0])
        relEresol_vs_energy_cells_graph.SetPoint(idx, energy_float, dict_energy_cells_relEresol_error[energy_str][0])
        relEresol_vs_energy_cells_graph.SetPointError(idx, 0, dict_energy_cells_relEresol_error[energy_str][1])
    energyResponseResol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_energyResponseResol_error[energy_str][0])
    energyResponseResol_vs_energy_graph.SetPointError(idx, 0, dict_energy_energyResponseResol_error[energy_str][1])

    # Phi resol
    phiResol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_phiResol_error[energy_str][0]*1000)
    phiResol_vs_energy_graph.SetPointError(idx, 0, dict_energy_phiResol_error[energy_str][1]*1000)

    # Theta resol
    thetaResol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_thetaResol_error[energy_str][0]*1000)
    thetaResol_vs_energy_graph.SetPointError(idx, 0, dict_energy_thetaResol_error[energy_str][1]*1000)

    #Efficiency
    efficiency_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_efficiency_error[energy_str][0])
    efficiency_vs_energy_graph.SetPointError(idx, 0, dict_energy_efficiency_error[energy_str][1])

    linearity_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_linearity_error[energy_str][0])
    linearity_vs_energy_cells_graph.SetPoint(idx, energy_float, dict_energy_cells_linearity_error[energy_str][0])

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
if args.cells:
    relEresol_cell_vs_energy_fit = plot_resolution_vs_energy_graph('relEresolcells', "", relEresol_vs_energy_cells_graph, write_formula = True)
    Ereco_cells_vs_energy_fit = plot_resolution_vs_energy_graph('Erecocells', "", Ereco_vs_energy_cells_graph, write_formula = True)

phiResol_vs_energy_fit = plot_resolution_vs_energy_graph('phiResol', "", phiResol_vs_energy_graph)
thetaResol_vs_energy_fit = plot_resolution_vs_energy_graph('thetaResol', "", thetaResol_vs_energy_graph)
efficiency_vs_energy_fit = plot_resolution_vs_energy_graph('efficiency', "", efficiency_vs_energy_graph)
linearity_vs_energy_fit = plot_resolution_vs_energy_graph('linearity', "", linearity_vs_energy_graph)
linearity_cells_vs_energy_fit = plot_resolution_vs_energy_graph('linearitycells', "", linearity_vs_energy_cells_graph)



#Fit energy resolution vs energy
#tf1_eResol_vs_e = TF1("tf1_eResol_vs_e", "sqrt([0] * [0] + pow([1] / sqrt(x), 2))", 5, 600)

#tf1_eResol_vs_e = ROOT.TF1("tf1_eResol_vs_e", "sqrt(pow([0]/x, 2) + pow([1]/sqrt(x), 2) + pow([2], 2))", energies_gev_float[0], energies_gev_float[-1])
#tf1_eResol_vs_e.SetLineColor(args.color)
#eResol_vs_e_fit = relEresol_vs_energy_graph.Fit(tf1_eResol_vs_e, "SQ")
#a = str(round(eResol_vs_e_fit.Get().Parameter(0), 2))
#b = str(round(eResol_vs_e_fit.Get().Parameter(1), 2))
#c = str(round(eResol_vs_e_fit.Get().Parameter(2), 2))
#
#tf1_phiResol_vs_e = ROOT.TF1("tf1_phiResol_vs_e", "sqrt(pow([0]/x, 2) + pow([1]/sqrt(x), 2) + pow([2], 2))", energies_gev_float[0], energies_gev_float[-1])
#tf1_phiResol_vs_e.SetLineColor(args.color)
#phiResol_vs_e_fit = relEresol_vs_energy_graph.Fit(tf1_phiResol_vs_e, "SQ")
#a = str(round(phiResol_vs_e_fit.Get().Parameter(0), 2))
#b = str(round(phiResol_vs_e_fit.Get().Parameter(1), 2))
#c = str(round(phiResol_vs_e_fit.Get().Parameter(2), 2))

#FIXME
#canvas_relEresol_vs_energy = ROOT.TCanvas("relEresolution_vs_energy", "relEresolution_vs_energy")
#canvas_relEresol_vs_energy.SetLogx(1)
#relEresol_vs_energy_graph.Draw("ap")
#relEresol_vs_energy_graph.Write()
#eResol_vs_e_formula = "#color[%d]{#frac{%s}{E} #oplus #frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, a, b, c) 
#eResol_vs_e_legend = gStyle.topRight_legend_relEresol
#eResol_vs_e_legend.AddEntry(ROOT.nullptr, eResol_vs_e_formula, "")
#eResol_vs_e_legend.Draw()
#canvas_relEresol_vs_energy.Print(os.path.join(plot_dir_name, "relEresolution_vs_energy.png"))
#canvas_relEresol_vs_energy.Write()

#if args.cells:
#    eResol_vs_e_fit = relEresol_vs_energy_cells_graph.Fit(tf1_eResol_vs_e, "SQ")
#    a = str(round(eResol_vs_e_fit.Get().Parameter(0), 2))
#    b = str(round(eResol_vs_e_fit.Get().Parameter(1), 2))
#    c = str(round(eResol_vs_e_fit.Get().Parameter(2), 2))
#
#    canvas_relEresol_vs_energy_cells = ROOT.TCanvas("relEresolution_vs_energy_cells", "relEresolution_vs_energy_cells")
#    canvas_relEresol_vs_energy_cells.SetLogx(1)
#    relEresol_vs_energy_cells_graph.Draw("ap")
#    relEresol_vs_energy_cells_graph.Write()
#    eResol_vs_e_formula = "#color[%d]{#frac{%s}{E} #oplus #frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, a, b, c) 
#    gStyle.topRight_legend_relEresol.Clear()
#    eResol_vs_e_legend_cells = gStyle.topRight_legend_relEresol
#    eResol_vs_e_legend_cells.AddEntry(ROOT.nullptr, eResol_vs_e_formula, "")
#    eResol_vs_e_legend_cells.Draw()
#    canvas_relEresol_vs_energy_cells.Print(os.path.join(plot_dir_name, "relEresolution_vs_energy_cells.png"))
#    canvas_relEresol_vs_energy_cells.Write()

output_rootfile_graph.Close()
