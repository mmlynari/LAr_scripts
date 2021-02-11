import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
import os, glob
import numpy as np
import argparse
from math import sqrt
from datetime import date

import gStyle

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")

parser = argparse.ArgumentParser()
#parser.add_argument("-inputFiles", default = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/eGun1GeV_10kevt_originalGeometry_allCells/output_fullCalo_SimAndDigi_withCluster_noMagneticField_1GeV_pythiaFalse.root", help = "Name of the input file.", type = str)
parser.add_argument("-inputFiles", default = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/210130_caloReco/fccsw_output_pdgID_22_pMin_*_pMax_*_thetaMin_90_thetaMax_90.root", help = "Name of the input file.", type = str)
parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d") + "_sqrtSum", help = "Postfix to append to the output folder.", type = str)
parser.add_argument("-color", default = 46, help = "Color of the graph", type = int)
parser.add_argument("-markerStyle", default = 34, help = "Style of the graph markers", type = int)
parser.add_argument("-cells", default = True, help = "Also produce plots with total energy from cells (no clustering)", type = bool)
args = parser.parse_args()

plot_dir_name = 'plots_performances_'+ args.outputPostfix
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

inputFiles = glob.glob(args.inputFiles)

max_evt = -1
#cutoff_dR = 0.015
#cutoff_relE = 0.2
cutoff_dR = 1000
cutoff_relE = 1000
dict_energy_relEresol_error = {}
dict_energy_cells_relEresol_error = {}
energies_gev_float = []
for inputFile in inputFiles:
    print "Treating %s..."%inputFile
    rootfile_path = inputFile
    energy = inputFile.split('_pMin_')[1].split("_")[0]
    energy = int(energy)/1000.0
    energies_gev_float.append(energy)
    energy = str(energy).replace(".","dot")
    print energy
    dict_energy_relEresol_error[energy] = []
    dict_energy_cells_relEresol_error[energy] = []
    prefix = energy + "GeV_"

    th1_phiresol = ROOT.TH1F(prefix + "phi_resolution", prefix + "phi_resolution", 100, -0.03, 0.03)
    th1_thetaresol = ROOT.TH1F(prefix + "theta_resolution", prefix + "theta_resolution", 100, -0.03, 0.03)
    th1_angularresol = ROOT.TH1F(prefix + "angular_resolution", prefix + "angular_resolution", 100, 0, 0.03)
    th1_Eresol = ROOT.TH1F(prefix + "energy_resolution", prefix + "energy_resolution", 500, -25, 20)
    th1_Eresol_cells = ROOT.TH1F(prefix + "energy_resolution_cells", prefix + "energy_resolution_cells", 500, -25, 20)
    th1_relEresol = ROOT.TH1F(prefix + "relative_energy_resolution", prefix + "relative_energy_resolution", 100, -0.6, 0.6)
    th1_relEresol_cells = ROOT.TH1F(prefix + "relative_energy_resolution_cells", prefix + "relative_energy_resolution_cells", 100, -0.6, 0.6)

    f = ROOT.TFile(rootfile_path)
    events = f.Get("events")

    theta_cells = []
    evt = 0
    for event in events:
        if evt >= max_evt and not max_evt == -1:
            break
        total_energy = 0
        gen_particle_energy = 0
        if args.cells:
            for cell_energy in event.ECalBarrelPositionedCells_energy:
                total_energy += cell_energy

        for genParticle_idx in range(len(event.genParticle_status)):
            if event.genParticle_status[genParticle_idx] != 1:
                continue
            best_genParticle_match_idx = -1
            best_d_R = 100000
            best_d_theta = 100000
            best_d_phi = 100000
            best_d_E = 100000
            best_d_relE = 10000
            for caloCluster_idx in range(len(event.CaloClusters_energy)):
                d_theta = event.CaloClusters_theta[caloCluster_idx] - event.genParticle_theta[genParticle_idx]
                d_phi = event.CaloClusters_phi[caloCluster_idx] - event.genParticle_phi[genParticle_idx]
                d_R = sqrt(d_theta * d_theta + d_phi * d_phi)
                d_E = event.CaloClusters_energy[caloCluster_idx] - event.genParticle_energy[genParticle_idx]
                d_relE = (event.CaloClusters_energy[caloCluster_idx] - event.genParticle_energy[genParticle_idx])/event.CaloClusters_energy[caloCluster_idx]
                #d_relE = (event.CaloClusters_energy[caloCluster_idx] - event.genParticle_energy[genParticle_idx])/event.genParticle_energy[genParticle_idx]

                if d_R < best_d_R:
                    best_genParticle_match_idx = caloCluster_idx
                    best_d_R = d_R
                    best_d_theta = d_theta
                    best_d_phi = d_phi
                    best_d_E = d_E
                    best_d_relE = d_relE
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
                if args.cells:
                    th1_Eresol_cells.Fill((total_energy - gen_particle_energy))
                    th1_relEresol_cells.Fill((total_energy - gen_particle_energy)/total_energy)
                    #th1_relEresol_cells.Fill((total_energy - gen_particle_energy)/gen_particle_energy)

        evt += 1
        if evt % 1000 == 0:
            print "Event processed: %d"%evt

    output_rootfile_path = os.path.join(plot_dir_name, prefix + "perfHistograms.root")
    output_rootfile = ROOT.TFile(output_rootfile_path, "recreate")
    output_rootfile.cd()

    canvas_phiresol = ROOT.TCanvas(prefix + "phi_resolution", prefix + "phi_resolution")
    th1_phiresol.Fit("gaus", "Q")
    th1_phiresol.Draw()
    th1_phiresol.GetXaxis().SetTitle("#Delta #Phi")
    canvas_phiresol.Print(os.path.join(plot_dir_name, prefix + "phiresol.png"))
    th1_phiresol.Write()

    canvas_thetaresol = ROOT.TCanvas(prefix + "theta_resolution", prefix + "theta_resolution")
    th1_thetaresol.Fit("gaus", "Q")
    th1_thetaresol.Draw()
    th1_thetaresol.GetXaxis().SetTitle("#Delta #Theta")
    canvas_thetaresol.Print(os.path.join(plot_dir_name, prefix + "thetaresol.png"))
    th1_thetaresol.Write()

    canvas_angularresol = ROOT.TCanvas(prefix + "angular_resolution", prefix + "angular_resolution")
    th1_angularresol.Draw()
    th1_angularresol.GetXaxis().SetTitle("#Delta R = #sqrt{#Delta #Phi^{2} + #Delta #Theta^{2}}")
    th1_angularresol.GetXaxis().SetTitleOffset(1.2)
    canvas_angularresol.Print(os.path.join(plot_dir_name, prefix + "angularresol.png"))
    th1_angularresol.Write()

    canvas_Eresol = ROOT.TCanvas(prefix + "E_resolution", prefix + "E_resolution")
    th1_Eresol.Fit("gaus", "Q")
    th1_Eresol.Draw()
    th1_Eresol.GetXaxis().SetTitle("E_{Reco} - E_{Gen}")
    canvas_Eresol.Print(os.path.join(plot_dir_name, prefix + "Eresol.png"))
    th1_Eresol.Write()

    canvas_relEresol = ROOT.TCanvas(prefix + "relE_resolution", prefix + "relE_resolution")
    fit_range_min = th1_relEresol.GetMean() - 1 * th1_relEresol.GetRMS()
    fit_range_max = th1_relEresol.GetMean() + 1 * th1_relEresol.GetRMS()
    #th1_relEresol.Fit("gaus", "Q", "", fit_range_min, fit_range_max)
    relEresolFit = th1_relEresol.Fit("gaus", "SQ", "")
    th1_relEresol.GetXaxis().SetTitle("(E_{Reco} - E_{Gen})/E_{Reco}")
    th1_relEresol.GetXaxis().SetTitleOffset(1.2)
    th1_relEresol.Draw()
    canvas_relEresol.Print(os.path.join(plot_dir_name, prefix + "relEresol.png"))
    th1_relEresol.Write()
    dict_energy_relEresol_error[energy].append(relEresolFit.Get().Parameter(2))
    dict_energy_relEresol_error[energy].append(relEresolFit.Get().ParError(2))

    if args.cells:
        canvas_relEresol_cells = ROOT.TCanvas(prefix + "relE_resolution_cells", prefix + "relE_resolution_cells")
        fit_range_min = th1_relEresol_cells.GetMean() - 1 * th1_relEresol_cells.GetRMS()
        fit_range_max = th1_relEresol_cells.GetMean() + 1 * th1_relEresol_cells.GetRMS()
        #th1_relEresol_cells.Fit("gaus", "Q", "", fit_range_min, fit_range_max)
        relEresolFit_cells = th1_relEresol_cells.Fit("gaus", "SQ", "")
        print relEresolFit_cells
        th1_relEresol_cells.GetXaxis().SetTitle("(E_{Reco} - E_{Gen})/E_{Reco}")
        th1_relEresol_cells.GetXaxis().SetTitleOffset(1.2)
        th1_relEresol_cells.Draw()
        canvas_relEresol_cells.Print(os.path.join(plot_dir_name, prefix + "relEresol_cells.png"))
        th1_relEresol_cells.Write()
        dict_energy_cells_relEresol_error[energy].append(relEresolFit_cells.Get().Parameter(2))
        dict_energy_cells_relEresol_error[energy].append(relEresolFit_cells.Get().ParError(2))

    output_rootfile.Close()

def create_relEresol_graph(postfix):
    relEresol_vs_energy_graph = ROOT.TGraphErrors("ECAL_Energy_Resolution" + postfix)
    relEresol_vs_energy_graph.SetMarkerSize(1.5)
    relEresol_vs_energy_graph.SetMarkerStyle(args.markerStyle)
    relEresol_vs_energy_graph.SetMarkerColor(args.color)
    relEresol_vs_energy_graph.SetTitle("ECAL Energy Resolution " + postfix)
    relEresol_vs_energy_graph.GetXaxis().SetTitle("E_{Gen} [GeV]")
    relEresol_vs_energy_graph.GetXaxis().SetTitleOffset(1.2)
    relEresol_vs_energy_graph.GetYaxis().SetTitle("#scale[1.9]{#sigma}#left(#frac{E_{Reco} - E_{Gen}}{E_{Reco}}#right)")
    relEresol_vs_energy_graph.GetYaxis().SetTitleOffset(1.4)
    #relEresol_vs_energy_graph = ROOT.TGraphErrors(len(energies_gev_float), energies_gev_float, energy_resolutions, 0, energy_resolution_errors)
    #relEresol_vs_energy_graph.GetYaxis().SetLabelSize(0.06)
    #relEresol_vs_energy_graph.GetYaxis().SetTitleSize(0.06)
    #relEresol_vs_energy_graph.GetYaxis().SetTitleSize(0.05)
    return relEresol_vs_energy_graph


energies_gev_float.sort()
idx = 0
relEresol_vs_energy_graph = create_relEresol_graph("")
relEresol_vs_energy_cells_graph = create_relEresol_graph("cells")

for energy_float in energies_gev_float:
    energy_str = str(energy_float).replace(".", "dot")
    relEresol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_relEresol_error[energy_str][0])
    relEresol_vs_energy_graph.SetPointError(idx, 0, dict_energy_relEresol_error[energy_str][1])
    if args.cells:
        relEresol_vs_energy_cells_graph.SetPoint(idx, energy_float, dict_energy_cells_relEresol_error[energy_str][0])
        relEresol_vs_energy_cells_graph.SetPointError(idx, 0, dict_energy_cells_relEresol_error[energy_str][1])
    idx +=1

relEresol_vs_energy_graph.GetXaxis().SetLimits(0.2, 300)
relEresol_vs_energy_cells_graph.GetXaxis().SetLimits(0.2, 300)

#Fit energy resolution vs energy
#tf1_eResol_vs_e = TF1("tf1_eResol_vs_e", "sqrt([0] * [0] + pow([1] / sqrt(x), 2))", 5, 600)
output_rootfile_graph = ROOT.TFile(os.path.join(plot_dir_name, "relResol_vs_energy.root"), "recreate")
ROOT.gStyle.SetOptFit(000)

tf1_eResol_vs_e = ROOT.TF1("tf1_eResol_vs_e", "sqrt(pow([0]/x, 2) + pow([1]/sqrt(x), 2) + pow([2], 2))", energies_gev_float[0], energies_gev_float[-1])
tf1_eResol_vs_e.SetLineColor(args.color)
eResol_vs_e_fit = relEresol_vs_energy_graph.Fit(tf1_eResol_vs_e, "SQ")
a = str(round(eResol_vs_e_fit.Get().Parameter(0), 2))
b = str(round(eResol_vs_e_fit.Get().Parameter(1), 2))
c = str(round(eResol_vs_e_fit.Get().Parameter(2), 2))


canvas_relEresol_vs_energy = ROOT.TCanvas("relEresolution_vs_energy", "relEresolution_vs_energy")
canvas_relEresol_vs_energy.SetLogx(1)
relEresol_vs_energy_graph.Draw("ap")
relEresol_vs_energy_graph.Write()
eResol_vs_e_formula = "#color[%d]{#frac{%s}{E} #oplus #frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, a, b, c) 
eResol_vs_e_legend = gStyle.topRight_legend
eResol_vs_e_legend.AddEntry(ROOT.nullptr, eResol_vs_e_formula, "")
eResol_vs_e_legend.Draw()
canvas_relEresol_vs_energy.Print(os.path.join(plot_dir_name, "relEresolution_vs_energy.png"))
canvas_relEresol_vs_energy.Write()

if args.cells:
    eResol_vs_e_fit = relEresol_vs_energy_cells_graph.Fit(tf1_eResol_vs_e, "SQ")
    a = str(round(eResol_vs_e_fit.Get().Parameter(0), 2))
    b = str(round(eResol_vs_e_fit.Get().Parameter(1), 2))
    c = str(round(eResol_vs_e_fit.Get().Parameter(2), 2))

    canvas_relEresol_vs_energy_cells = ROOT.TCanvas("relEresolution_vs_energy_cells", "relEresolution_vs_energy_cells")
    canvas_relEresol_vs_energy_cells.SetLogx(1)
    relEresol_vs_energy_cells_graph.Draw("ap")
    relEresol_vs_energy_cells_graph.Write()
    eResol_vs_e_formula = "#color[%d]{#frac{%s}{E} #oplus #frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, a, b, c) 
    gStyle.topRight_legend.Clear()
    eResol_vs_e_legend_cells = gStyle.topRight_legend
    eResol_vs_e_legend_cells.AddEntry(ROOT.nullptr, eResol_vs_e_formula, "")
    eResol_vs_e_legend_cells.Draw()
    canvas_relEresol_vs_energy_cells.Print(os.path.join(plot_dir_name, "relEresolution_vs_energy_cells.png"))
    canvas_relEresol_vs_energy_cells.Write()

output_rootfile_graph.Close()
