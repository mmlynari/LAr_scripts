import ROOT
import os
import numpy as np
import argparse
from math import sqrt

parser = argparse.ArgumentParser()
parser.add_argument("-input", default = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/eGun1GeV_10kevt_originalGeometry_allCells/output_fullCalo_SimAndDigi_withCluster_noMagneticField_1GeV_pythiaFalse.root", help = "Name of the input file.", type = str)
parser.add_argument("-outputPostfix", default = "", help = "Postfix to append to the output folder.", type = str)
args = parser.parse_args()

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")

ROOT.gROOT.SetBatch(ROOT.kTRUE)

max_evt = -1
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/201209/fccsw_output_pythia_ee_Z_ee.root'
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/201216/fccsw_output_pythia_ee_Z_ee_all131all.root'
rootfile_path = args.input

if args.outputPostfix != "":
    plot_dir_name = 'perfPlots_'+os.path.basename(rootfile_path).replace('.root','') + "_" + args.outputPostfix
else:
    plot_dir_name = 'perfPlots_'+os.path.basename(rootfile_path).replace('.root','')
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)
output_rootfile_path = os.path.join(plot_dir_name, "perfHistograms.root")

#cutoff_dR = 0.015 
#cutoff_relE = 0.2
cutoff_dR = 1000
cutoff_relE = 1000

f = ROOT.TFile(rootfile_path)
events = f.Get("events")

th1_phiresol = ROOT.TH1F("phi_resolution", "phi_resolution", 100, -0.03, 0.03)
th1_thetaresol = ROOT.TH1F("theta_resolution", "theta_resolution", 100, -0.03, 0.03)
th1_angularresol = ROOT.TH1F("angular_resolution", "angular_resolution", 100, 0, 0.03)
th1_Eresol = ROOT.TH1F("energy_resolution", "energy_resolution", 100, -2, 2)
th1_relEresol = ROOT.TH1F("relative_energy_resolution", "relative_energy_resolution", 100, -1, 1)

evt = 0
for event in events:
    if evt >= max_evt and not max_evt == -1:
        break
    for caloCluster_idx in range(len(event.CaloClusters_energy)):
        best_genParticle_match_idx = -1
        best_d_R = 100000
        best_d_theta = 100000
        best_d_phi = 100000
        best_d_E = 100000
        best_d_relE = 10000
        for genParticle_idx in range(len(event.genParticle_status)):
            if event.genParticle_status[genParticle_idx] != 1:
                continue
            d_theta = event.CaloClusters_theta[caloCluster_idx] - event.genParticle_theta[genParticle_idx]
            d_phi = event.CaloClusters_phi[caloCluster_idx] - event.genParticle_phi[genParticle_idx]
            d_R = sqrt(d_theta * d_theta + d_phi * d_phi)
            d_E = event.CaloClusters_energy[caloCluster_idx] - event.genParticle_energy[genParticle_idx]
            d_relE = (event.CaloClusters_energy[caloCluster_idx] - event.genParticle_energy[genParticle_idx])/event.genParticle_energy[genParticle_idx]

            if d_R < best_d_R:
                best_genParticle_match_idx = caloCluster_idx
                best_d_R = d_R
                best_d_theta = d_theta
                best_d_phi = d_phi
                best_d_E = d_E
                best_d_relE = d_relE

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

    evt += 1
    if evt % 1000 == 0:
        print "Event processed: %d"%evt

output_rootfile = ROOT.TFile(output_rootfile_path, "recreate")
output_rootfile.cd()

ROOT.gStyle.SetOptFit(111)

canvas_phiresol = ROOT.TCanvas("phi_resolution", "phi_resolution")
th1_phiresol.Fit("gaus", "Q")
th1_phiresol.Draw()
th1_phiresol.GetXaxis().SetTitle("#Delta #Phi")
canvas_phiresol.Print(os.path.join(plot_dir_name, "phiresol.png"))
th1_phiresol.Write()

canvas_thetaresol = ROOT.TCanvas("theta_resolution", "theta_resolution")
th1_thetaresol.Fit("gaus", "Q")
th1_thetaresol.Draw()
th1_thetaresol.GetXaxis().SetTitle("#Delta #Theta")
canvas_thetaresol.Print(os.path.join(plot_dir_name, "thetaresol.png"))
th1_thetaresol.Write()

canvas_angularresol = ROOT.TCanvas("angular_resolution", "angular_resolution")
th1_angularresol.Draw()
th1_angularresol.GetXaxis().SetTitle("#Delta R = #sqrt{#Delta #Phi^{2} + #Delta #Theta^{2}}")
th1_angularresol.GetXaxis().SetTitleOffset(1.2)
canvas_angularresol.Print(os.path.join(plot_dir_name, "angularresol.png"))
th1_angularresol.Write()

canvas_Eresol = ROOT.TCanvas("E_resolution", "E_resolution")
th1_Eresol.Fit("gaus", "Q")
th1_Eresol.Draw()
th1_Eresol.GetXaxis().SetTitle("E_{Reco} - E_{Gen}")
canvas_Eresol.Print(os.path.join(plot_dir_name, "Eresol.png"))
th1_Eresol.Write()

canvas_relEresol = ROOT.TCanvas("relE_resolution", "relE_resolution")
fit_range_min = th1_relEresol.GetMean() - 1 * th1_relEresol.GetRMS()
fit_range_max = th1_relEresol.GetMean() + 1 * th1_relEresol.GetRMS()
#th1_relEresol.Fit("gaus", "Q", "", fit_range_min, fit_range_max)
th1_relEresol.Fit("gaus", "Q", "")
th1_relEresol.GetXaxis().SetTitle("(E_{Reco} - E_{Gen})/E_{Gen}")
th1_relEresol.GetXaxis().SetTitleOffset(1.2)
th1_relEresol.Draw()
canvas_relEresol.Print(os.path.join(plot_dir_name, "relEresol.png"))
th1_relEresol.Write()

output_rootfile.Close()
