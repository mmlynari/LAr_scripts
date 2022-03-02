import ROOT
import os, sys
import numpy as np
from datetime import date
from math import sqrt

print("Careful: only works for eta == 0 at the moment")
#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")

#import gStyle
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetPadRightMargin(0.2)
ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
ROOT.gStyle.SetPadTickY(1)

ROOT.gROOT.SetBatch(ROOT.kTRUE)

#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/201209/fccsw_output_pythia_ee_Z_ee.root'
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/210217_caloReco/fccsw_output_pdgID_22_pMin_10000_pMax_10000_thetaMin_90_thetaMax_90.root'
#rootfile_path = sys.argv[1]
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/210401_caloReco_mip/output_fullCalo_SimAndDigi_withCluster_MagneticField_False_pMin_10GeV_ThetaMinMax_90.0_90.0_pdgId_13_pythiaFalse.root'
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/210607/output_fullCalo_SimAndDigi_withCluster_MagneticField_False_pMin_10000_MeV_ThetaMinMax_90.25_90.25_pdgId_13_pythiaFalse_phiFixed.root'
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/210607/FCCAnalyses/211001_energies_10kevt_caloReco/fccsw_output_pdgID_22_pMin_1000_pMax_1000_thetaMin_90_thetaMax_90.root'
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/211210/FCCAnalyses/220126_mip_10kevt_slidingWindow_noNoise_caloReco/fccsw_output_pdgID_13_pMin_10000_pMax_10000_thetaMin_90_thetaMax_90.root'
rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/211210/FCCAnalyses/mip_updatedSF/output_fullCalo_SimAndDigi_withCluster_MagneticField_False_pMin_20000_MeV_ThetaMinMax_90.25_90.25_pdgId_13_pythiaFalse_NoiseFalse.root'
postfix = 'mip'
#if sys.argv[2]:
#    postfix = sys.argv[2]
#else:
#    postfix = ''

#noise_root_file_path = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCSW_201207_geometry/LAr_scripts/geometry/noise_capa/elecNoise_ecalBarrelFCCee.root"
#noise_root_file_path = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/210927/LAr_scripts/geometry/noise_constant_vs_capa/elecNoise_ecalBarrelFCCee.root"
noise_root_file_path = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/210927/LAr_scripts/geometry/noise_capa_220216/elecNoise_ecalBarrelFCCee.root"
noise_root_file = ROOT.TFile(noise_root_file_path, 'r')
layer_radius_edges = [2160, 2173.6787835669493, 2206.0574015524803, 2239.0593055298277, 2272.657342921452, 2306.825466503636, 2341.5387200310156, 2376.7732184243205, 2412.5061235473704, 2448.71561648303, 2485.3808671070947, 2522.482001655514, 2560.000068884735] # needed to know from which layer to take the noise histogram, obtained from create_capacitance.py
SFfcc = [0.36504678560781667] * 1 + [0.09974087165838573] * 1 + [0.12392336840429007] * 1 + [0.1413266332223572] * 1 + [0.15415123193238958] * 1 + [0.1639900875460671] * 1 + [0.17156597031962592] * 1 + [0.17810674932424356] * 1 + [0.18340048249397345] * 1 + [0.18855877603870688] * 1 + [0.19307873042890955] * 1 + [0.21746137329706489] * 1
th1_noise_name_template = 'h_elecNoise_fcc_LAYERID'
th1_energy_name_template = 'h_energy_layer_LAYERID'
noise_histos = []
energy_histos = []
for layer_idx in range(len(layer_radius_edges) - 1) :
    noise_histo = noise_root_file.Get(th1_noise_name_template.replace('LAYERID', str(layer_idx + 1)))
    noise_histo.SetDirectory(0)
    noise_histos.append(noise_histo)
    energy_histo = ROOT.TH1F("th1_energy_layer_" + str(layer_idx + 1), "th1_energy_layer_" + str(layer_idx + 1), 50, 0, 30)
    energy_histo.SetDirectory(0)
    energy_histos.append(energy_histo)

if not os.path.isfile(rootfile_path):
    print("Provided roootfile does not exists")
    sys.exit(1)
plot_dir_name = 'plots_signalOverNoisePerLayer_' + postfix + date.today().strftime("%y%m%d")  #+ "_" + os.path.basename(rootfile_path).replace('.root','')
energy = int(rootfile_path.split('_pMin_')[1].split("_")[0].replace('GeV', ''))
#energy = int(energy)/1000.0
energy_str_gev = str(energy/1000.0).replace(".","dot")
print(energy_str_gev)

if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

f = ROOT.TFile(rootfile_path)
events = f.Get("events")
n_entries = events.GetEntries()

output_rootfile_path = os.path.join(plot_dir_name, "averaged_shower_shape_%s.root"%energy_str_gev)
output_rootfile = ROOT.TFile(output_rootfile_path, "recreate")

#th2_clusterCells_xy = ROOT.TH2F(energy_str_gev + "_th2_clusterCells_xy_energy_evt", energy_str_gev + "_th2_clusterCells_xy_energy_evt", 300, -3000, +3000, 300, -3000, +3000)
th2_cells_xy = ROOT.TH2F("th2_cells_xy_energy_" + energy_str_gev, "th2_cells_xy_energy_" + energy_str_gev, 150, -3000, +3000, 150, -3000, +3000)
#tprof_energy_vs_depth = ROOT.TProfile("tprof_energy_vs_depth_" + energy_str_gev, "tprof_energy_vs_depth_" + energy_str_gev, 20, 2150, 2550, "s")
tprof_energy_vs_depth = ROOT.TProfile("tprof_energy_vs_depth_" + energy_str_gev, "tprof_energy_vs_depth_" + energy_str_gev, len(layer_radius_edges) - 1, np.array(layer_radius_edges, 'd'), "s")
tprof_uncalibratedEnergy_vs_depth = ROOT.TProfile("tprof_uncalibratedEnergy_vs_depth_" + energy_str_gev, "tprof_uncalibratedEnergy_vs_depth_" + energy_str_gev, len(layer_radius_edges) - 1, np.array(layer_radius_edges, 'd'), "s")
tprof_signalOverNoise_vs_depth = ROOT.TProfile("tprof_signalOverNoise_vs_depth_" + energy_str_gev, "tprof_signalOverNoise_vs_depth_" + energy_str_gev, len(layer_radius_edges) - 1, np.array(layer_radius_edges, 'd'), "s")
th1_signalOverNoise_layer3 = ROOT.TH1F("th1_signalOverNoise_layer3_" + energy_str_gev, "th1_signalOverNoise_layer3_" + energy_str_gev, 50, 0, 30)
tprof_noise_vs_depth = ROOT.TProfile("tprof_noise_vs_depth_" + energy_str_gev, "tprof_noise_vs_depth_" + energy_str_gev, len(layer_radius_edges) - 1, np.array(layer_radius_edges, 'd'), "s")
tprof_nFiredCell_vs_depth = ROOT.TProfile("tprof_nFiredCell_vs_depth_" + energy_str_gev, "tprof_nFiredCell_vs_depth_" + energy_str_gev, len(layer_radius_edges) - 1, np.array(layer_radius_edges, 'd'), "s")

min_evt = 0
max_evt = n_entries * 10
evt = 0
for event in events:
    n_fired_cells = []
    for layer_idx in range(len(layer_radius_edges) - 1):
        n_fired_cells.append(0)
    evt += 1
    if evt % 100 == 0:
        print(evt)
    if evt < min_evt:
        continue
    if evt >= max_evt:
        break
    for cell_idx in range(len(event.ECalBarrelPositionedCells_energy)):
        hit_radius = sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx])
        #hit_radius = sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx] + event.ECalBarrelPositionedCells_z[cell_idx] * event.ECalBarrelPositionedCells_z[cell_idx])
        final_layer_idx = -1
        for layer_idx in range(len(layer_radius_edges) - 1):
            layer_radius_lower_edge = layer_radius_edges[layer_idx]
            layer_radius_higher_edge = layer_radius_edges[layer_idx + 1]
            if hit_radius > layer_radius_lower_edge and hit_radius < layer_radius_higher_edge:
                final_layer_idx = layer_idx 
        noise_value = noise_histos[final_layer_idx].GetBinContent(1)/1000.0
        th2_cells_xy.Fill(event.ECalBarrelPositionedCells_x[cell_idx], event.ECalBarrelPositionedCells_y[cell_idx], 1000*100*event.ECalBarrelPositionedCells_energy[cell_idx]/(float(energy)*n_entries))
        tprof_energy_vs_depth.Fill(sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx]), 1000*event.ECalBarrelPositionedCells_energy[cell_idx])
        tprof_uncalibratedEnergy_vs_depth.Fill(sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx]), 1000*event.ECalBarrelPositionedCells_energy[cell_idx]*SFfcc[final_layer_idx])
        tprof_signalOverNoise_vs_depth.Fill(sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx]), event.ECalBarrelPositionedCells_energy[cell_idx] / noise_value)
        energy_histos[final_layer_idx].Fill(event.ECalBarrelPositionedCells_energy[cell_idx] * 1000)
        if final_layer_idx == 2:
            #print sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx])
            th1_signalOverNoise_layer3.Fill(event.ECalBarrelPositionedCells_energy[cell_idx] / noise_value)
        n_fired_cells[final_layer_idx] += 1
    for layer_idx in range(len(layer_radius_edges) - 1):
        tprof_nFiredCell_vs_depth.Fill((layer_radius_edges[layer_idx] + layer_radius_edges[layer_idx + 1])/2.0, n_fired_cells[layer_idx])

    #for caloCluster_idx in range(len(event.CaloClusters_energy)):
    #    firstCell = event.CaloClusters_firstCell[caloCluster_idx]
    #    lastCell = event.CaloClusters_lastCell[caloCluster_idx]
    #    for PositionedCaloClusterCell_idx in range(len(event.PositionedCaloClusterCells_x[firstCell:lastCell])):
    #        th2_clusterCells_xy.Fill(event.PositionedCaloClusterCells_x[firstCell+PositionedCaloClusterCell_idx], event.PositionedCaloClusterCells_y[firstCell+PositionedCaloClusterCell_idx], 100*1000*event.PositionedCaloClusterCells_energy[firstCell+PositionedCaloClusterCell_idx]/float(energy))
    #    break


canvas_cells_xy = ROOT.TCanvas("th2_cells_xy_" + energy_str_gev, "th2_cells_xy_" + energy_str_gev)
th2_cells_xy.Smooth()
th2_cells_xy.GetZaxis().SetTitle("Energy (a.u.)")
th2_cells_xy.GetZaxis().SetTitleOffset(1.35)
th2_cells_xy.Draw('colz')
canvas_cells_xy.Write()
canvas_cells_xy.Print(os.path.join(plot_dir_name, "th2_cells_xy_" + energy_str_gev + ".png"))
#canvas_cells_xy.Print(os.path.join(plot_dir_name, 'th2_cells_xy.root'))

canvas_signalOverNoise_vs_depth = ROOT.TCanvas("tprof_signalOverNoise_vs_depth_canvas_" + energy_str_gev, "tprof_signalOverNoise_vs_depth_canvas_" + energy_str_gev)
#tprof_signalOverNoise_vs_depth.Scale(1/tprof_signalOverNoise_vs_depth.Integral())
#tprof_signalOverNoise_vs_depth.SetMarkerSize(4)
tprof_signalOverNoise_vs_depth.SetLineWidth(3)
tprof_signalOverNoise_vs_depth.GetXaxis().SetTitle("Radial depth (mm)")
tprof_signalOverNoise_vs_depth.GetYaxis().SetTitle("Signal over noise")
tprof_signalOverNoise_vs_depth.Draw()
tprof_signalOverNoise_vs_depth.Write()
canvas_signalOverNoise_vs_depth.Write()
canvas_signalOverNoise_vs_depth.Print(os.path.join(plot_dir_name, "tprof_signalOverNoise_vs_depth_" + energy_str_gev + ".png"))

canvas_energy_vs_depth = ROOT.TCanvas("tprof_energy_vs_depth_canvas_" + energy_str_gev, "tprof_energy_vs_depth_canvas_" + energy_str_gev)
print(("Integral of TProfile (average total energy deposit): ", tprof_energy_vs_depth.Integral()))
#tprof_energy_vs_depth.Scale(1/tprof_energy_vs_depth.Integral())
#tprof_energy_vs_depth.SetMarkerSize(4)
tprof_energy_vs_depth.SetLineWidth(3)
tprof_energy_vs_depth.GetXaxis().SetTitle("Radial depth (mm)")
tprof_energy_vs_depth.GetYaxis().SetTitle("Average energy deposit [MeV]")
tprof_energy_vs_depth.Draw()
tprof_energy_vs_depth.Write()
canvas_energy_vs_depth.Write()
canvas_energy_vs_depth.Print(os.path.join(plot_dir_name, "tprof_energy_vs_depth_" + energy_str_gev + ".png"))

canvas_uncalibratedEnergy_vs_depth = ROOT.TCanvas("tprof_uncalibratedEnergy_vs_depth_canvas_" + energy_str_gev, "tprof_uncalibratedEnergy_vs_depth_canvas_" + energy_str_gev)
print("Integral of TProfile (average total uncalibratedEnergy deposit): ", tprof_uncalibratedEnergy_vs_depth.Integral())
#tprof_uncalibratedEnergy_vs_depth.Scale(1/tprof_uncalibratedEnergy_vs_depth.Integral())
#tprof_uncalibratedEnergy_vs_depth.SetMarkerSize(4)
tprof_uncalibratedEnergy_vs_depth.SetLineWidth(3)
tprof_uncalibratedEnergy_vs_depth.GetXaxis().SetTitle("Radial depth (mm)")
tprof_uncalibratedEnergy_vs_depth.GetYaxis().SetTitle("Average energy deposit [MeV]")
tprof_uncalibratedEnergy_vs_depth.Draw()
tprof_uncalibratedEnergy_vs_depth.Write()
canvas_uncalibratedEnergy_vs_depth.Write()
canvas_uncalibratedEnergy_vs_depth.Print(os.path.join(plot_dir_name, "tprof_uncalibratedEnergy_vs_depth_" + energy_str_gev + ".png"))

canvas_noise_vs_depth = ROOT.TCanvas("tprof_noise_vs_depth_canvas_" + energy_str_gev, "tprof_noise_vs_depth_canvas_" + energy_str_gev)
#tprof_noise_vs_depth.Scale(1/tprof_noise_vs_depth.Integral())
#tprof_noise_vs_depth.SetMarkerSize(4)
tprof_noise_vs_depth.SetLineWidth(3)
tprof_noise_vs_depth.GetXaxis().SetTitle("Radial depth (mm)")
tprof_noise_vs_depth.GetYaxis().SetTitle("Noise [MeV]")
tprof_noise_vs_depth.Draw()
tprof_noise_vs_depth.Write()
canvas_noise_vs_depth.Write()
canvas_noise_vs_depth.Print(os.path.join(plot_dir_name, "tprof_noise_vs_depth_" + energy_str_gev + ".png"))

canvas_nFiredCell_vs_depth = ROOT.TCanvas("tprof_nFiredCell_vs_depth_canvas_" + energy_str_gev, "tprof_nFiredCell_vs_depth_canvas_" + energy_str_gev)
#tprof_nFiredCell_vs_depth.Scale(1/tprof_nFiredCell_vs_depth.Integral())
#tprof_nFiredCell_vs_depth.SetMarkerSize(4)
tprof_nFiredCell_vs_depth.SetLineWidth(3)
tprof_nFiredCell_vs_depth.GetXaxis().SetTitle("Radial depth (mm)")
tprof_nFiredCell_vs_depth.GetYaxis().SetTitle("Number of fired cells")
tprof_nFiredCell_vs_depth.Draw()
tprof_nFiredCell_vs_depth.Write()
canvas_nFiredCell_vs_depth.Write()
canvas_nFiredCell_vs_depth.Print(os.path.join(plot_dir_name, "tprof_nFiredCell_vs_depth_" + energy_str_gev + ".png"))

#canvas_clusterCells_xy = ROOT.TCanvas("th2_clusterCells_xy", "th2_clusterCells_xy")
##th2_clusterCells_xy.Smooth()
#th2_clusterCells_xy.GetZaxis().SetTitle("Energy (a.u.)")
#th2_clusterCells_xy.GetZaxis().SetTitleOffset(1.35)
#th2_clusterCells_xy.Draw('colz')
#canvas_cells_xy.Write()
#canvas_clusterCells_xy.Print(os.path.join(plot_dir_name, energy_str_gev + '_th2_clusterCells_xy.png'))
##canvas_clusterCells_xy.Print(os.path.join(plot_dir_name, 'th2_clusterCells_xy.root'))

ROOT.gStyle.SetOptStat(1101)
canvas_signalOverNoise_layer3 = ROOT.TCanvas("th1_signalOverNoise_layer3_canvas_" + energy_str_gev, "th1_signalOverNoise_layer3_canvas_" + energy_str_gev)
#th1_signalOverNoise_layer3.Scale(1/th1_signalOverNoise_layer3.Integral())
#th1_signalOverNoise_layer3.SetMarkerSize(4)
#th1_signalOverNoise_layer3.SetLineWidth(3)
th1_signalOverNoise_layer3.GetXaxis().SetTitle("Signal over noise")
th1_signalOverNoise_layer3.GetYaxis().SetTitle("count")
print(th1_signalOverNoise_layer3.GetEntries())
th1_signalOverNoise_layer3.Draw()
th1_signalOverNoise_layer3.Write()
canvas_signalOverNoise_layer3.Write()
canvas_signalOverNoise_layer3.Print(os.path.join(plot_dir_name, "th1_signalOverNoise_layer3_" + energy_str_gev + ".png"))

for layer_idx in range(len(layer_radius_edges) - 1) :
    canvas_energy_layer = ROOT.TCanvas("th1_energy_layer_" + str(layer_idx + 1) + "_canvas_" + energy_str_gev, "th1_energy_layer_" + str(layer_idx + 1) + "_canvas_" + energy_str_gev)
    #th1_energy_layer_.Scale(1/th1_energy_layer_.Integral())
    #th1_energy_layer_.SetMarkerSize(4)
    #th1_energy_layer_.SetLineWidth(3)
    energy_histos[layer_idx].GetXaxis().SetRange(1, energy_histos[layer_idx].GetNbinsX() + 1)
    energy_histos[layer_idx].GetXaxis().SetTitle("Energy deposit [MeV]")
    energy_histos[layer_idx].GetYaxis().SetTitle("Count")
    energy_histos[layer_idx].Draw()
    energy_histos[layer_idx].Write()
    canvas_energy_layer.Write()
    canvas_energy_layer.Print(os.path.join(plot_dir_name, "th1_energy_layer_" + str(layer_idx + 1) + "_" + energy_str_gev + ".png"))

output_rootfile.Close()
