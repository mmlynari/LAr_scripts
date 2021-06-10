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
rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/210401_caloReco_mip/output_fullCalo_SimAndDigi_withCluster_MagneticField_False_pMin_10GeV_ThetaMinMax_90.0_90.0_pdgId_13_pythiaFalse.root'
postfix = ''
#if sys.argv[2]:
#    postfix = sys.argv[2]
#else:
#    postfix = ''

noise_root_file_path = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCSW_201207_geometry/LAr_scripts/geometry/noise_capa/elecNoise_ecalBarrelFCCee.root"
noise_root_file = ROOT.TFile(noise_root_file_path, 'r')
th1_name_template = 'h_elecNoise_fcc_LAYERID'
layer_radius_edges = [2160, 2173.6787835669493, 2206.0574015524803, 2239.0593055298277, 2272.657342921452, 2306.825466503636, 2341.5387200310156, 2376.7732184243205, 2412.5061235473704, 2448.71561648303, 2485.3808671070947, 2522.482001655514, 2560.000068884735] # needed to know from which layer to take the noise histogram, obtained from create_capacitance.py
noise_histos = []
for layer_idx in range(len(layer_radius_edges) - 1) :
    noise_histo = noise_root_file.Get(th1_name_template.replace('LAYERID', str(layer_idx + 1)))
    noise_histo.SetDirectory(0)
    noise_histos.append(noise_histo)

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
tprof_signalOverNoise_vs_depth = ROOT.TProfile("tprof_signalOverNoise_vs_depth_" + energy_str_gev, "tprof_signalOverNoise_vs_depth_" + energy_str_gev, len(layer_radius_edges) - 1, np.array(layer_radius_edges, 'd'), "s")
tprof_noise_vs_depth = ROOT.TProfile("tprof_noise_vs_depth_" + energy_str_gev, "tprof_noise_vs_depth_" + energy_str_gev, len(layer_radius_edges) - 1, np.array(layer_radius_edges, 'd'), "s")

min_evt = 0
max_evt = n_entries * 10
evt = 0
for event in events:
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
        noise_value = noise_histos[final_layer_idx].GetBinContent(1)
        th2_cells_xy.Fill(event.ECalBarrelPositionedCells_x[cell_idx], event.ECalBarrelPositionedCells_y[cell_idx], 1000*100*event.ECalBarrelPositionedCells_energy[cell_idx]/(float(energy)*n_entries))
        tprof_energy_vs_depth.Fill(sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx]), 1000*event.ECalBarrelPositionedCells_energy[cell_idx])
        tprof_signalOverNoise_vs_depth.Fill(sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx]), event.ECalBarrelPositionedCells_energy[cell_idx] / noise_value)
        tprof_noise_vs_depth.Fill(sqrt(event.ECalBarrelPositionedCells_x[cell_idx] * event.ECalBarrelPositionedCells_x[cell_idx] + event.ECalBarrelPositionedCells_y[cell_idx] * event.ECalBarrelPositionedCells_y[cell_idx]), 1000*noise_value)

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

canvas_energy_vs_depth = ROOT.TCanvas("tprof_energy_vs_depth_canvas_" + energy_str_gev, "tprof_energy_vs_depth_canvas_" + energy_str_gev)
print("Integral of TProfile (average total energy deposit): ", tprof_energy_vs_depth.Integral())
#tprof_energy_vs_depth.Scale(1/tprof_energy_vs_depth.Integral())
#tprof_energy_vs_depth.SetMarkerSize(4)
tprof_energy_vs_depth.SetLineWidth(3)
tprof_energy_vs_depth.GetXaxis().SetTitle("Radial depth (mm)")
tprof_energy_vs_depth.GetYaxis().SetTitle("Average energy deposit [MeV]")
tprof_energy_vs_depth.Draw()
tprof_energy_vs_depth.Write()
canvas_energy_vs_depth.Write()
canvas_energy_vs_depth.Print(os.path.join(plot_dir_name, "tprof_energy_vs_depth_" + energy_str_gev + ".png"))

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

#canvas_clusterCells_xy = ROOT.TCanvas("th2_clusterCells_xy", "th2_clusterCells_xy")
##th2_clusterCells_xy.Smooth()
#th2_clusterCells_xy.GetZaxis().SetTitle("Energy (a.u.)")
#th2_clusterCells_xy.GetZaxis().SetTitleOffset(1.35)
#th2_clusterCells_xy.Draw('colz')
#canvas_cells_xy.Write()
#canvas_clusterCells_xy.Print(os.path.join(plot_dir_name, energy_str_gev + '_th2_clusterCells_xy.png'))
##canvas_clusterCells_xy.Print(os.path.join(plot_dir_name, 'th2_clusterCells_xy.root'))

output_rootfile.Close()
