import ROOT
import os, sys
import numpy as np
from datetime import date

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")

#import gStyle
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetPadRightMargin(0.2)
ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)

ROOT.gROOT.SetBatch(ROOT.kTRUE)

#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/201209/fccsw_output_pythia_ee_Z_ee.root'
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/210217_caloReco/fccsw_output_pdgID_22_pMin_10000_pMax_10000_thetaMin_90_thetaMax_90.root'
rootfile_path = sys.argv[1]
if not os.path.isfile(rootfile_path):
    print("Provided roootfile does not exists")
    sys.exit(1)
plot_dir_name = 'plots_showerShape_' + date.today().strftime("%y%m%d")  + "_" + os.path.basename(rootfile_path).replace('.root','')
energy = int(rootfile_path.split('_pMin_')[1].split("_")[0])
#energy = int(energy)/1000.0
energy_str_gev = str(energy/1000.0).replace(".","dot")
print(energy_str_gev)

if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)


f = ROOT.TFile(rootfile_path)
events = f.Get("events")

output_rootfile_path = os.path.join(plot_dir_name, "shower_shape_%s.root"%energy_str_gev)
output_rootfile = ROOT.TFile(output_rootfile_path, "recreate")

min_evt = 10
max_evt = 20
evt = 0
for event in events:
    evt += 1
    if evt % 100 == 0:
        print(evt)
    if evt < min_evt:
        continue
    if evt >= max_evt:
        break
    th2_clusterCells_xy = ROOT.TH2F(energy_str_gev + "th2_clusterCells_xy_energy_evt%d"%evt, energy_str_gev + "th2_clusterCells_xy_energy_evt%d"%evt, 300, -3000, +3000, 300, -3000, +3000)
    th2_cells_xy = ROOT.TH2F(energy_str_gev + "th2_cells_xy_energy%d"%evt, energy_str_gev + "th2_cells_xy_energy%d"%evt, 300, -3000, +3000, 300, -3000, +3000)
    for cell_idx in range(len(event.ECalBarrelPositionedCells_energy)):
        th2_cells_xy.Fill(event.ECalBarrelPositionedCells_x[cell_idx], event.ECalBarrelPositionedCells_y[cell_idx], 1000*100*event.ECalBarrelPositionedCells_energy[cell_idx]/float(energy))

    for caloCluster_idx in range(len(event.CaloClusters_energy)):
        firstCell = event.CaloClusters_firstCell[caloCluster_idx]
        lastCell = event.CaloClusters_lastCell[caloCluster_idx]
        for PositionedCaloClusterCell_idx in range(len(event.PositionedCaloClusterCells_x[firstCell:lastCell])):
            th2_clusterCells_xy.Fill(event.PositionedCaloClusterCells_x[firstCell+PositionedCaloClusterCell_idx], event.PositionedCaloClusterCells_y[firstCell+PositionedCaloClusterCell_idx], 100*1000*event.PositionedCaloClusterCells_energy[firstCell+PositionedCaloClusterCell_idx]/float(energy))
        #print event.CaloClusters_x[caloCluster_idx]
        #print event.CaloClusters_y[caloCluster_idx]
        #print "------------------"
        break


    canvas_cells_xy = ROOT.TCanvas("th2_cells_xy", "th2_cells_xy")
    th2_cells_xy.Smooth()
    th2_cells_xy.GetZaxis().SetTitle("Energy (%)")
    th2_cells_xy.GetZaxis().SetTitleOffset(1.35)
    th2_cells_xy.Draw('colz')
    canvas_cells_xy.Write()
    canvas_cells_xy.Print(os.path.join(plot_dir_name, energy_str_gev + 'th2_cells_xy_evt%d.png'%evt))
    #canvas_cells_xy.Print(os.path.join(plot_dir_name, 'th2_cells_xy.root'))

    canvas_clusterCells_xy = ROOT.TCanvas("th2_clusterCells_xy", "th2_clusterCells_xy")
    th2_clusterCells_xy.Smooth()
    th2_clusterCells_xy.GetZaxis().SetTitle("Energy (%)")
    th2_clusterCells_xy.GetZaxis().SetTitleOffset(1.35)
    th2_clusterCells_xy.Draw('colz')
    canvas_cells_xy.Write()
    canvas_clusterCells_xy.Print(os.path.join(plot_dir_name, energy_str_gev + 'th2_clusterCells_xy_evt%d.png'%evt))
    #canvas_clusterCells_xy.Print(os.path.join(plot_dir_name, 'th2_clusterCells_xy.root'))

output_rootfile.Close()
