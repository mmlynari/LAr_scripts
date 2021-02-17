import ROOT
import os, sys
import numpy as np

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")

#import gStyle
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetPadRightMargin(0.2)

#ROOT.gROOT.SetBatch(ROOT.kTRUE)

max_evt = 1
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/201209/fccsw_output_pythia_ee_Z_ee.root'
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/210217_caloReco/fccsw_output_pdgID_22_pMin_10000_pMax_10000_thetaMin_90_thetaMax_90.root'
rootfile_path = sys.argv[1]
plot_dir_name = 'plots_showerShape_' + os.path.basename(rootfile_path).replace('.root','')
energy = int(rootfile_path.split('_pMin_')[1].split("_")[0])
#energy = int(energy)/1000.0
energy_str_gev = str(energy/1000.0).replace(".","dot")
print energy_str_gev

if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

f = ROOT.TFile(rootfile_path)
events = f.Get("events")

th2_clusterCells_xy = ROOT.TH2F("th2_clusterCells_xy", "th2_clusterCells_xy", 200, -3000, +3000, 200, -3000, +3000)
th2_cells_xy = ROOT.TH2F("th2_cells_xy", "th2_cells_xy", 200, -3000, +3000, 200, -3000, +3000)

evt = 0
for event in events:
    if evt >= max_evt:
        break
    for cell_idx in xrange(len(event.ECalBarrelPositionedCells_energy)):
        th2_cells_xy.Fill(event.ECalBarrelPositionedCells_x[cell_idx], event.ECalBarrelPositionedCells_y[cell_idx], 100*event.ECalBarrelPositionedCells_energy[cell_idx]/float(energy))

    for caloCluster_idx in range(len(event.CaloClusters_energy)):
        firstCell = event.CaloClusters_firstCell[caloCluster_idx]
        lastCell = event.CaloClusters_lastCell[caloCluster_idx]
        for PositionedCaloClusterCell_idx in range(len(event.PositionedCaloClusterCells_x[firstCell:lastCell])):
            th2_clusterCells_xy.Fill(event.PositionedCaloClusterCells_x[firstCell+PositionedCaloClusterCell_idx], event.PositionedCaloClusterCells_y[firstCell+PositionedCaloClusterCell_idx], 100*event.PositionedCaloClusterCells_energy[firstCell+PositionedCaloClusterCell_idx]/float(energy))
        #print event.CaloClusters_x[caloCluster_idx]
        #print event.CaloClusters_y[caloCluster_idx]
        #print "------------------"
        break

    evt += 1
    if evt % 100 == 0:
        print(evt)

canvas_cells_xy = ROOT.TCanvas("th2_cells_xy", "th2_cells_xy")
th2_cells_xy.GetZaxis().SetTitle("Energy (%)")
th2_cells_xy.Draw('colz')
canvas_cells_xy.Print(os.path.join(plot_dir_name, 'th2_cells_xy.png'))
canvas_cells_xy.Print(os.path.join(plot_dir_name, 'th2_cells_xy.root'))

canvas_clusterCells_xy = ROOT.TCanvas("th2_clusterCells_xy", "th2_clusterCells_xy")
th2_clusterCells_xy.Draw('colz')
canvas_clusterCells_xy.Print(os.path.join(plot_dir_name, 'th2_clusterCells_xy.png'))
canvas_clusterCells_xy.Print(os.path.join(plot_dir_name, 'th2_clusterCells_xy.root'))
