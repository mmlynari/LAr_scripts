import ROOT
import os
import numpy as np

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")

ROOT.gROOT.SetBatch(ROOT.kTRUE)

#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/201209/fccsw_output_pythia_ee_Z_ee.root'
#rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/outputs/201216/fccsw_output_pythia_ee_Z_ee_all131all.root'
rootfile_path = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/210607/FCCAnalyses/outputs/211020/output_fullCalo_SimAndDigi_withCluster_MagneticField_False_pMin_3000_MeV_ThetaMinMax_90.25_90.25_pdgId_111_pythiaFalse.root'
plot_dir_name = 'plots_clusterShape_' + os.path.basename(rootfile_path).replace('.root','')
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

f = ROOT.TFile(rootfile_path)
events = f.Get("events")
#print "Chaning tree to matrix..."
#events_asmatrix = events.AsMatrix()
#print("Tree converted to a numpy array:\n{}\n".format(array))

th2_clusterCell1_xy = ROOT.TH2F("th2_clusterCell1_xy", "th2_clusterCell1_xy", 3000, -3000, +3000, 3000, -3000, +3000)
th2_clusterCell2_xy = ROOT.TH2F("th2_clusterCell2_xy", "th2_clusterCell2_xy", 3000, -3000, +3000, 3000, -3000, +3000)
evt = 0
n_skip = 1
max_evt = 1
for event in events:
    print(evt)
    if evt < n_skip:
        evt += 1
        max_evt += 1
        continue
    if evt >= max_evt:
        break
    if len(event.CaloClusters_energy) != 2:
        print("we dont have two clusters")
        continue
    for caloCluster_idx in range(len(event.CaloClusters_energy)):
        firstCell = event.CaloClusters_firstCell[caloCluster_idx]
        lastCell = event.CaloClusters_lastCell[caloCluster_idx]
        print(firstCell, lastCell)
        for PositionedCaloClusterCell_idx in range(lastCell - firstCell):
            if caloCluster_idx == 0:
                th2_clusterCell1_xy.Fill(event.PositionedCaloClusterCells_x[firstCell+PositionedCaloClusterCell_idx], event.PositionedCaloClusterCells_y[firstCell+PositionedCaloClusterCell_idx])
            if caloCluster_idx == 1:
                th2_clusterCell2_xy.Fill(event.PositionedCaloClusterCells_x[firstCell+PositionedCaloClusterCell_idx], event.PositionedCaloClusterCells_y[firstCell+PositionedCaloClusterCell_idx])
        print(event.CaloClusters_x[caloCluster_idx])
        print(event.CaloClusters_y[caloCluster_idx])
        print("------------------")

    evt += 1
    if evt % 100 == 0:
        print(evt)

canvas_clusterCell1_xy = ROOT.TCanvas("th2_clusterCell1_xy", "th2_clusterCell1_xy")
th2_clusterCell1_xy.Draw()
canvas_clusterCell1_xy.Print('th2_clusterCell1_xy.png')

canvas_clusterCell2_xy = ROOT.TCanvas("th2_clusterCell2_xy", "th2_clusterCell2_xy")
th2_clusterCell2_xy.Draw()
canvas_clusterCell2_xy.Print('th2_clusterCell2_xy.png')
