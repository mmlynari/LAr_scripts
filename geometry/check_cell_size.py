import ROOT
from math import degrees, sqrt

inputFileName = "/eos/user/b/brfranco/rootfile_storage/fcc_analysis_ouput/220618_gamma_flat_1_100_noNoise_caloReco/fccsw_output_pdgID_22_pMin_1000_pMax_100000_thetaMin_50_thetaMax_130.root"

input_rootfile = ROOT.TFile(inputFileName)

tree = input_rootfile.Get('events')
tree.GetEntry(4)
for reference_cell_index in range(0, 100): 
    etaBin = tree.PositionedCaloClusterCells_etaBin[reference_cell_index]
    layer = tree.PositionedCaloClusterCells_layer[reference_cell_index]
    phiBin = tree.PositionedCaloClusterCells_phiBin[reference_cell_index]
    #if layer != 11:
    #    continue
    for cell_idx in range(len(tree.PositionedCaloClusterCells_phiBin)):
        #if abs(tree.PositionedCaloClusterCells_phiBin[cell_idx] - phiBin) == 1 and tree.PositionedCaloClusterCells_etaBin[cell_idx] == etaBin and tree.PositionedCaloClusterCells_layer[cell_idx] == layer:
        #    print("Theta: ", degrees(tree.PositionedCaloClusterCells_theta[cell_idx]))
        #    print("Layer: ", tree.PositionedCaloClusterCells_layer[cell_idx])
        #    print("Radius: ", sqrt(tree.PositionedCaloClusterCells_x[cell_idx]*tree.PositionedCaloClusterCells_x[cell_idx] + tree.PositionedCaloClusterCells_y[cell_idx]*tree.PositionedCaloClusterCells_y[cell_idx]))
        #    print("Phi difference: ", tree.PositionedCaloClusterCells_phi[cell_idx] - tree.PositionedCaloClusterCells_phi[reference_cell_index])
        #    print("Cell size in phi (mm): ", (tree.PositionedCaloClusterCells_phi[cell_idx] - tree.PositionedCaloClusterCells_phi[reference_cell_index]) * sqrt(tree.PositionedCaloClusterCells_x[cell_idx]*tree.PositionedCaloClusterCells_x[cell_idx] + tree.PositionedCaloClusterCells_y[cell_idx]*tree.PositionedCaloClusterCells_y[cell_idx]))
        #    print("-------------------------")
        if abs(tree.PositionedCaloClusterCells_etaBin[cell_idx] - etaBin) == 1 and tree.PositionedCaloClusterCells_phiBin[cell_idx] == phiBin and tree.PositionedCaloClusterCells_layer[cell_idx] == layer:
            print("Theta: ", degrees(tree.PositionedCaloClusterCells_theta[cell_idx]))
            print("Layer: ", tree.PositionedCaloClusterCells_layer[cell_idx])
            print("Radius: ", sqrt(tree.PositionedCaloClusterCells_x[cell_idx]*tree.PositionedCaloClusterCells_x[cell_idx] + tree.PositionedCaloClusterCells_y[cell_idx]*tree.PositionedCaloClusterCells_y[cell_idx]))
            print("Theta difference: ", tree.PositionedCaloClusterCells_theta[cell_idx] - tree.PositionedCaloClusterCells_theta[reference_cell_index])
            print("Cell size in Theta (mm): ", sqrt(tree.PositionedCaloClusterCells_x[cell_idx]*tree.PositionedCaloClusterCells_x[cell_idx] + tree.PositionedCaloClusterCells_y[cell_idx]*tree.PositionedCaloClusterCells_y[cell_idx]) * (tree.PositionedCaloClusterCells_theta[cell_idx] - tree.PositionedCaloClusterCells_theta[reference_cell_index]))
            print("-------------------------")
