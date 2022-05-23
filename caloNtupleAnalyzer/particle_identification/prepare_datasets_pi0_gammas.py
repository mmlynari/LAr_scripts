import os
import ROOT

#input_dir = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/210607/FCCAnalyses/210927_energies_10kevt_PARTICLE_caloReco/'
#input_dir = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/210607/FCCAnalyses/210927_energies_10kevt_pi0_caloReco/'
input_dir = '/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/210607/FCCAnalyses/211020_energies_10kevt_both_caloReco/'
file_name_template = 'fccsw_output_pdgID_PDGID_pMin_ENERGY_pMax_ENERGY_thetaMin_90_thetaMax_90.root'
energy_list = [1000, 2000, 5000, 10000, 20000, 30000, 45000, 100000]
#energy_list = [1000]

pdgIDs = ['111', '22']

max_evt_spatial_extent = 1000

for pdgID in pdgIDs:
    for energy in energy_list:
        max_eta_extent = -1
        max_phi_extent = -1
        derive_spatial_extent = True

        input_file = os.path.join(input_dir, file_name_template.replace('PDGID', pdgID).replace('ENERGY', str(energy)))
        print("Treating ", energy, " GeV with PID ", pdgID)
        input_rootfile = ROOT.TFile(input_file)
        tree = input_rootfile.Get('events')
        #print(tree.GetEntries())

        #FIXME we need the cluster cells to derive the spatial extent ortherwise on single fired cell far away from energy deposit spoils it all... 
        if derive_spatial_extent:
            evt_count_spatial_extent = 0
            for entry in tree:
                min_eta_bin = 11111111111
                max_eta_bin = -111111111111
                max_layer = -11111111111
                min_phi_bin = 11111111111
                max_phi_bin = -11111111111

                for caloCluster_idx in range(len(entry.CaloClusters_energy)):
                    firstCell = entry.CaloClusters_firstCell[caloCluster_idx]
                    lastCell = entry.CaloClusters_lastCell[caloCluster_idx]
                    for cell_idx in range(len(entry.PositionedCaloClusterCells_x[firstCell:lastCell])):
                        if entry.ECalBarrelPositionedCells_phiBin[cell_idx] > max_phi_bin:
                            max_phi_bin = entry.ECalBarrelPositionedCells_phiBin[cell_idx]
                        elif entry.ECalBarrelPositionedCells_phiBin[cell_idx] < min_phi_bin:
                            min_phi_bin = entry.ECalBarrelPositionedCells_phiBin[cell_idx]
                        if entry.ECalBarrelPositionedCells_etaBin[cell_idx] < min_eta_bin:
                            min_eta_bin = entry.ECalBarrelPositionedCells_etaBin[cell_idx]
                        elif entry.ECalBarrelPositionedCells_etaBin[cell_idx] > max_eta_bin:
                            max_eta_bin = entry.ECalBarrelPositionedCells_etaBin[cell_idx]
                        if entry.ECalBarrelPositionedCells_layer[cell_idx] > max_layer:
                            max_layer = entry.ECalBarrelPositionedCells_layer[cell_idx]
                if max_eta_bin - min_eta_bin + 1 > max_eta_extent:
                    max_eta_extent = max_eta_bin - min_eta_bin + 1
                if max_phi_bin - min_phi_bin + 1 > max_phi_extent:
                    print (max_phi_bin, " to ", min_phi_bin)
                    max_phi_extent = max_phi_bin - min_phi_bin + 1
                    print(max_phi_extent)
                if evt_count_spatial_extent >= max_evt_spatial_extent:
                    break
                evt_count_spatial_extent += 1
            print("Eta fired cells: ", max_eta_extent)
            print("Phi fired cells: ", max_phi_extent)
            print("Layers: ", max_layer + 1)

            #n_eta_bins = max_eta_bin - min_eta_bin
            #n_phi_bins = max_phi_bin - min_phi_bin
            #n_layers = max_layer

            #n_eta_bins = max_eta_bin - min_eta_bin
            #n_phi_bins = max_phi_bin - min_phi_bin
            n_layers = 12
