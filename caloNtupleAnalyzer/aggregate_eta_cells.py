import ROOT
import pprint
import argparse
import re
from math import sqrt, atan, tan, exp, pi
import sys
#import numpy as np

ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser = argparse.ArgumentParser()
#parser.add_argument("-inputFile", default = "/eos/user/b/brfranco/rootfile_storage/fcc_analysis_ouput/220810_gamma_flat_1_100_noNoise_caloReco/fccsw_output_pdgID_22_pMin_1000_pMax_100000_thetaMin_50_thetaMax_130.root", help = "Path to the input rootfile. Output rootfile path will be based on it.", type = str)
parser.add_argument("-inputFile", default = "/eos/user/b/brfranco/rootfile_storage/fcc_analysis_ouput/220810_gamma_flat_1_100_noNoise_caloReco/fccsw_output_pdgID_22_pMin_1000_pMax_100000_thetaMin_50_thetaMax_130.root", help = "Path to the input rootfile. Output rootfile path will be based on it.", type = str)
parser.add_argument("-postfix", default = "test", help = "Postfix to append to the output rootfile.", type = str)
parser.add_argument("-startEvt", default = 0, help = "First event to process in the input rootfile (to allow splitting in many jobsi when running on a single rootfile). Event count starts at 0.", type = int)
parser.add_argument("-endEvt", default = 10, help = "Last event to process in the input rootfile (to allow splitting in many jobsi when running on a single rootfile). Event count starts at 0. Set to -1 for all events.", type = int)

args = parser.parse_args()

file_path = args.inputFile
#file_path = "/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/develop_strip_layer/LAr_scripts/caloNtupleAnalyzer/test.root"
#file_path = "old_test.root"
eta_rebinning_factor = 4
output_file_path = file_path.replace(".root", "_withAggregatedCells_" + args.postfix + ".root")
rootfile = ROOT.TFile(file_path)
output_rootfile = ROOT.TFile(output_file_path, "recreate")

events = rootfile.Get('events')
output_tree = ROOT.TTree('events', 'events')
#output_tree = events.CloneTree()
#events.GetListOfClones().Remove(output_tree)

# split the list of branches in the one with Cells info (to be use for the new eta segmentation), cluster info (we have to altered FirstCell and LastCell) and other branches (we copy them as is)
listOfBranches = [events.GetListOfBranches().At(i).GetName() for i in range(len(events.GetListOfBranches()))]
#print(listOfBranches)
cell_branch_families = set()
cell_branch_suffixes = set()
cluster_branches = set()
other_branches = set()
dict_branchName_vectorForFilling = {}
for branch in listOfBranches:
    #print(branch)
    dict_branchName_vectorForFilling[branch] = ROOT.RVec('float')()
    output_tree.Branch(branch, dict_branchName_vectorForFilling[branch])
    if "Cells" in branch:
        cell_branch_families.add(branch.split("_")[0])
        cell_branch_suffixes.add(branch.split("_")[1])
    elif "Cluster" in branch:
        cluster_branches.add(branch)
    else:
        other_branches.add(branch)
#print(cell_branch_families)
#print(cell_branch_suffixes)

# use the input rootfile info to get old eta segmentation from whatever cell collection
events.GetEntry(0)
random_cell_family = next(iter(cell_branch_families))
anEtaBin = getattr(events, random_cell_family + "_etaBin")[0]
anEta = getattr(events, random_cell_family + "_eta")[0]
anotherEtaBin = getattr(events, random_cell_family + "_etaBin")[1]
anotherEta = getattr(events, random_cell_family + "_eta")[1]
tmp_bin = 2
while(anotherEtaBin == anEtaBin):
    anotherEtaBin = getattr(events, random_cell_family + "_etaBin")[tmp_bin]
    anotherEta = getattr(events, random_cell_family + "_eta")[tmp_bin]
    tmp_bin += 1
binDifference = abs(anEtaBin - anotherEtaBin)
etaDifference = abs(anEta - anotherEta)
etaStep = round(etaDifference / binDifference, 6)
new_etaStep = eta_rebinning_factor * etaStep
eta_etaBinEqualZero = round(anEta - anEtaBin * etaStep, 6)
eta_etaBinEqualMax = abs(eta_etaBinEqualZero)
new_eta_etaBinEqualZero = round(((eta_etaBinEqualZero + (eta_rebinning_factor - 1) * etaStep) + eta_etaBinEqualZero) / 2, 6)
etaBinMax = (eta_etaBinEqualMax - eta_etaBinEqualZero) / float(etaStep)
new_etaBinMax = (eta_etaBinEqualMax - eta_etaBinEqualZero) / float(new_etaStep)
print("Old eta step ", etaStep)
print("New eta step ", new_etaStep)
print("Eta corresponding to etaBin == 0 ", eta_etaBinEqualZero)
print("New eta corresponding to etaBin == 0 ", new_eta_etaBinEqualZero)
print("Eta bin max ", etaBinMax)
print("New eta bin max ", new_etaBinMax)
# define the new eta, theta, z for the new cell size
old_eta_bin_groups = []
new_etas = []
new_thetas = []
i = 0
idx_for_new_binning = 0
while i < etaBinMax:
    old_eta_bin_groups.append([i + idx for idx in range(eta_rebinning_factor)])
    new_eta = round(new_eta_etaBinEqualZero + idx_for_new_binning * new_etaStep, 4)
    new_etas.append(new_eta)
    new_thetas.append(2 * atan(exp(-new_eta)))
    i += eta_rebinning_factor
    idx_for_new_binning += 1
#print("new etas:",  new_etas)
#print("new thetas:", new_thetas)
#print(old_eta_bin_groups)

n_evt = 0
for event in events:
    if(n_evt < args.startEvt):
        continue
    if(n_evt > args.endEvt and not args.endEvt == -1):
        break
    if n_evt % 10 == 0:
        print("Treating event %d"%n_evt)
    # deal with branches that we do not touch
    for branch in other_branches:
        dict_branchName_vectorForFilling[branch].clear()
        for idx in range(len(getattr(event, branch))):
            dict_branchName_vectorForFilling[branch].push_back(getattr(event, branch)[idx])

    # deal with cluster branches and store the first/last cell to be modified later
    dict_clusterType_firstCells = {}
    dict_clusterType_lastCells = {}
    dict_clusterType_newfirstCells = {}
    dict_clusterType_newlastCells = {}
    for branch in cluster_branches:
        if "firstCell" in branch:
            dict_clusterType_firstCells[branch] = []
            dict_clusterType_newfirstCells[branch] = []
            for idx in range(len(getattr(event, branch))):
                dict_clusterType_firstCells[branch].append(getattr(event, branch)[idx])
                dict_clusterType_newfirstCells[branch].append(getattr(event, branch)[idx])
        elif "lastCell" in branch:
            dict_clusterType_lastCells[branch] = []
            dict_clusterType_newlastCells[branch] = []
            for idx in range(len(getattr(event, branch))):
                dict_clusterType_lastCells[branch].append(getattr(event, branch)[idx])
                dict_clusterType_newlastCells[branch].append(getattr(event, branch)[idx])

    # deal with cell branches that we have to aggregate in eta everywhere except in strip layer 
    for cell_branch_family in cell_branch_families:
        n_cell_merged = 0
        for cell_branch_suffix in cell_branch_suffixes:
            dict_branchName_vectorForFilling[cell_branch_family + "_" + cell_branch_suffix].clear()
        dict_layer_phiBin_etaGroup_cellIdx = {} # used to see to which coarse eta bin the old cell belong
        layers_in_evt = set(getattr(event, cell_branch_family + "_layer"))
        phiBins_in_evt = set(getattr(event, cell_branch_family + "_phiBin"))
        # prepare the disctionary
        for layer in layers_in_evt:
            dict_layer_phiBin_etaGroup_cellIdx[str(layer)] = {}
            for phiBin in phiBins_in_evt:
                dict_layer_phiBin_etaGroup_cellIdx[str(layer)][str(phiBin)] = {}
                for etaGroup in old_eta_bin_groups:
                    dict_layer_phiBin_etaGroup_cellIdx[str(layer)][str(phiBin)][str(etaGroup)] = -1
        # fill the new cell collection
        for cell_idx in range(len(getattr(event, cell_branch_family + "_layer"))):
            need_to_decrement_cluster_index = False
            if getattr(event, cell_branch_family + "_layer")[cell_idx] == 1:# do not alter the strip layer cells
                for cell_branch_suffix in cell_branch_suffixes:
                    dict_branchName_vectorForFilling[cell_branch_family + "_" + cell_branch_suffix].push_back(getattr(event, cell_branch_family + "_" + cell_branch_suffix)[cell_idx])
            else:
                # find to which coarse eta bin the cell belongs
                current_etaGroup = ""
                for etaGroup in old_eta_bin_groups:
                    if getattr(event, cell_branch_family + "_etaBin")[cell_idx] in etaGroup:
                        current_etaGroup = etaGroup
                if current_etaGroup == "":
                    print("Error: found a cell not belonging to any coarse eta bin, consider enlarging etaBinMax")
                    sys.exit(1)
                # check if we had already one cell in this eta bin and fill things accordingly
                if dict_layer_phiBin_etaGroup_cellIdx[str(getattr(event, cell_branch_family + "_layer")[cell_idx])][str(getattr(event, cell_branch_family + "_phiBin")[cell_idx])][str(current_etaGroup)] == -1:
                    # save the current index in the output vector to be able to add later cell energy to the right eta bini(if size is called before the push_back it gives directly the index
                    dict_layer_phiBin_etaGroup_cellIdx[str(getattr(event, cell_branch_family + "_layer")[cell_idx])][str(getattr(event, cell_branch_family + "_phiBin")[cell_idx])][str(current_etaGroup)] = dict_branchName_vectorForFilling[cell_branch_family + "_layer"].size() # use layer here but any other member would work
                    for cell_branch_suffix in cell_branch_suffixes:
                        if cell_branch_suffix == 'etaBin':
                            dict_branchName_vectorForFilling[cell_branch_family + "_" + cell_branch_suffix].push_back(old_eta_bin_groups.index(current_etaGroup))
                        elif cell_branch_suffix == 'eta':
                            dict_branchName_vectorForFilling[cell_branch_family + "_" + cell_branch_suffix].push_back(new_etas[old_eta_bin_groups.index(current_etaGroup)])
                        elif cell_branch_suffix == 'theta':
                            dict_branchName_vectorForFilling[cell_branch_family + "_" + cell_branch_suffix].push_back(new_thetas[old_eta_bin_groups.index(current_etaGroup)])
                        elif cell_branch_suffix == 'z':
                            radius = sqrt(getattr(event, cell_branch_family + "_x")[cell_idx] * getattr(event, cell_branch_family + "_x")[cell_idx] + getattr(event, cell_branch_family + "_y")[cell_idx] * getattr(event, cell_branch_family + "_y")[cell_idx])
                            dict_branchName_vectorForFilling[cell_branch_family + "_" + cell_branch_suffix].push_back(radius / tan(new_thetas[old_eta_bin_groups.index(current_etaGroup)]))
                        else:
                            dict_branchName_vectorForFilling[cell_branch_family + "_" + cell_branch_suffix].push_back(getattr(event, cell_branch_family + "_" + cell_branch_suffix)[cell_idx])
                else:# if we had already one cell populating this coarse eta bin, just add present cell energy to the existing entry
                    dict_branchName_vectorForFilling[cell_branch_family + "_energy"][dict_layer_phiBin_etaGroup_cellIdx[str(getattr(event, cell_branch_family + "_layer")[cell_idx])][str(getattr(event, cell_branch_family + "_phiBin")[cell_idx])][str(current_etaGroup)]] += getattr(event, cell_branch_family + "_energy")[cell_idx]
                    n_cell_merged += 1
                    if "Cluster" in cell_branch_family: # because it  makes sense only for cell collections that are attached to a cluster
                        need_to_decrement_cluster_index = True
            if need_to_decrement_cluster_index:
                cluster_type = re.search('Positioned(.*)Cells', cell_branch_family).group(1) + "s"
                for clusterIdx in range(len(dict_clusterType_firstCells[cluster_type + "_firstCell"])):
                    tmp_firstCell = dict_clusterType_firstCells[cluster_type + "_firstCell"][clusterIdx]
                    tmp_lastCell = dict_clusterType_lastCells[cluster_type + "_lastCell"][clusterIdx]
                    if tmp_firstCell <= cell_idx < tmp_lastCell:
                        dict_clusterType_newlastCells[cluster_type + "_lastCell"][clusterIdx] -= 1
                        for second_clusterIdx in range(clusterIdx+1, len(dict_clusterType_firstCells[cluster_type + "_firstCell"])):
                            dict_clusterType_newfirstCells[cluster_type + "_firstCell"][second_clusterIdx] -= 1
                            dict_clusterType_newlastCells[cluster_type + "_lastCell"][second_clusterIdx] -= 1
        #if "CaloCluster" in cell_branch_family:
        #    print("Cell to be merged: ", n_cell_merged)
        #end of loop on cells
    #end of loop on cell branch family

    # store cluster info
    for branch in cluster_branches:
        if "Corrected" in branch or not "CaloCluster" in branch:
            continue
        dict_branchName_vectorForFilling[branch].clear()
        if "firstCell" in branch:
            for clusterIdx in range(len(dict_clusterType_firstCells[branch])):
                #print("Old first cell: ", getattr(event, branch)[clusterIdx])
                #print("New first cell: ", dict_clusterType_newfirstCells[branch][clusterIdx])
                dict_branchName_vectorForFilling[branch].push_back(dict_clusterType_newfirstCells[branch][clusterIdx])
        elif "lastCell" in branch:
            for clusterIdx in range(len(dict_clusterType_lastCells[branch])):
                #print("Old last cell: ", getattr(event, branch)[clusterIdx])
                #print("New last cell: ", dict_clusterType_newlastCells[branch][clusterIdx])
                dict_branchName_vectorForFilling[branch].push_back(dict_clusterType_newlastCells[branch][clusterIdx])
        else:
            for idx in range(len(getattr(event, branch))):
                dict_branchName_vectorForFilling[branch].push_back(getattr(event, branch)[idx])
    output_tree.Fill()
    print("##################")
    n_evt += 1
    #end of loop on events
output_rootfile.Write()
output_rootfile.Close()
print(output_file_path, " written.")
rootfile.Close()
