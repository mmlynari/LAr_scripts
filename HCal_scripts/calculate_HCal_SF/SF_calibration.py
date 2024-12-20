## ABOUT SF_calibration script
## This script performs the HCal sampling fraction (SF) calculation as the ratio of the energy deposited in active material devided by the generated particle energy
## SF is used to bring back the reconstructed energy to the level of generated energy = correct for losses in the passive material
## Additionally, this script calculates the inverse SF that is used by the reconstruction code run_reco_HCal.py
## As an input, you need to simulate several thousands of events, shooting electrons (for the calibration at electromagnetic scale)
## or charged pions (for the calibration at hadron scale) at a given theta, BUT with setting invSF=1 in run_reco_HCal.py (no sampling fraction used for the calibration)
## NOTE: Sampling fraction depends on theta (eta), so every time you want to perform studies at different angle, you need to recalculate the SF

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
import os, sys, glob
import numpy as np
import argparse
from math import sqrt
from datetime import date
#from shutil import copy
from copy import copy
from calibration_functions import get_sampling_fraction
import aux_colors

import gStyle

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")

parser = argparse.ArgumentParser()
# parser.add_argument("-inputFiles", default = "/eos/user/m/mmlynari/FCC_rootfile_storage/HCal_v26May23/HCal_SF_v26May23/230905/output_fullCalo_SimAndDigi_withCluster_MagneticField_False_pMin_100000_MeV_ThetaMinMax_69.805_69.805_pdgId_211_pythiaFalse_NoiseFalse.root", help = "Regex for input files.", type = str)
parser.add_argument("-inputFiles", default = "/afs/cern.ch/user/m/mmlynari/workspace/ALLEGRO_PandoraPFA/HCal_standalone/outputs/241208/ALLEGRO_reco_e.root", help = "Regex for input files.", type = str)
parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d"), help = "Postfix to append to the output folder.", type = str)
parser.add_argument("-calibrate", help = "Calculate HCal scale factors", action='store_true')
parser.add_argument("-check", help = "Check whether obtained SF is correct", action='store_true')
args = parser.parse_args()


if args.check:
	calib = 'check'
else: 
	calib = 'SF'

plot_dir_name = 'calib_'+ calib +'_'+ args.outputPostfix
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

inputFiles = glob.glob(args.inputFiles)
if not inputFiles:
    print("No file found")
    sys.exit(1)

for inputFile in inputFiles:
    print("Treating %s..."%inputFile)
    rootfile_path = inputFile
    ## OUTDATED splitting the input file name to extract some information
    '''
    energy = inputFile.split('_pMin_')[1].split("_")[0]
    energy_gev_float = float(energy)
    energy = str(energy_gev_float).replace(".","dot")
    pdgID = inputFile.split('_pdg_')[1].split("_")[0]
    theta = inputFile.split('_ThetaMinMax_')[1].split("_")[0] 
    theta = str(theta).replace(".","dot") 
    prefix = energy + "GeV_pdgID_" + pdgID + "_theta_" + theta
    ''' 
    prefix = "HCal_SF"

    ## book the histogram for the scale factor 

    if args.check:
        th1_SF_HCal = ROOT.TH1F(prefix + "_" + calib, prefix + "_" + calib, 200, 0., 200.)
    else:
        th1_SF_HCal = ROOT.TH1F(prefix + "_" + calib, prefix + "_" + calib, 100, 0.015, 0.05)

    ## open ROOT file, read branches with variables
    f = ROOT.TFile(rootfile_path)
    events = f.Get("events")
    events.SetBranchStatus("*", 0) # //disable all branches
    events.SetBranchStatus("HCalBarrelReadoutPositioned_energy", 1)
    events.SetBranchStatus("genParticle_*", 1)

    evt = 0
    max_evt = 10000

## loop over all generated events, fill the histogram for each event
    for event in events:
        if evt >= max_evt and not max_evt == -1:
            break
        evt += 1
        total_energy_HCal = 0
        gen_particle_energy = 0
        for cell_energy in event.HCalBarrelReadoutPositioned_energy:
            total_energy_HCal += cell_energy
            print(total_energy_HCal)

        for genParticle_idx in range(len(event.genParticle_status)):
            if event.genParticle_status[genParticle_idx] != 1:
                continue
            gen_particle_energy = event.genParticle_energy[genParticle_idx]
            if args.check:
                th1_SF_HCal.Fill(total_energy_HCal)
            else:
                th1_SF_HCal.Fill(total_energy_HCal/gen_particle_energy) 

        if evt % 1000 == 0:
            print("Event processed: %d"%evt)

    output_rootfile_path = os.path.join(plot_dir_name, prefix +".root")
    output_rootfile = ROOT.TFile(output_rootfile_path, "recreate")
    output_rootfile.cd()

## fit the histogram with get_sampling_fraction function which returns parameters of a gaussian fit
    SF_HCal_fit = get_sampling_fraction(th1_SF_HCal, prefix, plot_dir_name)
    SF_HCal_fitted = SF_HCal_fit.Get().Parameter(1)
    sigma_SF_HCal_fitted = SF_HCal_fit.Get().Parameter(2)
    invSF_HCal_fitted = 1/SF_HCal_fitted
    print("File processed: energy_", prefix)
    if args.check:
        print("\033[31mThis text is red\033[0m")
        print("\033[32mHCal mean energy is ",SF_HCal_fitted, " this was a check\033[0m")
    else: 
        print(f"\033[32mHCal SF is {round(SF_HCal_fitted*100, 2)}% , this value is not to be used in the simulation code.\033[0m")
        print("\033[31mHCal invSF is ",invSF_HCal_fitted, ", put this number into run_reco_HCal.py and rerun the simulation.\033[0m")
    output_rootfile.Close()

