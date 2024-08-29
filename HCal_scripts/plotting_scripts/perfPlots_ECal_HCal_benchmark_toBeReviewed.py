import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
import os, sys, glob
import numpy as np
import argparse
from math import sqrt, floor
from datetime import date
#from shutil import copy
from copy import copy

import gStyle

print("Launch without sourcing k4 environment!!")

#ROOT.gROOT.ProcessLine(".L FCCAnalysesDict.C+")
ROOT.gROOT.ProcessLine(".L /afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/211210/FCCAnalyses/install/lib/libFCCAnalyses.C+")

parser = argparse.ArgumentParser()
#parser.add_argument("-inputFiles", default = "/eos/user/m/mmlynari/FCC_rootfile_storage/ECal_HCal_v24Apr24/240502_energies_10kevt_topoAndSW_noNoise_EMscale/ntuples/240504/fccsw_output_pdgID_211_pMin_*.root", help = "Regex for input files.", type = str)
#parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d") + "_ECal_HCal_thetaSegmentation_pi_EMscale_SW", help = "Postfix to append to the output folder.", type = str)

parser.add_argument("-inputFiles", default = "/eos/user/m/mmlynari/FCC_rootfile_storage/ECal_HCal_v24Apr24/240429_energies_10kevt_topoAndSW_noNoise_3MeV/ntuples/240506/fccsw_output_pdgID_211_pMin_*.root", help = "Regex for input files.", type = str)
parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d") + "_ECal_HCal_thetaSegmentation_pi_benchmark_SW", help = "Postfix to append to the output folder.", type = str)

#parser.add_argument("-inputFiles", default = "/eos/user/m/mmlynari/FCC_rootfile_storage/ECal_HCal_v24Apr24/fcc_analysis_ouput/240508_energies_10kevt_topoAndSW_noNoise_benchmark_0dot01MeV/fccsw_output_pdgID_211_pMin_*.root", help = "Regex for input files.", type = str)
#parser.add_argument("-outputPostfix", default = date.today().strftime("%y%m%d") + "_ECal_HCal_thetaSegmentation_pi_benchmark_0dot01MeV_Topo_highestEneCluster", help = "Postfix to append to the output folder.", type = str)


parser.add_argument("-color", default = 46, help = "Color of the graph", type = int)
parser.add_argument("-markerStyle", default = 21, help = "Style of the graph markers", type = int)
parser.add_argument("-cells", default = True, help = "Also produce plots with total energy from cells (no clustering) -- needs some refreshing", type = bool)
parser.add_argument("-inclEcal", default = True, help = "Include ECal in cells plots", type = bool)
parser.add_argument("-benchmark", default = False, help = "Run benchmark calibration on cells", type = bool)
parser.add_argument("-clusterCollection", default = "CorrectedCaloClusters", help = "Name of the cluster collection to use", type = str)
args = parser.parse_args()

plot_dir_name = 'plots_performances_'+ args.outputPostfix
if not os.path.isdir(plot_dir_name):
    os.mkdir(plot_dir_name)

individual_plots = os.path.join(plot_dir_name, "resolution") 
if not os.path.isdir(individual_plots):
    os.mkdir(individual_plots)

inputFiles = glob.glob(args.inputFiles)
if not inputFiles:
    print("No file found")
    sys.exit(1)

if args.inclEcal:
        calo_label='ECal+HCal'
else:
    calo_label='HCal'

scale = 'benchmark'
if args.benchmark:
    scale_cluster = 'benchmark'
else:
    scale_cluster=scale

# Single resolution histogram 

def draw_resol_canvas(th1, prefix, variable): 
    canvas_resol = ROOT.TCanvas(prefix + variable + "_resolution", prefix + variable + "_resolution")
    fit_range_min = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin()) - 2 * th1.GetRMS()
    fit_range_max = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin()) + 2 * th1.GetRMS()
    fit_result = th1.Fit("gaus", "SQ", "", fit_range_min, fit_range_max)
    th1.Draw()
    if 'rel' in variable:
        name = variable.replace("rel", "")
        th1.GetXaxis().SetTitle("%s_{Reco} - %s_{True}/%s_{Reco}"%(name, name, name))
        th1.GetXaxis().SetTitleOffset(1.2)
        if 'eGenDivided' in variable:
            th1.GetXaxis().SetTitle("E_{Reco} - E_{True}/E_{True}")
    elif 'esponse' in variable:
        th1.GetXaxis().SetTitle("#frac{E_{Reco}}{E_{True}}")
        th1.GetXaxis().SetTitleOffset(1.2)
    elif 'recoE' in variable:
        if 'cells' in variable: 
            th1.GetXaxis().SetTitle("#scale[1.4]{E_{rec}^{cells}}")
        else: 
            th1.GetXaxis().SetTitle("#scale[1.4]{E_{rec}^{cluster}}")
        th1.GetXaxis().SetTitleOffset(1.2)
    else:
        #th1.GetXaxis().SetTitle("#{0}_{Reco} - #{0}_{True}".format(variable))
        th1.GetXaxis().SetTitle("#%s_{Reco} - #%s_{True}"%(variable, variable))

    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptFit(1)

    canvas_resol.Print(os.path.join(individual_plots, prefix + variable + "resol.png"))
    th1.Write()
    return fit_result

max_evt = 10000

# benchmark method parameters obtained from 3000 evts of 100 GeV pion
#[1.29167, 1., 0.915777, -0.00199074],
## for 1000 GeV pion and 2000 events --> [1.14678, 1, 0.803669, -0.000154067]

last_ECal_layer = 11

p0_approx = 1.2909
p1 = 1.
p2_approx = 0.91
p3_approx = -0.0019

a0 = 2.27
b0 = 1.065
a2 = -4.14
b2 = 0.915
a3 = -0.26
b3 = 0.0004

a0_low_ene = 0.34
b0_low_ene = -0.0906
c0_low_ene = 0.0059
d0_low_ene = 2.0184
a2_low_ene = -1.62
b2_low_ene = 0.2694
a3_low_ene = -0.24
b3_low_ene = 0.0044

cutoff_dR = 0.03
cutoff_dPhi = 0.1
cutoff_dTheta = 0.1
cutoff_relE = 0.5

dict_energy_Eresol_error = {}
dict_energy_Eresol_error_clusters_sum = {}
dict_energy_Eresol_error_clusters_all = {}
dict_energy_Eresol_error_clusters_dRcut = {}
dict_energy_Eresol_error_clusters_dRcut_divTrue = {}

dict_energy_Eresol_error_cells = {}
dict_energy_Eresol_error_cells_all = {}
dict_energy_relEresol_error = {}
dict_energy_relEresol_eGenDivided_error = {}
dict_energy_energyResponseResol_error = {}#
dict_energy_efficiency_error = {}
dict_energy_phiResol_error = {}
dict_energy_thetaResol_error = {}
dict_energy_cells_relEresol_error = {}

dict_energy_linearity_error = {}
dict_energy_linearity_error_clusters_sum = {}
dict_energy_linearity_error_clusters_all = {}
dict_energy_linearity_error_clusters_dRcut = {}

dict_energy_linearity_error_cells = {}
dict_energy_linearity_error_cells_all = {}
dict_energy_benchmark_Eresol_error = {}
dict_energy_linearity_benchmark = {}
dict_energy_linearity_benchmark_approx = {}
energies_gev_float = []

for inputFile in inputFiles:
    print("Treating %s..."%inputFile)
    rootfile_path = inputFile
    energy = inputFile.split('_pMin_')[1].split("_")[0]
    energy_gev_float = int(energy)/1000.0
    energies_gev_float.append(energy_gev_float)
    energy = str(energy_gev_float).replace(".","dot")
    print(energy)
    dict_energy_Eresol_error[energy] = []
    dict_energy_Eresol_error_clusters_all[energy] = []
    dict_energy_Eresol_error_clusters_sum[energy] = []
    dict_energy_Eresol_error_clusters_dRcut[energy] = []
    dict_energy_Eresol_error_clusters_dRcut_divTrue[energy] = []
    dict_energy_Eresol_error_cells[energy] = []
    dict_energy_Eresol_error_cells_all[energy] = [] 
    dict_energy_relEresol_error[energy] = []
    dict_energy_relEresol_eGenDivided_error[energy] = []
    dict_energy_energyResponseResol_error[energy] = []
    dict_energy_phiResol_error[energy] = []
    dict_energy_thetaResol_error[energy] = []
    dict_energy_cells_relEresol_error[energy] = []
    dict_energy_efficiency_error[energy] = []
    dict_energy_linearity_error[energy] = []
    dict_energy_linearity_error_clusters_all[energy] = []
    dict_energy_linearity_error_clusters_sum[energy] = []
    dict_energy_linearity_error_clusters_dRcut[energy] = []

    dict_energy_linearity_error_cells[energy] = []
    dict_energy_linearity_error_cells_all[energy] = []
    dict_energy_benchmark_Eresol_error[energy] = []
    dict_energy_linearity_benchmark[energy] = []
    dict_energy_linearity_benchmark_approx[energy] = []

    prefix = energy + "GeV_"

    phiresol_range = 0.1
    if energy_gev_float < 6:
        phiresol_range = 0.01
    if energy_gev_float < 6:
        phiresol_range = 0.03

    th1_phiresol = ROOT.TH1F(prefix + "phi_resolution", prefix + "phi_resolution", 100, -1 * phiresol_range, phiresol_range)
    th1_thetaresol = ROOT.TH1F(prefix + "theta_resolution", prefix + "theta_resolution", 100, -0.5, 0.5)
    th1_angularresol = ROOT.TH1F(prefix + "angular_resolution", prefix + "angular_resolution", 100, 0, 0.03)
    th1_Eresol = ROOT.TH1F(prefix + "energy_resolution", prefix + "energy_resolution", 500, -25, 20)
    th1_Eresol_cells = ROOT.TH1F(prefix + "energy_resolution_cells", prefix + "energy_resolution_cells", 500, -25, 20)
    th1_Ereco = ROOT.TH1F(prefix + "reconstructed_cluster_energy", prefix + "reconstructed_cluster_energy", 200, 0, energy_gev_float*2)
    th1_Ereco_clusters_all = ROOT.TH1F(prefix + "reco_cluster_ene_all", prefix + "reco_cluster_ene_all", 200, 0, energy_gev_float*2)
    th1_Ereco_clusters_sum = ROOT.TH1F(prefix + "reco_cluster_ene_sum_per_event", prefix + "reco_cluster_ene_sum_per_event", 200, 0, energy_gev_float*2)
    th1_Ereco_clusters_dRcut = ROOT.TH1F(prefix + "reco_cluster_ene_sum_in_dR", prefix + "reco_cluster_ene_sum_in_dR", 200, 0, energy_gev_float*2)

    th1_Ereco_cells = ROOT.TH1F(prefix + "reconstructed_cells_energy", prefix + "reconstructed_cells_energy", 200, 0, energy_gev_float*2)
    th1_Ereco_cells_all = ROOT.TH1F(prefix + "all_reconstructed_cells_energy", prefix + "all_reconstructed_cells_energy", 200, 0, energy_gev_float*2)
    th1_Ereco_cells_ECal = ROOT.TH1F(prefix + "ECal_reconstructed_cells_energy", prefix + "ECal_reconstructed_cells_energy", 200, 0, energy_gev_float*2)
    th1_Ereco_cells_HCal = ROOT.TH1F(prefix + "HCal_reconstructed_cells_energy", prefix + "HCal_reconstructed_cells_energy", 200, 0, energy_gev_float*2)

    th1_relEresol = ROOT.TH1F(prefix + "relative_energy_resolution", prefix + "relative_energy_resolution", 100, -0.6, 0.6)
    th1_relEresol_eGenDivided = ROOT.TH1F(prefix + "relative_energy_resolution_toEGen", prefix + "relative_energy_resolution_toEgen", 150, -0.6, 0.9)
    th1_energy_response = ROOT.TH1F(prefix + "energy_response", prefix + "energy_response", 100, 0, 2)
    th1_relEresol_cells = ROOT.TH1F(prefix + "relative_energy_resolution_cells", prefix + "relative_energy_resolution_cells", 100, -0.6, 0.6)
    th1_phi_cells = ROOT.TH1F(prefix + "phi_cells", prefix + "phi_cells", 20, -0.25, 0.25)
    th1_eta_cells = ROOT.TH1F(prefix + "eta_cells", prefix + "eta_cells", 20, 0.1, 0.6)
    th2_eta_phi_cells = ROOT.TH2F(prefix + "eta_phi_cells", prefix + "eta_phi_cells", 20, 0.1, 0.6, 20, -0.25, 0.25)
    th1_Ene_benchmark = ROOT.TH1F(prefix + "energy_benchmark", prefix + "energy_benchmark", 200, 0, energy_gev_float*2)
    th1_Ene_benchmark_approx = ROOT.TH1F(prefix + "energy_benchmark_approx", prefix + "energy_benchmark_approx", 200, 0, energy_gev_float*2)

    f = ROOT.TFile(rootfile_path)
    events = f.Get("events")
    events.SetBranchStatus("*", 0) # //disable all branches
    events.SetBranchStatus("genParticle_*", 1)
    if args.cells:
        events.SetBranchStatus("HCalBarrelPositionedCells2_energy", 1)
        events.SetBranchStatus("HCalBarrelPositionedCells2_eta", 1)
        events.SetBranchStatus("HCalBarrelPositionedCells2_phi", 1)
        if args.inclEcal:
            events.SetBranchStatus("ECalBarrelPositionedCells2_energy", 1)
        print("set cell branch")

    if args.benchmark:
        events.SetBranchStatus("HCalBarrelPositionedCells2_x", 1)
        events.SetBranchStatus("HCalBarrelPositionedCells2_y", 1)
        events.SetBranchStatus("ECalBarrelPositionedCells2_layer", 1)
        #args.clusterCollection="CaloTopoClusters"
       
    events.SetBranchStatus(args.clusterCollection + "_*", 1)

    theta_cells = []
    evt = 0
    n_gen_particles = 0
    evt_with_cluster_matching_genParticle = 0
    for event in events:
        if evt >= max_evt and not max_evt == -1:
            break
        evt += 1
        total_energy = 0
        total_energy_HCal = 0
        energy_HCal_first = 0
        total_energy_ECal = 0
        energy_ECal_last = 0
        total_energy_benchmark = 0
        gen_particle_energy = 0

        if args.cells:
            for cell in range(len(getattr(event, "HCalBarrelPositionedCells2_energy"))):
                cell_energy = event.HCalBarrelPositionedCells2_energy[cell]
                total_energy += cell_energy
                eta_cell = -99.
                phi_cell = -99.
                eta_cell = event.HCalBarrelPositionedCells2_eta[cell]
                phi_cell = event.HCalBarrelPositionedCells2_phi[cell]
                #th2_eta_phi_cells.SetBinContent(eta_bin, phi_bin, cell_energy)
                th2_eta_phi_cells.Fill(eta_cell,phi_cell)
                total_energy_HCal += cell_energy

            if args.inclEcal:
                for ECal_cell in range(len(getattr(event, "ECalBarrelPositionedCells2_energy"))):
                    ECal_cell_energy = event.ECalBarrelPositionedCells2_energy[ECal_cell]
                    total_energy += ECal_cell_energy
                    total_energy_ECal += ECal_cell_energy

        if args.benchmark:
            x=-99
            y=-99

            for i_HCal in range(len(event.HCalBarrelPositionedCells2_energy)):
                total_energy_HCal += event.HCalBarrelPositionedCells2_energy[i_HCal]
                ## do not have layer info in the ntuples for HCal, calculate radius instead
                x=event.HCalBarrelPositionedCells2_x[i_HCal]
                y=event.HCalBarrelPositionedCells2_y[i_HCal]
                if ((sqrt(x*x+y*y)>2810.5) and (sqrt(x*x+y*y)<2860.5)):
                    energy_HCal_first += event.HCalBarrelPositionedCells2_energy[i_HCal]

            for i_ECal in range(len(event.ECalBarrelPositionedCells2_energy)):
                total_energy_ECal += event.ECalBarrelPositionedCells2_energy[i_ECal]
                if (event.ECalBarrelPositionedCells2_layer[i_ECal]==last_ECal_layer):
                    energy_ECal_last += event.ECalBarrelPositionedCells2_energy[i_ECal]

            ## first estimate of the benchmark energy with average parameters
            total_energy_benchmark_approx = total_energy_ECal*p0_approx + total_energy_HCal*p1 + p2_approx*sqrt(abs(energy_ECal_last*p0_approx*energy_HCal_first*p1)) + p3_approx*pow(total_energy_ECal*p0_approx,2)
            #if total_energy_benchmark_approx<7.:
            #    p0 = a0_low_ene*total_energy_benchmark_approx + b0_low_ene*total_energy_benchmark_approx**2 + c0_low_ene*total_energy_benchmark_approx**3 + d0_low_ene
            #    p2 = a2_low_ene/total_energy_benchmark_approx + b2_low_ene
            #    p3 = a3_low_ene/total_energy_benchmark_approx + b3_low_ene
            #else: 
            p0 = a0/sqrt(total_energy_benchmark_approx) + b0 
            p2 = a2/total_energy_benchmark_approx + b2
            p3 = a3/total_energy_benchmark_approx + b3

            total_energy_benchmark = total_energy_ECal*p0 + total_energy_HCal*p1 + p2*sqrt(abs(energy_ECal_last*p0*energy_HCal_first*p1)) + p3*pow(total_energy_ECal*p0,2)
            
            
        for genParticle_idx in range(len(event.genParticle_status)):
            if event.genParticle_status[genParticle_idx] != 1:
                continue
            n_gen_particles += 1
            best_genParticle_match_idx = -1
            best_d_R = 100000
            best_d_theta = 100000
            best_d_phi = 100000
            best_d_E = 100000
            best_d_relE = 10000
            best_d_relE_eGenDivided = 10000
            best_energyResponse = 10000
            cluster_energy_sum = 0
            cluster_energy_dRcut = 0 
            highest_cluster_energy = 0


            gen_particle_energy = event.genParticle_energy[genParticle_idx]

            for caloCluster_idx in range(len(getattr(event, args.clusterCollection + "_energy"))):
                ##print("mici cluster idx : ", caloCluster_idx) 
                d_theta = getattr(event, args.clusterCollection + "_theta")[caloCluster_idx] - event.genParticle_theta[genParticle_idx]
                d_phi = getattr(event, args.clusterCollection + "_phi")[caloCluster_idx] - event.genParticle_phi[genParticle_idx]
                d_R = sqrt(d_theta * d_theta + d_phi * d_phi)
                d_E = getattr(event, args.clusterCollection + "_energy")[caloCluster_idx] - event.genParticle_energy[genParticle_idx]
                d_relE = d_E/getattr(event, args.clusterCollection + "_energy")[caloCluster_idx]
                d_relE_eGenDivided = d_E/event.genParticle_energy[genParticle_idx]
                energy_response =  getattr(event, args.clusterCollection + "_energy")[caloCluster_idx] / event.genParticle_energy[genParticle_idx]
                #d_relE = (event.CorrectedCaloClusters_energy[caloCluster_idx] - event.genParticle_energy[genParticle_idx])/event.genParticle_energy[genParticle_idx]
                cluster_energy_all = getattr(event, args.clusterCollection + "_energy")[caloCluster_idx]
                cluster_energy_sum += getattr(event, args.clusterCollection + "_energy")[caloCluster_idx]
                cluster_energy = getattr(event, args.clusterCollection + "_energy")[caloCluster_idx]

                if cluster_energy > highest_cluster_energy:
                    highest_cluster_energy = cluster_energy


                ##print("d_R : ", d_R)
                ##print("best_d_R : ", best_d_R)

                if d_R < best_d_R:
                    best_genParticle_match_idx = caloCluster_idx
                    best_d_R = d_R
                    best_d_theta = d_theta
                    best_d_phi = d_phi
                    best_d_E = d_E
                    best_d_relE = d_relE
                    best_d_relE_eGenDivided = d_relE_eGenDivided
                    best_energyResponse = energy_response
                    ##cluster_energy = getattr(event, args.clusterCollection + "_energy")[caloCluster_idx]

                if d_R < 0.4:
                    cluster_energy_dRcut += getattr(event, args.clusterCollection + "_energy")[caloCluster_idx]

            # for energy resol, ask for phi and theta closeness and vice versa
            if best_genParticle_match_idx == -1:
                continue
            if abs(best_d_relE) < cutoff_relE:
                th1_phiresol.Fill(best_d_phi)
                th1_thetaresol.Fill(best_d_theta)
                th1_angularresol.Fill(best_d_R)

            if args.cells:
                th1_Ereco_cells_all.Fill(total_energy)
                th1_Ereco_cells_ECal.Fill(total_energy_ECal)
                th1_Ereco_cells_HCal.Fill(total_energy_HCal)

            th1_Ereco_clusters_all.Fill(cluster_energy_all)
            th1_Ereco_clusters_sum.Fill(cluster_energy_sum)
            th1_Ereco_clusters_dRcut.Fill(cluster_energy_dRcut)
            th1_Ereco.Fill(highest_cluster_energy)


                
            if best_d_R < cutoff_dR:
                th1_Eresol.Fill(best_d_E)
                th1_relEresol.Fill(best_d_relE)
                th1_relEresol_eGenDivided.Fill(best_d_relE_eGenDivided)
                th1_energy_response.Fill(best_energyResponse)
                ## th1_Ereco.Fill(cluster_energy)

                if args.cells:
                    th1_Eresol_cells.Fill((total_energy - gen_particle_energy))
                    th1_relEresol_cells.Fill((total_energy - gen_particle_energy)/total_energy)
                    th1_Ereco_cells.Fill(total_energy)
                    #print(total_energy)

                if args.benchmark:
                    th1_Ene_benchmark.Fill(total_energy_benchmark)
                    th1_Ene_benchmark_approx.Fill(total_energy_benchmark_approx) 


            # Ask both for dR and dE matching for the efficiency
            if best_d_R < cutoff_dR and abs(best_d_relE) < cutoff_relE:
                evt_with_cluster_matching_genParticle += 1


        if evt % 1000 == 0:
            print("Event processed: %d"%evt)

    output_rootfile_path = os.path.join(individual_plots, prefix + "perfHistograms.root")
    output_rootfile = ROOT.TFile(output_rootfile_path, "recreate")
    output_rootfile.cd()

    print('mici si supis zas')

    phiResolFit = draw_resol_canvas(th1_phiresol, prefix, 'Phi')
    thetaResolFit = draw_resol_canvas(th1_thetaresol, prefix, 'Theta')
    EresolFit = draw_resol_canvas(th1_Eresol, prefix, 'E')
    relEresolFit = draw_resol_canvas(th1_relEresol, prefix, 'relE')
    relEresolFit_eGenDivided = draw_resol_canvas(th1_relEresol_eGenDivided, prefix, 'relE_eGenDivided')
    eResponseFit = draw_resol_canvas(th1_energy_response, prefix, 'energyResponse')
    recoEneFit = draw_resol_canvas(th1_Ereco, prefix, 'recoE')
    recoEneFit_clusters_sum = draw_resol_canvas(th1_Ereco_clusters_sum, prefix, 'recoEclustersSum')
    recoEneFit_clusters_all = draw_resol_canvas(th1_Ereco_clusters_all, prefix, 'recoEclustersAll')
    recoEneFit_clusters_dRcut = draw_resol_canvas(th1_Ereco_clusters_dRcut, prefix, 'recoEclustersdRcut')

    if args.cells:
        relEresolFit_cells = draw_resol_canvas(th1_relEresol_cells, prefix, 'relE')
        dict_energy_cells_relEresol_error[energy].append(relEresolFit_cells.Get().Parameter(2))
        dict_energy_cells_relEresol_error[energy].append(relEresolFit_cells.Get().ParError(2))
        recoEneFit_cells = draw_resol_canvas(th1_Ereco_cells, prefix, 'recoEcells')
        recoEneFit_cells_all = draw_resol_canvas(th1_Ereco_cells_all, prefix, 'recoEcellsall')
        recoEneFit_cells_ECal = draw_resol_canvas(th1_Ereco_cells_ECal, prefix, 'recoEcellsECal')
        recoEneFit_cells_HCal = draw_resol_canvas(th1_Ereco_cells_HCal, prefix, 'recoEcellsHCal')

        ## for the linearity, mean energy divided by the true energy
        dict_energy_linearity_error_cells[energy].append(recoEneFit_cells.Get().Parameter(1)/gen_particle_energy)
        dict_energy_linearity_error_cells_all[energy].append(recoEneFit_cells_all.Get().Parameter(1)/gen_particle_energy)

        ## for the resolution, fitted sigma divided by the mean energy
        dict_energy_Eresol_error_cells[energy].append(recoEneFit_cells.Get().Parameter(2)/recoEneFit_cells.Get().Parameter(1))
        dict_energy_Eresol_error_cells_all[energy].append(recoEneFit_cells_all.Get().Parameter(2)/recoEneFit_cells_all.Get().Parameter(1))


    if args.benchmark:
        Ene_benchmark_fit = draw_resol_canvas(th1_Ene_benchmark, prefix, 'recoEbenchmarkcells')
        Ene_benchmark_fit_approx = draw_resol_canvas(th1_Ene_benchmark_approx, prefix, 'recoEbenchmarkapproxcells')
        dict_energy_benchmark_Eresol_error[energy].append(Ene_benchmark_fit.Get().Parameter(2)/Ene_benchmark_fit.Get().Parameter(1))
        dict_energy_linearity_benchmark[energy].append(Ene_benchmark_fit.Get().Parameter(1)/gen_particle_energy)
        dict_energy_linearity_benchmark_approx[energy].append(Ene_benchmark_fit_approx.Get().Parameter(1)/gen_particle_energy)

    dict_energy_efficiency_error[energy].append(evt_with_cluster_matching_genParticle * 100 / float(n_gen_particles))
    dict_energy_efficiency_error[energy].append(sqrt(evt_with_cluster_matching_genParticle) * 100 / float(n_gen_particles))

    dict_energy_relEresol_error[energy].append(relEresolFit.Get().Parameter(2))
    dict_energy_relEresol_error[energy].append(relEresolFit.Get().ParError(2))

    dict_energy_relEresol_eGenDivided_error[energy].append(relEresolFit_eGenDivided.Get().Parameter(2))
    dict_energy_relEresol_eGenDivided_error[energy].append(relEresolFit_eGenDivided.Get().ParError(2))

    dict_energy_phiResol_error[energy].append(phiResolFit.Get().Parameter(2))
    dict_energy_phiResol_error[energy].append(phiResolFit.Get().ParError(2))

    dict_energy_thetaResol_error[energy].append(thetaResolFit.Get().Parameter(2))
    dict_energy_thetaResol_error[energy].append(thetaResolFit.Get().ParError(2))

    ### for the linearity, mean energy divided by the true energy
    dict_energy_linearity_error[energy].append(recoEneFit.Get().Parameter(1)/gen_particle_energy)
    dict_energy_linearity_error_clusters_sum[energy].append(recoEneFit_clusters_sum.Get().Parameter(1)/gen_particle_energy)
    dict_energy_linearity_error_clusters_all[energy].append(recoEneFit_clusters_all.Get().Parameter(1)/gen_particle_energy)
    dict_energy_linearity_error_clusters_dRcut[energy].append(recoEneFit_clusters_dRcut.Get().Parameter(1)/gen_particle_energy)


    ### for the resolution, fitted sigma divided by the mean energy
    dict_energy_Eresol_error[energy].append(recoEneFit.Get().Parameter(2)/recoEneFit.Get().Parameter(1))
    dict_energy_Eresol_error_clusters_sum[energy].append(recoEneFit_clusters_sum.Get().Parameter(2)/recoEneFit_clusters_sum.Get().Parameter(1))
    dict_energy_Eresol_error_clusters_all[energy].append(recoEneFit_clusters_all.Get().Parameter(2)/recoEneFit_clusters_all.Get().Parameter(1))
    dict_energy_Eresol_error_clusters_dRcut[energy].append(recoEneFit_clusters_dRcut.Get().Parameter(2)/recoEneFit_clusters_dRcut.Get().Parameter(1))
    dict_energy_Eresol_error_clusters_dRcut_divTrue[energy].append(recoEneFit_clusters_dRcut.Get().Parameter(2)/gen_particle_energy)



    output_rootfile.Close()

postfix = ""

def plot_resolution_vs_energy_graph(variable_name, postfix, relEresol_vs_energy_graph, write_formula = False):
    setGridx = False
    setGridy = False
    if 'relEresol' in variable_name:
        plot_title = "#%s energy resolution"%(calo_label)
        y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco} - E_{True}}{E_{Reco}}#right)"
        if 'eGenDivided' in variable_name:
            y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco} - E_{True}}{E_{True}}#right)"
    elif 'esponse' in variable_name:
        plot_title = '#%s energy response'%(calo_label)
        y_axis_label = "#scale[1.9]{#sigma}#left(#frac{E_{Reco}}{E_{True}}#right)"
    elif variable_name == 'phiResol':
        plot_title = '#%s #Phi resolution'%(calo_label)
        y_axis_label = "#scale[1.9]{#sigma}#left(#Phi_{Reco} - #Phi_{True}#right) [mrad]"
    elif variable_name == 'thetaResol':
        plot_title = '#%s #theta resolution'%(calo_label)
        y_axis_label = "#scale[1.9]{#sigma}#left(#theta_{Reco} - #theta_{True}#right) [mrad]"
    elif variable_name == "efficiency":
        plot_title = '#%s efficiency'%(calo_label)
        y_axis_label = "Efficiency [%]"
        setGrid = True
    elif 'linearity' in variable_name:
        plot_title = ' '
        y_axis_label = "#scale[1.4]{#LTE_{rec}#GT/E_{true}}"
        setGridx = False
        setGridy = True
    elif 'Eresol' in variable_name: 
        if 'EresolDivTrue' in variable_name: 
            plot_title = "#%s energy resolutiin"%(calo_label)
            y_axis_label = "#scale[1.4]{#sigma_{E_{Reco}}/E_{True}}"
        else:
            plot_title = "#%s energy resolutiin"%(calo_label)
            y_axis_label = "#scale[1.4]{#sigma_{E_{Reco}}/#LTE_{Reco}#GT}"
    else: 
        plot_title = 'undefined '
        y_axis_label = "undefined "
    print(y_axis_label)

    x_axis_label = "E_{True} [GeV]"
    #relEresol_vs_energy_graph.SetMarkerSize(1.5)
    relEresol_vs_energy_graph.SetMarkerStyle(args.markerStyle)
    relEresol_vs_energy_graph.SetMarkerColor(args.color)
    relEresol_vs_energy_graph.SetTitle(plot_title)
    relEresol_vs_energy_graph.GetXaxis().SetTitle(x_axis_label)
    relEresol_vs_energy_graph.GetXaxis().SetTitleOffset(1.2)
    relEresol_vs_energy_graph.GetXaxis().SetLimits(1, 500)
    #relEresol_vs_energy_graph.GetXaxis().SetLimits(0.7, 300)
    relEresol_vs_energy_graph.GetYaxis().SetTitle(y_axis_label)
    relEresol_vs_energy_graph.GetYaxis().SetTitleOffset(1.45)
    #relEresol_vs_energy_graph = ROOT.TGraphErrors(len(energies_gev_float), energies_gev_float, energy_resolutions, 0, energy_resolution_errors)
    #relEresol_vs_energy_graph.GetYaxis().SetLabelSize(0.06)
    #relEresol_vs_energy_graph.GetYaxis().SetTitleSize(0.06)
    #relEresol_vs_energy_graph.GetYaxis().SetTitleSize(0.05)
    tf1_relEresol_vs_e = ROOT.TF1("tf1_" + variable_name, "sqrt(pow([0]/x, 2) + pow([1]/sqrt(x), 2) + pow([2], 2))", energies_gev_float[0], energies_gev_float[-1])
    tf1_relEresol_vs_e.SetLineColor(args.color)
    y_min=0.9
    y_max=1.1
    if  'linearity' in variable_name:
        if variable_name =='linearitycells':
            relEresol_vs_energy_graph.SetMinimum(y_min)
            relEresol_vs_energy_graph.SetMaximum(y_max)
        else: 
            if args.benchmark:
                relEresol_vs_energy_graph.SetMinimum(y_min)
                relEresol_vs_energy_graph.SetMaximum(y_max)
            else: 
                relEresol_vs_energy_graph.SetMinimum(y_min)
                relEresol_vs_energy_graph.SetMaximum(y_max)
    else:
        relEresol_vs_energy_graph.SetMinimum(0)

    canvas_relEresol_vs_energy = ROOT.TCanvas(variable_name + "_vs_energy", variable_name + "_vs_energy")
    canvas_relEresol_vs_energy.SetLogx(1)
    canvas_relEresol_vs_energy.SetGrid(setGridx, setGridy)
    print('mici si supis')

    ## LEGENDary
    if variable_name == "efficiency":
        eResol_vs_e_legend = copy(gStyle.bottomRight_legend)
    else: 
        eResol_vs_e_legend = copy(gStyle.topRight_legend_relEresol)
    eResol_vs_e_legend.AddEntry(ROOT.nullptr,"#bf{FCC-ee simulation}","")
    eResol_vs_e_legend.AddEntry(ROOT.nullptr, "#pi^{#minus} @ #theta=69, B=0T", "")
    eResol_vs_e_legend.AddEntry(ROOT.nullptr,  "%s Barrel"%(calo_label), "")

    if 'cells' in variable_name:
        if 'benchmark' in variable_name:
            eResol_vs_e_legend.AddEntry(ROOT.nullptr, "cells, benchmark", "")
        else:
            eResol_vs_e_legend.AddEntry(ROOT.nullptr, "cells, %s"%(scale), "")
    else:
        eResol_vs_e_legend.AddEntry(ROOT.nullptr, "Topo clusters, %s"%(scale_cluster), "")

    if write_formula:
        relEresol_vs_e_fit = relEresol_vs_energy_graph.Fit(tf1_relEresol_vs_e, "SM")
        print(relEresol_vs_e_fit.Get().Parameter(0))
        #a = str(round(relEresol_vs_e_fit.Get().Parameter(0), 4))
        b = str(round(relEresol_vs_e_fit.Get().Parameter(1)*100, 1))
        c = str(round(relEresol_vs_e_fit.Get().Parameter(2)*100, 1))
        #eResol_vs_e_formula = "#color[%d]{#frac{%s}{E} #oplus #frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, a, b, c) 
        eResol_vs_e_formula = "#color[%d]{#frac{%s}{#sqrt{E}} #oplus %s}"%(args.color, b, c) 
        #eResol_vs_e_legend = copy(gStyle.topRight_legend_relEresol)
        eResol_vs_e_legend.AddEntry(ROOT.nullptr, "", "")
        eResol_vs_e_legend.AddEntry(ROOT.nullptr, eResol_vs_e_formula, "")
    relEresol_vs_energy_graph.Draw("ap")
    eResol_vs_e_legend.Draw()
    relEresol_vs_energy_graph.Write()
    canvas_relEresol_vs_energy.Print(os.path.join(plot_dir_name, variable_name + "_vs_energy.png"))
    canvas_relEresol_vs_energy.Write()
    #return relEresol_vs_e_fit


energies_gev_float.sort()
idx = 0
relEresol_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_relEresol_vs_energy_graph")
relEresol_eGenDivided_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_relEresol_eGenDivided_vs_energy_graph")
relEresol_vs_energy_cells_graph = ROOT.TGraphErrors(prefix + postfix)
#relEresol_vs_energy_cells_graph, relEerol_vs_energy_fit = create_resolution_vs_energy_graph('relEresol', "cells")
Eresol_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_Eresol_vs_energy_graph")
Eresol_vs_energy_graph_clusters_sum = ROOT.TGraphErrors(prefix + postfix + "_Eresol_vs_energy_graph_clusters_sum")
Eresol_vs_energy_graph_clusters_all = ROOT.TGraphErrors(prefix + postfix + "_Eresol_vs_energy_graph_clusters_all")
Eresol_vs_energy_graph_clusters_dRcut = ROOT.TGraphErrors(prefix + postfix + "_Eresol_vs_energy_graph_clusters_dRcut")
Eresol_vs_energy_graph_clusters_dRcut_divTrue = ROOT.TGraphErrors(prefix + postfix + "_EresolDivTrue_vs_energy_graph_clusters_dRcut")



Eresol_vs_energy_graph_cells = ROOT.TGraphErrors(prefix + postfix + "_Eresol_vs_energy_graph_cells")
Eresol_vs_energy_graph_cells_all = ROOT.TGraphErrors(prefix + postfix + "_Eresol_vs_energy_graph_cells_all")
Eresol_vs_energy_benchmark_graph = ROOT.TGraphErrors(prefix + postfix + "resolution_vs_energy_benchmark")

efficiency_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_efficiency_vs_energy_graph")

phiResol_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_phiResol_vs_energy_graph")
thetaResol_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "_thetaResol_vs_energy_graph")

linearity_vs_energy_graph = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_graph")
linearity_vs_energy_graph_clusters_all = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_graph_clusters_all")
linearity_vs_energy_graph_clusters_sum = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_graph_clusters_sum")
linearity_vs_energy_graph_clusters_dRcut = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_graph_clusters_dRcut")


linearity_vs_energy_graph_cells = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_graph_cells")
linearity_vs_energy_benchmark_graph = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_benchmark")
linearity_vs_energy_benchmark_approx_graph = ROOT.TGraphErrors(prefix + postfix + "linearity_vs_energy_benchmark_approx")


for energy_float in energies_gev_float:
    energy_str = str(energy_float).replace(".", "dot")

    # E resol
    relEresol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_relEresol_error[energy_str][0])
    relEresol_vs_energy_graph.SetPointError(idx, 0, dict_energy_relEresol_error[energy_str][1])
    relEresol_eGenDivided_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_relEresol_eGenDivided_error[energy_str][0])
    relEresol_eGenDivided_vs_energy_graph.SetPointError(idx, 0, dict_energy_relEresol_eGenDivided_error[energy_str][1])
    if args.cells:
        relEresol_vs_energy_cells_graph.SetPoint(idx, energy_float, dict_energy_cells_relEresol_error[energy_str][0])
        relEresol_vs_energy_cells_graph.SetPointError(idx, 0, dict_energy_cells_relEresol_error[energy_str][1])
        linearity_vs_energy_graph_cells.SetPoint(idx, energy_float, dict_energy_linearity_error_cells[energy_str][0])

    if args.benchmark: 
        Eresol_vs_energy_benchmark_graph.SetPoint(idx, energy_float, dict_energy_benchmark_Eresol_error[energy_str][0])
        linearity_vs_energy_benchmark_graph.SetPoint(idx, energy_float, dict_energy_linearity_benchmark[energy_str][0])
        linearity_vs_energy_benchmark_approx_graph.SetPoint(idx, energy_float, dict_energy_linearity_benchmark_approx[energy_str][0])

    Eresol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_Eresol_error[energy_str][0])
    Eresol_vs_energy_graph_clusters_sum.SetPoint(idx, energy_float, dict_energy_Eresol_error_clusters_sum[energy_str][0])
    Eresol_vs_energy_graph_clusters_all.SetPoint(idx, energy_float, dict_energy_Eresol_error_clusters_all[energy_str][0])
    Eresol_vs_energy_graph_clusters_dRcut.SetPoint(idx, energy_float, dict_energy_Eresol_error_clusters_dRcut[energy_str][0])
    Eresol_vs_energy_graph_clusters_dRcut_divTrue.SetPoint(idx, energy_float, dict_energy_Eresol_error_clusters_dRcut_divTrue[energy_str][0])

    Eresol_vs_energy_graph_cells.SetPoint(idx, energy_float, dict_energy_Eresol_error_cells[energy_str][0])
    Eresol_vs_energy_graph_cells_all.SetPoint(idx, energy_float, dict_energy_Eresol_error_cells_all[energy_str][0])

    # Phi resol
    phiResol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_phiResol_error[energy_str][0]*1000)
    phiResol_vs_energy_graph.SetPointError(idx, 0, dict_energy_phiResol_error[energy_str][1]*1000)

    # Theta resol
    thetaResol_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_thetaResol_error[energy_str][0]*1000)
    thetaResol_vs_energy_graph.SetPointError(idx, 0, dict_energy_thetaResol_error[energy_str][1]*1000)

    #Efficiency
    efficiency_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_efficiency_error[energy_str][0])
    efficiency_vs_energy_graph.SetPointError(idx, 0, dict_energy_efficiency_error[energy_str][1])

    #Linearity
    linearity_vs_energy_graph.SetPoint(idx, energy_float, dict_energy_linearity_error[energy_str][0])
    linearity_vs_energy_graph_clusters_all.SetPoint(idx, energy_float, dict_energy_linearity_error_clusters_all[energy_str][0])
    linearity_vs_energy_graph_clusters_sum.SetPoint(idx, energy_float, dict_energy_linearity_error_clusters_sum[energy_str][0])
    linearity_vs_energy_graph_clusters_dRcut.SetPoint(idx, energy_float, dict_energy_linearity_error_clusters_dRcut[energy_str][0])

    idx +=1

# Writing the VARIABLE_RESOL_VS_ENERGY
ROOT.gStyle.SetOptFit(000)

output_rootfile_graph_name = os.path.join(individual_plots, "relResol_vs_energy.root")
print("Writing %s"%output_rootfile_graph_name)
output_rootfile_graph = ROOT.TFile(output_rootfile_graph_name, "recreate")

Eresol_vs_energy_fit = plot_resolution_vs_energy_graph('Eresol', "", Eresol_vs_energy_graph, write_formula = True)
Eresol_vs_energy_fit_clusters_sum = plot_resolution_vs_energy_graph('Eresolclusterssum', "", Eresol_vs_energy_graph_clusters_sum, write_formula = True)
Eresol_vs_energy_fit_clusters_all = plot_resolution_vs_energy_graph('Eresolclustersall', "", Eresol_vs_energy_graph_clusters_all, write_formula = True)
Eresol_vs_energy_fit_clusters_dRcut = plot_resolution_vs_energy_graph('EresolclustersdRcut', "", Eresol_vs_energy_graph_clusters_dRcut, write_formula = True)
Eresol_vs_energy_fit_clusters_dRcut_divTrue = plot_resolution_vs_energy_graph('EresolDivTrueclustersdRcut', "", Eresol_vs_energy_graph_clusters_dRcut_divTrue, write_formula = True)




##relEresol_vs_energy_fit = plot_resolution_vs_energy_graph('relEresol', "", relEresol_vs_energy_graph, write_formula = True)
##relEresol_eGenDivided_vs_energy_fit = plot_resolution_vs_energy_graph('relEresol_eGenDivided', "", relEresol_eGenDivided_vs_energy_graph, write_formula = True)

if args.cells:
    #relEresol_cell_vs_energy_fit = plot_resolution_vs_energy_graph('relEresolCell', "", relEresol_vs_energy_cells_graph, write_formula = True)
    Eresol_vs_energy_fit_cells = plot_resolution_vs_energy_graph('Eresolcells', "", Eresol_vs_energy_graph_cells, write_formula = True)
    Eresol_vs_energy_fit_cells_all = plot_resolution_vs_energy_graph('Eresolcellsall', "", Eresol_vs_energy_graph_cells_all, write_formula = True)
    linearity_vs_energy_fit_cells = plot_resolution_vs_energy_graph('linearitycells', "", linearity_vs_energy_graph_cells)

if args.benchmark:
    Eresol_benchmark_vs_energy_fit = plot_resolution_vs_energy_graph('Eresolbenchmarkcells', "", Eresol_vs_energy_benchmark_graph, write_formula = True)
    linearity_benchmark_vs_energy_fit = plot_resolution_vs_energy_graph('linearitybenchmarkcells', "", linearity_vs_energy_benchmark_graph)
    linearity_benchmark_approx_vs_energy_fit = plot_resolution_vs_energy_graph('linearitybenchmarkapproxcells', "", linearity_vs_energy_benchmark_approx_graph)



#phiResol_vs_energy_fit = plot_resolution_vs_energy_graph('phiResol', "", phiResol_vs_energy_graph)
#thetaResol_vs_energy_fit = plot_resolution_vs_energy_graph('thetaResol', "", thetaResol_vs_energy_graph)
efficiency_vs_energy_fit = plot_resolution_vs_energy_graph('efficiency', "", efficiency_vs_energy_graph)

linearity_vs_energy_fit = plot_resolution_vs_energy_graph('linearity', "", linearity_vs_energy_graph)
linearity_vs_energy_fit_clusters_all = plot_resolution_vs_energy_graph('linearityclustersall', "", linearity_vs_energy_graph_clusters_all)
linearity_vs_energy_fit_clusters_sum = plot_resolution_vs_energy_graph('linearityclusterssum', "", linearity_vs_energy_graph_clusters_sum)
linearity_vs_energy_fit_clusters_dRcut = plot_resolution_vs_energy_graph('linearityclustersdRcut', "", linearity_vs_energy_graph_clusters_sum)


output_rootfile_graph.Close()
