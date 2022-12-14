import os, sys
import ROOT
from math import cos, pi, tan

ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load("libDelphes")

try:
    ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
    #ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
    pass

input_file_path = "/eos/user/b/brfranco/rootfile_storage/delphes_ee_to_z_to_all/ee_z.root"
#input_file_path = "/eos/user/b/brfranco/rootfile_storage/delphes_ee_to_z_to_all/ee_ee.root"
#input_file_path = "/eos/user/b/brfranco/rootfile_storage/delphes_ee_to_z_to_all/ee_aa.root"
output_prefix = os.path.basename(input_file_path).split(".")[0]
input_rootfile = ROOT.TFile(input_file_path)
events = input_rootfile.Get("Delphes")
print(events.GetEntries())

flux_vs_cosTheta = ROOT.TH1F("flux_vs_cosTheta", "flux_vs_cosTheta", 50, -1, 1)
flux_vs_theta = ROOT.TH1F("flux_vs_theta", "flux_vs_theta", 50, 0, pi)
flux_vs_eta = ROOT.TH1F("flux_vs_eta", "flux_vs_eta", 50, -5, 5)
flux_vs_z = ROOT.TH1F("flux_vs_z", "flux_vs_z", 50, -300, 300)

detector_inner_radius = 216

n_particles = 0
for event in events:
    for particle in event.Particle:
        #print(particle.Status, " ", particle.PID, " ", particle.E)
        if particle.Status != 1:
            continue
        pdgId = particle.PID
        acceptedPID = [22, 11, -11]
        if pdgId in acceptedPID:
            if particle.E > 0.3:
                p4 = ROOT.TLorentzVector(particle.Px, particle.Py, particle.Pz, particle.E)
                flux_vs_cosTheta.Fill(cos(p4.Theta()))
                flux_vs_theta.Fill(p4.Theta())
                flux_vs_eta.Fill(p4.Eta())
                flux_vs_z.Fill(detector_inner_radius/tan(p4.Theta()))
                n_particles += 1

print("Number of particles accepted: ", n_particles)

cos_theta_canvas = ROOT.TCanvas("cos_theta_canvas", "cos_theta_canvas")
flux_vs_cosTheta.SetTitle("")
flux_vs_cosTheta.GetYaxis().SetTitle("#Particles")
flux_vs_cosTheta.GetXaxis().SetTitle("cos(#theta)")
flux_vs_cosTheta.Draw()
cos_theta_canvas.Print("cos_theta_canvas_" + output_prefix + ".png")

cos_theta_canvas_log = ROOT.TCanvas("cos_theta_canvas_log", "cos_theta_canvas_log")
cos_theta_canvas_log.SetLogy()
flux_vs_cosTheta.SetTitle("")
flux_vs_cosTheta.GetYaxis().SetTitle("#Particles")
flux_vs_cosTheta.GetXaxis().SetTitle("cos(#theta)")
flux_vs_cosTheta.Draw()
cos_theta_canvas_log.Print("cos_theta_canvas_log_" + output_prefix + ".png")

theta_canvas = ROOT.TCanvas("theta_canvas", "theta_canvas")
flux_vs_theta.SetTitle("")
flux_vs_theta.GetYaxis().SetTitle("#Particles")
flux_vs_theta.GetXaxis().SetTitle("#theta (rad)")
flux_vs_theta.Draw()
theta_canvas.Print("theta_canvas_" + output_prefix + ".png")

theta_canvas_log = ROOT.TCanvas("theta_canvas_log", "theta_canvas_log")
theta_canvas_log.SetLogy()
flux_vs_theta.SetTitle("")
flux_vs_theta.GetYaxis().SetTitle("#Particles")
flux_vs_theta.GetXaxis().SetTitle("#theta (rad)")
flux_vs_theta.Draw()
theta_canvas_log.Print("theta_canvas_log_" + output_prefix + ".png")

eta_canvas = ROOT.TCanvas("eta_canvas", "eta_canvas")
flux_vs_eta.SetTitle("")
flux_vs_eta.GetYaxis().SetTitle("#Particles")
flux_vs_eta.GetXaxis().SetTitle("#eta")
flux_vs_eta.Draw()
eta_canvas.Print("eta_canvas_" + output_prefix + ".png")

z_canvas = ROOT.TCanvas("z_canvas", "z_canvas")
flux_vs_z.SetTitle("")
flux_vs_z.GetYaxis().SetTitle("#Particles")
flux_vs_z.GetXaxis().SetTitle("z (cm)")
flux_vs_z.Draw()
z_canvas.Print("z_canvas_" + output_prefix + ".png")
