from ROOT import TH1F, TCanvas, TLegend, TFile, gStyle, gPad
import ROOT
import itertools
from datetime import date
import os
from math import ceil, sin, cos, atan, exp, log, tan, pi, sqrt, asin, degrees, radians

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetPadTickY(1)

# you know the noise current rms from Martin
# you need to further get, for each layer, the peak current corresponding to an energy deposit of 1 MeV in the cell (cell considered including the energy in absorber and PCB), the cell merging strategy does not matter yet here to first approximation because the 1 MeV current equivalent will be the same in any merging scenarios (the current will be shared in more gaps when merging many cells, but all these current will then bu 'summed' before to reach the readout). Merging many cells will pay back later, when more signal will be collected per read out channel compared to merging less cells, for the same noise values
# Parameters impacting the peak current induce by 1 MeV deposit (susceptible to change): sampling fraction, LAr gap size (from ramo shockley, the closer the plate the more the charge drift induce a high current), drift velocity (the electric field is smaller when large gap size). Not changing: critical energy LAr, 

noise_current_rms = 10 # nA, FIXME in a next iteration, should be defined as a linear function linking cell capacitance and noise ENI

SFfcc = [0.36571381189697705] * 1 + [0.09779064189677973] * 1 + [0.12564152224404024] * 1 + [0.14350599973146283] * 1 + [0.1557126972314961] * 1 + [0.16444759076233928] * 1 + [0.17097165096847836] * 1 + [0.17684775359805122] * 1 + [0.18181154293837265] * 1 + [0.18544247938196395] * 1 + [0.18922747431624687] * 1 + [0.21187001375505543] * 1

readoutLayerRadialLengths = [1.500000] * 1 + [3.500000] * 11 # cm

qe = 1.602*pow(10, -19) # Coulomb
r_recomb = 0.04
v_lar = 4.75 # mm/microsec, at 87 Kelvin and 1 kV/mm (neglect for now the fact that we might not have 1 kV/mm on every layer due to the lar_gap widening, if we have two HV suply for top and bottom, the electric field will vary by 1/3 which lead to up 13% variation in drift velocity - drift velocity not linear with electric field)
w_lar = 23.6 # eV needed to create a ion/electron pair

nLayers = len(SFfcc)
rmin = 2160
Nplanes = 1536
inclination_degree = 50
angle = inclination_degree / 180. * pi # inclination angle in radian
passiveThickness = 2.0 # mm
pcbThickness = 1.2 # mm

activeTotal = 400.0
inclinedTotal = 564.964
dilution_factor = inclinedTotal / activeTotal

#Segmentation
deltaEta = 0.01
maxEta = 0.881 # 45 degrees
nbins = int(ceil(maxEta/deltaEta))

def get_ref_current(SF, d_lar, E_dep = 1 * pow(10, 6)): #E_dep en eV, choose 1 MeV
    #return E_dep * SF * (1 - r_recomb) * qe * v_lar / (w_lar * d_lar) # Coulomb/microsec
    return E_dep * SF * (1 - r_recomb) * qe * v_lar * pow(10, 6) * pow(10, 9) / (w_lar * d_lar) # nA

# compute LAr gap size for each layer
real_radial_separation = [rmin]
real_radial_depth = []
inclinations_wrt_radial_dir_at_middleRadialDepth = []
lar_gap_sizes_perp = []
ref_current_1mev = []
noise_per_layer = []
# first get the radius of the middle of each layer (PCB cell length is constant parallel to the PCB --> it is not constant w.r.t. the radial direction)
# It will be used later to get the LAr gap size in the direction perpendicular to the plates
current_electrode_length = 0
for idx in range(nLayers):
    readoutLayerRadialLengths[idx] *= 10 # to get mm
    parallel_length = readoutLayerRadialLengths[idx] * dilution_factor
    # Tricky point: in the xml geo, you define 'radial'segmentation, but these depths will be the one parallel to the plates after scaling by the dilution factor --> even when setting constant radial depth, the geoemtry builder will make constant parallel length step, not constant radial steps
    current_electrode_length += parallel_length
    real_radial_separation.append(sqrt(rmin * rmin + current_electrode_length * current_electrode_length + 2 * rmin * current_electrode_length * cos(angle)))
    real_radial_depth.append(real_radial_separation[idx+1] - real_radial_separation[idx])
    # treating the fact that radial angle decreases when radial depth increase
    # angle comprise by lines from  1) Interaction point to inner right edge of a cell, 2) Interaction point to outer left edge of the considered cell (useful to get the plate angle with radial direction that changes with increasing R)
    # based on scalene triangle sine law A/sin(a) = B/sin(b) = C/sin(c) (outer left edge aligned on the Y axis)
    inclinations_wrt_radial_dir_at_middleRadialDepth.append(asin(rmin * sin(angle) / (real_radial_separation[idx] + ((real_radial_separation[idx+1] - real_radial_separation[idx]) / 2))))

for idx in range(nLayers):
    # get the cell size perpendicular to the plate direction from the cell size on the circle at given radius and the inclination w.r.t. radial dir, then remove the PCB and lead thickness (no need for any factor here because we are perpendicular to the PCB and lead plates) --> gives the LAr gap sizxe perpendicular
    lar_gap_sizes_perp.append((2 * pi * (real_radial_separation[idx+1] + real_radial_separation[idx]) / 2. / Nplanes * cos (inclinations_wrt_radial_dir_at_middleRadialDepth[idx]) - pcbThickness - passiveThickness) / 2.) # divided by two because two lar gap per cell
    ref_current_1mev.append(get_ref_current(SFfcc[idx], lar_gap_sizes_perp[idx]))
    noise_per_layer.append(noise_current_rms / get_ref_current(SFfcc[idx], lar_gap_sizes_perp[idx]))


SF_rounded_forPrint = []
for SF in SFfcc:
    SF_rounded_forPrint.append(round(SF,2))
print('SF:', SF_rounded_forPrint)

gap_rounded_forPrint = []
for gap in lar_gap_sizes_perp:
    gap_rounded_forPrint.append(round(gap, 2))
print("lar_gap_sizes: ", gap_rounded_forPrint)
print("ref_current_1mev: ", ref_current_1mev)
print("noise_per_layer [MeV]: ", noise_per_layer)

#output_folder = "noise_capa_" + date.today().strftime("%y%m%d") 
output_folder = "noise_constant_vs_capa" 
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)
fSave = TFile(os.path.join(output_folder, "elecNoise_ecalBarrelFCCee.root"),"RECREATE")


GeV = 1000.

gStyle.SetOptStat(0)

# electronic noise histograms
h_elecNoise_fcc = [] # default total noise shield + detector capacitance (without trace capacitance) -> to be used in FCCSW as noise estimation
h_1MevEquivCurrent_fcc = [] # default total noise shield + detector capacitance (without trace capacitance) -> to be used in FCCSW as noise estimation

maximumNoise = 0.

line_color_number = 1
line_style_number = 1
for i in range (0, nLayers):  
    if line_color_number == 10:
        line_color_number = 28
    if line_style_number > 10:
        line_style_number = 1
    #Prepare electronic noise histograms    
    h_elecNoise_fcc.append( TH1F() )
    h_elecNoise_fcc[i].SetLineWidth(3)
    h_elecNoise_fcc[i].SetLineColor(line_color_number)
    h_elecNoise_fcc[i].SetLineStyle(line_style_number)
    h_elecNoise_fcc[i].SetBins(nbins, 0., maxEta)
    h_elecNoise_fcc[i].SetTitle("Default electronic noise; |#eta|; Electronic noise [GeV]")
    h_elecNoise_fcc[i].SetName("h_elecNoise_fcc_"+str(i+1))

    h_1MevEquivCurrent_fcc.append( TH1F() )
    h_1MevEquivCurrent_fcc[i].SetLineWidth(3)
    h_1MevEquivCurrent_fcc[i].SetLineColor(line_color_number)
    h_1MevEquivCurrent_fcc[i].SetLineStyle(line_style_number)
    h_1MevEquivCurrent_fcc[i].SetBins(nbins, 0., maxEta)
    h_1MevEquivCurrent_fcc[i].SetTitle("1 MeV equivalent current; |#eta|; 1 MeV current [nA]")
    h_1MevEquivCurrent_fcc[i].SetName("h_1MevEquivCurrent_fcc_"+str(i+1))

    for ibin in range(0, nbins+1):
        noise = noise_per_layer[i] / GeV
        #find maximum for drawing of histograms
        if noise > maximumNoise:
            maximumNoise = noise
        #fill histogram
        h_elecNoise_fcc[i].SetBinContent(ibin, noise)
        h_1MevEquivCurrent_fcc[i].SetBinContent(ibin, ref_current_1mev[i])
    line_color_number += 1
    line_style_number += 1

#legend = TLegend(0.135,0.573,0.466,0.872)
legend = TLegend(0.135,0.693,0.8,0.892)
legend.SetBorderSize(0)
legend.SetLineColor(0)
legend.SetLineStyle(0)
legend.SetLineWidth(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetHeader("Longitudinal layers")
legend.SetNColumns(4)

i = 0
for h in h_elecNoise_fcc:
    print(h)
    h.SetMinimum(0.)
    h.SetMaximum(maximumNoise*1.5)
    h.GetYaxis().SetTitleOffset(1.4)
    h.Write()
    legend.AddEntry(h, "Layer " + str(i+1),"l")
    i += 1

i = 0
for h in h_1MevEquivCurrent_fcc:
    print(h)
    h.SetMinimum(0)
    h.SetMaximum(13)
    h.GetYaxis().SetTitleOffset(1.1)
    h.Write()
    i += 1

cNoise = TCanvas("cNoise","Electronic noise per cell",800,600)
cNoise.cd()
for i, h in enumerate(h_elecNoise_fcc):
    if i == 0:
        h.Draw("")
    else:
        h.Draw("same")
legend.Draw()
cNoise.Update()
cNoise.Write()
cNoise.Print(os.path.join(output_folder, "cNoise.png"))

cEquivCurrent = TCanvas("cEquivCurrent","Electronic noise per cell",800,600)
cEquivCurrent.cd()
for i, h in enumerate(h_1MevEquivCurrent_fcc):
    if i == 0:
        h.Draw("")
    else:
        h.Draw("same")

legend.Draw()
cEquivCurrent.Update()
cEquivCurrent.Write()
cEquivCurrent.Print(os.path.join(output_folder, "cEquivCurrent.png"))
