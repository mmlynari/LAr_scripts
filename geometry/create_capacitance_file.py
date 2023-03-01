from ROOT import TH1F, TF1, TF2, TCanvas, TLegend, TFile, gStyle
import ROOT
from math import ceil, sin, cos, atan, exp, log, tan, pi, sqrt, asin, degrees, radians

ROOT.gROOT.SetBatch(ROOT.kTRUE)

gStyle.SetPadTickY(1)

#Dimensions
### FCCee
#ECalConstruction     INFO ECAL cryostat: front: rmin (cm) = 210 rmax (cm) = 215 dz (cm) = 226
#ECalConstruction     INFO ECAL cryostat: back: rmin (cm) = 260 rmax (cm) = 270 dz (cm) = 226
#ECalConstruction     INFO ECAL cryostat: side: rmin (cm) = 215 rmax (cm) = 260 dz (cm) = 10
#ECalConstruction     INFO ECAL services: front: rmin (cm) = 215 rmax (cm) = 216 dz (cm) = 216
#ECalConstruction     INFO ECAL services: back: rmin (cm) = 256 rmax (cm) = 260 dz (cm) = 216
#ECalConstruction     INFO ECAL bath: material = LAr rmin (cm) =  216 rmax (cm) = 256 thickness in front of ECal (cm) = 1 thickness behind ECal (cm) = 4
#ECalConstruction     INFO ECAL calorimeter volume rmin (cm) =  216 rmax (cm) = 256
#ECalConstruction     INFO passive inner material = Lead and outer material = lArCaloSteel thickness of inner part (cm) =  0.14 thickness of outer part (cm) =  0.04 thickness of total (cm) =  0.2 rotation angle = 0.872665
#ECalConstruction     INFO number of passive plates = 1536 azim. angle difference =  0.00409062
#ECalConstruction     INFO  distance at inner radius (cm) = 0.883573 distance at outer radius (cm) = 1.0472
#ECalConstruction     INFO readout material = PCB thickness of readout planes (cm) =  0.12 number of readout layers = 12
#ECalConstruction     INFO thickness of calorimeter (cm) = 40 length of passive or readout planes (cm) =  56.4964
#ECalConstruction     INFO active material = LAr active layers thickness at inner radius (cm) = 0.247949 thickness at outer radious (cm) = 0.479746 making 93.4857 % increase.
#ECalConstruction     INFO active passive initial overlap (before subtraction) (cm) = 0.1 = 50 %


filename = "capacitances_perSource_ecalBarrelFCCee.root"
#filename = "capacitances_ecalBarrelFCCee_nLayer_%d_fromAnalytical.root"%numLayers

# layer 2 require special care as it is separated in several cells and that the shield run beneath the etch: cell 2 signal pad top capa: 0.68 + 0.20 = 0.88, cell 2 signal pad bot: 0.56 + 0.21 = 0.77, cell 3: 0.34 + 2.4 = 2.74, cell 4: 1 + 0.25 = 1.25, cell 5: 1.85 + 0.28 = 2.13 

activeTotal = 400.0
inclinedTotal = 564.964
tracesPerLayer = [6, 1, 1, 0, 0, 1, 2, 3, 4, 5, 6, 7] # only one trace for strip layer because 4 cells instead of one
ncells_strip_layer = 4.0
# careful, this is not really the radial spacing, it is, after dilution, the spacing in the parallel direction --> radial depth spacing will not be constant
readoutLayerRadialLengths = [1.500000] * 1 + [3.500000] * 11
#Detector
rmin = 2160
Nplanes = 1536
inclination_degree = 50
angle = inclination_degree / 180. * pi #inclination angle inn radian
passiveThickness = 2.0 #mm
#Segmentation
deltaEta = 0.01
maxEta = 1.2 # 3.1 m z extent
numEta = int(ceil(maxEta/deltaEta))
print(numEta)
#PCB dimensions [mm]
hhv = 0.1
hs = 0.17
t = 0.035
w = 0.127
ws = 0.250
#hm = 0.250
hm = 0.2075
pcbThickness = 7 * t + 2 * hhv + 2 * hs + 2 * hm #mm
print("pcbThickness: %f"%pcbThickness)
#constants:
# distance from signal trace to shield (HS) - from impedance vs. trace width vs. distance to ground layer 2D plot (Z = 50 Ohm)
# trace width (W) - min value
# trace thickness (T) - min value
# distance from shield to the edge of PCB
#http://www.analog.com/media/en/training-seminars/design-handbooks/Basic-Linear-Design/Chapter12.pdf, page 40
#signal trace
epsilonR = 4.4 # PCB
#conversion factor: 1 inch = 25.4 mm
inch2mm = 25.4
#capa per length from maxwel
capa_per_mm = 0.123 # pF/mm
capa_per_mm_stripLayer = 0.062 # pF/mm
# multiplicative factor
# factor two because we merge two phi cells together, another factor 2 becasue we have two 1) signal pad / shield capa  2) HV plate / absorber capa per cell
nmult = 2
nmult_trace = 1 # for the trace only the number of phi cell merged playes a role
epsilonRLAr = 1.5 # LAr at 88 K
epsilon0 = 8.854/1000. #pF/mm


# Fill the layer length, trace length, etc
readoutLayerParallelLengths = []
real_radial_separation = [rmin]
real_radial_depth = []
inclinations_wrt_radial_dir_at_middleRadialDepth = []
trace_length = []
numLayers = len(readoutLayerRadialLengths)
dilution_factor = inclinedTotal / activeTotal
trace_length_inner = 0
trace_length_outer = 0
outer = False
current_electrode_length = 0
for idx in range(numLayers):# first pass to get all length parallel to the readout, real radial separation, inclination at the middle of the layer
    readoutLayerRadialLengths[idx] *= 10
    parallel_length = readoutLayerRadialLengths[idx] * dilution_factor
    # Tricky point: in the xml geo, you define 'radial'segmentation, but these depths will be the one parallel to the plates after scaling by the dilution factor --> even when setting constant radial depth, the geoemtry builder will make constant parallel length step, not constant radial steps
    readoutLayerParallelLengths.append(parallel_length)
    if outer: # prepare the starting trace length when starting to extract by the back of the PCB
        trace_length_outer += parallel_length
    if tracesPerLayer[idx] == 0 and tracesPerLayer[idx - 1] == 0:
        outer = True
    # sqrt(r**2+(L1+i*L2)**2+2*r*(L1+i*L2)*cos(alpha)) where L1 = 2.68, L2=12.09, r=192, alpha=50)
    current_electrode_length += parallel_length
    real_radial_separation.append(sqrt(rmin * rmin + current_electrode_length * current_electrode_length + 2 * rmin * current_electrode_length * cos(angle)))
    real_radial_depth.append(real_radial_separation[idx+1] - real_radial_separation[idx])
    # treating the fact that radial angle decreases when radial depth increase
    # angle comprise by lines from  1) Interaction point to inner right edge of a cell, 2) Interaction point to outer left edge of the considered cell (useful to get the plate angle with radial direction that changes with increasing R)
    # based on scalene triangle sine law A/sin(a) = B/sin(b) = C/sin(c) (outer left edge aligned on the Y axis)
    inclinations_wrt_radial_dir_at_middleRadialDepth.append(asin(rmin * sin(angle) / (real_radial_separation[idx] + ((real_radial_separation[idx+1] - real_radial_separation[idx]) / 2))))
    

# second pass to get trace lengths
outer = False
for idx in range(numLayers):
    if tracesPerLayer[idx] == 0 and tracesPerLayer[idx - 1] == 0: # we change direction
        outer = True
    if outer:
        trace_length.append(trace_length_outer)
        if idx == numLayers - 1:
            trace_length_outer == 0
            continue
        trace_length_outer -= readoutLayerParallelLengths[idx+1]
    else:
        trace_length.append(trace_length_inner)
        trace_length_inner += readoutLayerParallelLengths[idx]

print('Readout radial lengths originally asked: ', readoutLayerRadialLengths)
print('Readout parallel lengths: ', readoutLayerParallelLengths)
print("Real radial separation: ", real_radial_separation)
print("Real radial depth: ", real_radial_depth)
print("inclinations_wrt_radial_dir_at_middleRadialDepth: ", [degrees(inclinations) for inclinations in inclinations_wrt_radial_dir_at_middleRadialDepth])
print("Signal trace length per layer: ", trace_length)


gStyle.SetOptStat(0)

cImpedance = TCanvas("cImpedance","",600,800)
cImpedance.Divide(1,2)
cImpedance.cd(1)
fImpedance = TF2("fImpedance","60/sqrt([0])*log(1.9*(2*x+[1])/(0.8*y+[1]))",0.04,0.2,0.04,0.2)
fImpedance.SetTitle("Impedance vs trace width and distance to ground")
fImpedance.SetParameters(epsilonR, t)
fImpedance.Draw("colz")
fImpedance.GetXaxis().SetTitle("Distance to ground [mm]")
fImpedance.GetYaxis().SetTitle("Trace width [mm]")
cImpedance.cd(2)
fImpedance1D = TF1("fImpedance1D","60/sqrt([0])*log(1.9*(2*x+[1])/(0.8*[2]+[1]))",0.04,0.2)
fImpedance1D.SetTitle("Impedance vs distance to ground")
fImpedance1D.SetParameters(epsilonR, t, w)
fImpedance1D.Draw()
fImpedance1D.GetXaxis().SetTitle("Distance to ground [mm]")
fImpedance1D.GetYaxis().SetTitle("Impedance [#Omega]")

#prepare the TH1
hCapTrace = []
hCapShield = []
hCapDetector = []
line_color_number = 1
line_style_number = 1
for i in range (0, len(readoutLayerRadialLengths)):
    if line_color_number == 8:
        line_color_number = 22
    if line_style_number == 8:
        line_style_number = 1
    #traces 
    hCapTrace.append(TH1F())
    hCapTrace[i].SetBins(numEta, 0, maxEta)
    hCapTrace[i].SetLineColor(line_color_number)
    hCapTrace[i].SetLineStyle(line_style_number)
    hCapTrace[i].SetLineWidth(2)
    hCapTrace[i].SetTitle("Stripline capacitance; |#eta|; Capacitance [pF]")
    hCapTrace[i].SetName("hCapacitance_traces"+str(i))
    #shields
    hCapShield.append(TH1F())
    hCapShield[i].SetBins(numEta, 0, maxEta)
    hCapShield[i].SetLineColor(line_color_number)
    hCapShield[i].SetLineStyle(line_style_number)
    hCapShield[i].SetLineWidth(2)
    hCapShield[i].SetTitle("Signal pads - ground shields capacitance; |#eta|; Capacitance [pF]")
    hCapShield[i].SetName("hCapacitance_shields"+str(i))
    #area
    hCapDetector.append(TH1F())
    hCapDetector[i].SetBins(numEta, 0, maxEta)
    hCapDetector[i].SetLineColor(line_color_number)
    hCapDetector[i].SetLineStyle(line_style_number)
    hCapDetector[i].SetLineWidth(2)
    hCapDetector[i].SetTitle("Signal pad - absorber capacitance; |#eta|; Capacitance [pF]")
    hCapDetector[i].SetName("hCapacitance_detector"+str(i))
    if line_color_number > 8:
        line_color_number += 10
    else:
        line_color_number += 1
    line_style_number += 1

cTrace = TCanvas("cTrace","",600,400)
cShield = TCanvas("cShield","",600,400)
cDetector = TCanvas("cDetector","",600,400)

legend = TLegend(0.1,0.693,0.8,0.9)
#legend = TLegend(0.135,0.693,0.8,0.892)
legend.SetHeader("Longitudinal layers")
legend.SetNColumns(4)
capa_shield_max = 0
capa_det_max = 0
for i in range (0, len(readoutLayerParallelLengths)):
    print("--------------")
    for index in range(0, numEta):
        eta = index * deltaEta
        # take into account the inclination in eta
        traceLength = trace_length[i] / (sin(2. * atan(exp(-eta))))
        #print("Layer %d trace length %f"%(i+1, traceLength))
        #Trace capacitance (stripline)
        logStripline = log(3.1 * hs / (0.8 * w + t))
        capacitanceTrace = nmult_trace * 1 / inch2mm * 1.41 * epsilonR / logStripline * traceLength
        hCapTrace[i].SetBinContent(index+1, capacitanceTrace)
    
        #Shield capacitance (microstrip)
        cellLength = readoutLayerParallelLengths[i] / (sin(2. * atan(exp(-eta))))
        logMicrostrip = log(5.98 * hm / (0.8 * ws + t))
        # analytical formula 
        #capacitanceShield = nmult * cellLength * tracesPerLayer[i] * 1 / inch2mm * 0.67 * (epsilonR + 1.41) / logMicrostrip
        # from maxwell 
        capacitanceShield = nmult * cellLength * tracesPerLayer[i] * capa_per_mm
        if i == 1: #strip layer has smaller capacitance due to traces running beneath the anti-etch
            capacitanceShield = nmult * cellLength * tracesPerLayer[i] * capa_per_mm_stripLayer
        if capacitanceShield > capa_shield_max:
            capa_shield_max = capacitanceShield
        hCapShield[i].SetBinContent(index+1, capacitanceShield)

        ##Detector area (C = epsilon*A/d)
        #area = ( radius[i] * ( 1 / (tan(2. * atan(exp(- (index + 1) * deltaEta)))) -  1 / (tan(2. * atan(exp(- index * deltaEta))) ) )
        #         + radius[i + 1] * ( 1 / (tan(2. * atan(exp(- (index + 1) * deltaEta)))) -  1 / (tan(2. * atan(exp(- index * deltaEta))) ) )
        #         ) / 2. * (radius[i+1] - radius[i])
        #distance = (radius[i+1] + radius[i]) / 2. * pi / Nplanes * cos (angle) - pcbThickness / 2. - passiveThickness / 2.

        #Detector area (C = epsilon*A/d)
        area = ( real_radial_separation[i] * ( 1 / (tan(2. * atan(exp(- (index + 1) * deltaEta)))) -  1 / (tan(2. * atan(exp(- index * deltaEta))) ) )
                 + real_radial_separation[i + 1] * ( 1 / (tan(2. * atan(exp(- (index + 1) * deltaEta)))) -  1 / (tan(2. * atan(exp(- index * deltaEta))) ) )
                 ) / 2. * (real_radial_separation[i+1] - real_radial_separation[i])
        # get the cell size perpendicular to the plate direction from the cell size on the circle at given radius and the inclination w.r.t. radial dir, then remove the PCB and lead thickness (no need for any factor here because we are perpendicular to the PCB and lead plates) --> gives the LAr gap sizxe perpendicular
        distance = (2 * pi * (real_radial_separation[i+1] + real_radial_separation[i]) / 2. / Nplanes * cos (inclinations_wrt_radial_dir_at_middleRadialDepth[i]) - pcbThickness - passiveThickness) / 2. # divided by two because two lar gap per cell
        distance += hhv #the capa is between signal plate and absorber --> need to add distance between HV plate and signal pad
        distance += t #the capa is between signal plate and absorber --> need to add distance between HV plate and signal pad
        if eta == 0:
            print("LAr gap size (perpendicular) + hhv + t: %f mm"%distance)
        capacitanceDetector = nmult * epsilon0 * epsilonRLAr * area / distance
        if i == 1: #strip layer has smaller capacitance because it is divided in 4 smaller cells
            capacitanceDetector /= ncells_strip_layer
        hCapDetector[i].SetBinContent(index+1, capacitanceDetector)
        if capacitanceDetector > capa_det_max:
            capa_det_max = capacitanceDetector
        if index==0:
            print("layer %d" %(i+1), "eta==0: capacitanceTrace %.0f pF," %capacitanceTrace, "capacitanceShield %.0f pF" %capacitanceShield, "capacitanceDetector %.0f pF" %capacitanceDetector)
            #, "distance %.1f mm" %distance

    #Draw
    cTrace.cd()
    if i == 0:
        hCapTrace[i].Draw()
    else:
        hCapTrace[i].Draw("same")
    legend.AddEntry(hCapTrace[i],"layer %d"%(i+1),"l")
    cShield.cd()
    if i == 0:
        hCapShield[i].Draw()
    else:
        hCapShield[i].Draw("same")
    cDetector.cd()
    if i == 0:
        hCapDetector[i].Draw()
    else:
        hCapDetector[i].Draw("same")

maximum = capa_shield_max

plots = TFile(filename,"RECREATE")

for i in range (0, len(readoutLayerParallelLengths)):
    hCapTrace[i].SetMinimum(0)
    hCapTrace[i].SetMaximum(maximum*1.8)
    hCapTrace[i].Write()
    hCapShield[i].SetMinimum(0)
    hCapShield[i].SetMaximum(capa_shield_max*1.5)
    hCapShield[i].Write()
    hCapDetector[i].SetMinimum(0)
    hCapDetector[i].SetMaximum(capa_det_max*1.5)
    hCapDetector[i].Write()

cTrace.cd()
legend.Draw()
cTrace.Update()
cTrace.Write()
cTrace.Print("capa_trace.png")
cTrace.Print("capa_trace.pdf")
cShield.cd()
legend.Draw()
cShield.Update()
cShield.Write()
cShield.Print("capa_shield.png")
cShield.Print("capa_shield.pdf")
cDetector.cd()
legend.Draw()
cDetector.Update()
cDetector.Write()
cDetector.Print("capa_detector.png")
cDetector.Print("capa_detector.pdf")

fImpedance.Write()
fImpedance1D.Write()

#closeInput = raw_input("Press ENTER to exit")
