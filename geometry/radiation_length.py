from math import *

rhoU = 18.95;   X0U = 0.3166;   RmU = 1.009;    EcU = 6.65
rhoW = 19.30;   X0W = 0.3504;   RmW = 0.9377;   EcW = 7.97
rhoPb = 11.35;  X0Pb = 0.5612;  RmPb = 1.602;   EcPb = 7.43
rhoFe = 7.874;  X0Fe = 1.757;   RmFe = 1.719;   EcFe = 21.68
rhoG10 = 1.800; X0G10 = 17.87;  RmG10 = 6.050;  EcG10 = 62.84
rhoLAr = 1.396; X0LAr = 14.0;   RmLAr = 9.043;  EcLAr = 32.84
rhoLKr = 2.418; X0LKr = 4.703;  RmLKr = 5.857;  EcLKr = 17.03
rhoLXe = 2.953; X0LXe = 2.872;  RmLXe = 5.224;  EcLXe = 11.66
rhoPcb = 1; X0Pcb = 15.92; RmPcb = 1.; EcPcb = 1.;
rhoCryo = 2.699; X0Cryo = 8.896; RmCryo = 1; EcCryo = 1;
# -------- All values set to 1 are just dummy fillers

#Global Cryo and Service space dimensions
cryostatInRadius = 210 #cm
cryostatThicknessIn = 5 #cm
cryostatThicknessOut = 10 #cm
extraMarginIn = 1 #cm
extraMarginOut = 4 #cm

#ratio between pre-sampler radial length and other cells
preSamplerLength = 1.5     #27.5*(1.5/40) # cm
cellLength = preSamplerLength*2.333333    #27.5*(3.5/40) # cm

#dimensions
plateInclination = 40 * pi/180 # first number in degrees
absorberThickness = 0.18 #cm
glueThickness = 0.01 #cm
steelThickness = 0.01 #cm
activeThickness = 0.1239749*2 #cm at Rmin
pcb_thickness = 0.12 # cm

#Set the active and absorber material
absorberX0 = X0Pb
activeX0 = X0LKr



####################################################################
factorThicknessLength = 1/sin(pi/2 - plateInclination)
factorThicknessWidth = 1/cos(plateInclination)


####### section for LAr_gap scaling ###################### Not in use for now
activeX0calc = activeThickness/X0LAr
passiveX0calc = absorberThickness/X0Pb + glueThickness/X0G10 + steelThickness/X0Fe + pcb_thickness/X0Pcb
passiveEqX0 = absorberThickness/X0W + glueThickness/X0G10 + steelThickness/X0Fe + pcb_thickness/X0Pcb

activeLayerThicknessX0calc = passiveEqX0 * (activeX0calc / passiveX0calc)
print("absX0 ratio: %f"%(passiveX0calc/passiveEqX0))
newCellThickness = (activeThickness/(passiveX0calc/passiveEqX0))
print("new active thickness: %f"%(newCellThickness))

# Calculate new LAr_gap size
cellThickness = (absorberThickness + glueThickness + steelThickness + activeThickness + pcb_thickness)*factorThicknessWidth
innerCircumferance = 2.*pi*(cryostatInRadius + cryostatThicknessIn + extraMarginIn)
n_cells_per_module = 64.0
n_contrcution_modules = 16.0
NCells = innerCircumferance/(cellThickness)
print(""); print(""); print("")
print("N cells needed: %f"%NCells)
print("N modules needed: %f"%(NCells/n_cells_per_module))
NModulesRounded = round(NCells/n_cells_per_module)
NCellsRounded = NModulesRounded * n_cells_per_module
print("N modules rounded: %f"%NModulesRounded)
print("N Cells rounded: %f"%NCellsRounded)
new_gap_size_inner = ((innerCircumferance/(NCellsRounded*factorThicknessWidth))- absorberThickness - glueThickness - steelThickness - pcb_thickness)
print("New Gap size = %f*2"%(new_gap_size_inner/2.))
####################





#uncomment when using pre-defined active gap
new_gap_size_inner = activeThickness

totalX0 = 0; NLayersNeeded = 0; totalLength = 0; NLayersPreSamp = 0; preSamplerLengthTotal = 0; new_gap_size = 0; gap_growth = 0;
layerThickness = (absorberThickness + glueThickness + steelThickness + new_gap_size_inner + pcb_thickness)*factorThicknessLength
nonActiveThickness = (absorberThickness + glueThickness + steelThickness + pcb_thickness)*factorThicknessLength

dPhi = 2*pi/NCellsRounded

print("")
print("Radial length layer 0 = %f"%(layerThickness)); print("")

#include the cryo wall and service space in the X0 calculation
totalX0 += cryostatThicknessIn/X0Cryo
totalX0 += cryostatThicknessOut/X0Cryo
totalX0 += extraMarginIn/activeX0
totalX0 += extraMarginOut/activeX0


#Start the calculation loop
preSamplerLength *= 0.95 # so the pre sampler doesn't get an extra layer when totalLenght = ~1.49

while totalX0 < 22. :
  if totalLength > 100:
    break
  if totalLength < preSamplerLength and NLayersNeeded < 2:
    new_gap_size = new_gap_size_inner
    if NLayersNeeded > 0:
      gap_growth = ((NLayersNeeded*layerThickness)/factorThicknessLength)*tan(dPhi)*factorThicknessLength
      layerThickness += gap_growth
      new_gap_size = gap_growth + new_gap_size_inner*factorThicknessLength
    totalX0 += (factorThicknessLength*absorberThickness)/X0G10
    totalX0 += (factorThicknessLength*new_gap_size)/activeX0
    totalX0 += (factorThicknessLength*steelThickness)/X0G10
    totalX0 += (factorThicknessLength*glueThickness)/X0G10
    totalX0 += (factorThicknessLength*pcb_thickness)/X0Pcb
    totalLength += (nonActiveThickness + new_gap_size)
    preSamplerLengthTotal += layerThickness
    NLayersPreSamp += 1
    NLayersNeeded += 1
    print("layer %i LAr_gap = %f*2      totalX0 = %f totalRadialLength = %f pre-sampler "%(NLayersNeeded,(gap_growth+new_gap_size_inner)/2,totalX0, totalLength))
  else:
    if NLayersNeeded > 0:
      gap_growth = (NLayersNeeded*layerThickness/factorThicknessLength)*tan(dPhi)*factorThicknessLength
      layerThickness += gap_growth
      new_gap_size = gap_growth + new_gap_size_inner*factorThicknessLength
    totalX0 += (factorThicknessLength*absorberThickness)/absorberX0
    totalX0 += (factorThicknessLength*new_gap_size)/activeX0
    totalX0 += (factorThicknessLength*steelThickness)/X0Fe
    totalX0 += (factorThicknessLength*glueThickness)/X0G10
    totalX0 += (factorThicknessLength*pcb_thickness)/X0Pcb
    totalLength += (nonActiveThickness + new_gap_size)
    NLayersNeeded += 1
    print("layer %i LAr_gap = %f*2      totalX0 = %f totalRadialLength = %f"%(NLayersNeeded,(gap_growth+new_gap_size_inner)/2,totalX0, totalLength))

print(""); print(""); print("");print(""); print("");
print("----------------")

print("Length pre-sampler: %f"%(preSamplerLengthTotal))
print("Length total Normal layers: %f"%(totalLength-preSamplerLengthTotal))
print("Length single normal layer: %f"%(((totalLength-preSamplerLengthTotal)/11.)))

print('Total X0: %f'%totalX0)
print('Layers Needed: %d'%NLayersNeeded)
print('Radial distance needed: %f [cm]'%(totalLength))


print(""); print(""); print("");print(""); print("");
print("----------------")

print("AbsoberX0 = %f"%((factorThicknessLength*absorberThickness)/absorberX0))
print("Steel X0 = %f"%((factorThicknessLength*steelThickness)/X0Fe))
print("Glue X0 = %f"%((factorThicknessLength*glueThickness)/X0G10))
print("PCB X0 = %f"%((factorThicknessLength*pcb_thickness)/X0Pcb))
print("Active inner X0 = %f"%((factorThicknessLength*new_gap_size_inner)/activeX0))
print("Active outer X0 = %f"%((factorThicknessLength*new_gap_size)/activeX0))
print("Total Layer inner X0 = %f"%(((factorThicknessLength*absorberThickness)/absorberX0) + ((factorThicknessLength*activeThickness)/activeX0) + ((factorThicknessLength*steelThickness)/X0Fe) + ((factorThicknessLength*glueThickness)/X0G10) + ((factorThicknessLength*pcb_thickness)/X0Pcb)))
print("Total Layer inner Thickness = %f"%((absorberThickness + glueThickness + steelThickness + new_gap_size_inner + pcb_thickness)*factorThicknessLength))

print(""); print("");
