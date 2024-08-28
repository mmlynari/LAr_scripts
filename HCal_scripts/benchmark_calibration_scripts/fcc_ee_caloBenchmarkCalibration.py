from Configurables import ApplicationMgr
from Configurables import EventCounter
from Configurables import AuditorSvc, ChronoAuditor
from Configurables import PodioOutput
from Configurables import CaloTowerToolFCCee
from Configurables import CreateEmptyCaloCellsCollection
from Configurables import CreateCaloCellPositionsFCCee
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
from Configurables import RedoSegmentation
from Configurables import CreateCaloCells
from Configurables import CalibrateCaloHitsTool
from Configurables import CalibrateInLayersTool
from Configurables import SimG4Alg
from Configurables import SimG4PrimariesFromEdmTool
from Configurables import SimG4SaveCalHits
from Configurables import SimG4ConstantMagneticFieldTool
from Configurables import SimG4Svc
from Configurables import SimG4FullSimActions
from Configurables import SimG4SaveParticleHistory
from Configurables import GeoSvc
from Configurables import HepMCToEDMConverter
from Configurables import GenAlg
from Configurables import FCCDataSvc
#from Gaudi.Configuration import INFO
# , VERBOSE, DEBUG
from Gaudi.Configuration import *

import os

from GaudiKernel.SystemOfUnits import GeV, tesla, mm
from GaudiKernel.PhysicalConstants import pi, halfpi, twopi
from math import cos, sin, tan

use_pythia = False

## script to obtain parameters using the benchmark method 
## to be used when shooting charged pions into ECal+HCal;

# Input for simulations (momentum is expected in GeV!)
# Parameters for the particle gun simulations, dummy if use_pythia is set to True
# theta from 80 to 100 degrees corresponds to -0.17 < eta < 0.17 
Nevts = 10
momentum = 50 # in GeV
#thetaMin = 90.25 # degrees
#thetaMax = 90.25 # degrees
thetaMin = 69. # degrees corresponds to eta = 0.36
thetaMax = 69. # degrees
#thetaMin = 50 # degrees
#thetaMax = 130 # degrees
pdgCode = 211 # 11 electron, 13 muon, 22 photon, 111 pi0, 211 pi+
magneticField = False


podioevent  = FCCDataSvc("EventDataSvc")

################## Particle gun setup
_pi = 3.14159

genAlg = GenAlg()

from Configurables import  MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("ParticleGun_Electron")
pgun.PdgCodes = [pdgCode]
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin = 0
#pgun.PhiMax = 0
pgun.PhiMax = 2 * _pi
pgun.ThetaMin = thetaMin * _pi / 180.
pgun.ThetaMax = thetaMax * _pi / 180.
genAlg.SignalProvider = pgun

genAlg.hepmc.Path = "hepmc"

hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path="hepmc"
genParticlesOutputName = "genParticles"
hepmc_converter.GenParticles.Path = genParticlesOutputName
hepmc_converter.hepmcStatusList = []

################## Simulation setup
# Detector geometry
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use=[
                    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml',
                  ]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions

#from Configurables import SimG4FullSimActions, SimG4Alg, SimG4PrimariesFromEdmTool, SimG4SaveParticleHistory
actions = SimG4FullSimActions()
# Uncomment if history from Geant4 decays is needed (e.g. to get the photons from pi0) and set actions=actions in SimG4Svc + uncomment saveHistTool in SimG4Alg
#actions.enableHistory=True
#actions.energyCut = 0.2 * GeV
#saveHistTool = SimG4SaveParticleHistory("saveHistory")

#from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions=actions)

# Fixed seed to have reproducible results, change it for each job if you split one production into several jobs
# Mind that if you leave Gaudi handle random seed and some job start within the same second (very likely) you will have duplicates
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 4242

# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field
#from Configurables import SimG4ConstantMagneticFieldTool
if magneticField == 1:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool", FieldComponentZ=-2*tesla, FieldOn=True,IntegratorStepper="ClassicalRK4")
else:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits

# Detector readouts
# ECAL
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"
ecalBarrelReadoutName2 = "ECalBarrelModuleThetaMerged2"
# HCAL
hcalBarrelReadoutName = "HCalBarrelReadout"
hcalBarrelReadoutName2 = "BarHCal_Readout_phitheta"

# Configure saving of calorimeter hits
ecalBarrelHitsName = "ECalBarrelPositionedHits"
#from Configurables import SimG4SaveCalHits
saveECalBarrelTool = SimG4SaveCalHits("saveECalBarrelHits", readoutNames = [ecalBarrelReadoutName])
saveECalBarrelTool.CaloHits.Path = ecalBarrelHitsName

hcalBarrelHitsName = "HCalBarrelPositionedHits"
saveHCalTool = SimG4SaveCalHits("saveHCalBarrelHits", readoutNames = [hcalBarrelReadoutName])
saveHCalTool.CaloHits.Path = hcalBarrelHitsName

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
#from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = genParticlesOutputName

#from Configurables import SimG4Alg
geantsim = SimG4Alg("SimG4Alg",
                       outputs= [saveECalBarrelTool,
                                 saveHCalTool,
                                 #saveHistTool
                       ],
                       eventProvider=particle_converter,
                       OutputLevel=INFO)

############## Digitization (Merging hits into cells, EM scale calibration)
# EM scale calibration (sampling fraction)
#from Configurables import CalibrateInLayersTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=[0.3775596654349802] * 1 + [0.13400227700041234] * 1 + [0.14390509963164044] * 1 + [0.14998482026270935] * 1 + [0.15457673722531148] * 1 + [0.15928098152159675] * 1 + [0.1635367867767212] * 1 + [0.16801070646031507] * 1 + [0.1713409944779989] * 1 + [0.17580195406064622] * 1 + [0.17966699467772812] * 1,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")

#from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="35.26") 
# HCal at EM scale 30.3953
# HCal at HAD scale 35.2556

# Create cells in ECal barrel
# 1. step - merge hits into cells with theta and module segmentation
# (module is a 'physical' cell i.e. lead + LAr + PCB + LAr +lead)
# 2. step - rewrite the cellId using the merged theta-module segmentation
# (merging several modules and severla theta readout cells).
# Add noise at this step if you derived the noise already assuming merged cells

# Step 1: merge hits into cells according to initial segmentation
ecalBarrelCellsName = "ECalBarrelCells"
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalBarrel,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        addPosition=True,
                                        OutputLevel=INFO,
                                        hits=ecalBarrelHitsName,
                                        cells=ecalBarrelCellsName)

# Step 2a: compute new cellID of cells based on new readout
# (merged module-theta segmentation with variable merging vs layer)
resegmentEcalBarrel = RedoSegmentation("ReSegmentationEcal",
                                       # old bitfield (readout)
                                       oldReadoutName=ecalBarrelReadoutName,
                                       # specify which fields are going to be altered (deleted/rewritten)
                                       oldSegmentationIds=["module", "theta"],
                                       # new bitfield (readout), with new segmentation (merged modules and theta cells)
                                       newReadoutName=ecalBarrelReadoutName2,
                                       OutputLevel=INFO,
                                       debugPrint=200,
                                       inhits=ecalBarrelCellsName,
                                       outhits="ECalBarrelCellsMerged")

# Step 2b: merge new cells with same cellID together
# do not apply cell calibration again since cells were already
# calibrated in Step 1
ecalBarrelCellsName2 = "ECalBarrelCells2"
createEcalBarrelCells2 = CreateCaloCells("CreateECalBarrelCells2",
                                         doCellCalibration=False,
                                         addCellNoise=False,
                                         filterCellNoise=False,
                                         OutputLevel=INFO,
                                         hits="ECalBarrelCellsMerged",
                                         cells=ecalBarrelCellsName2)

cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
ecalBarrelPositionedCellsName = "ECalBarrelPositionedCells"
createEcalBarrelPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateECalBarrelPositionedCells",
    OutputLevel=INFO
)
createEcalBarrelPositionedCells.positionsTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCells.hits.Path = ecalBarrelCellsName
createEcalBarrelPositionedCells.positionedHits.Path = ecalBarrelPositionedCellsName

cellPositionEcalBarrelTool2 = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel2",
    readoutName=ecalBarrelReadoutName2,
    OutputLevel=INFO
)
createEcalBarrelPositionedCells2 = CreateCaloCellPositionsFCCee(
    "CreateECalBarrelPositionedCells2",
    OutputLevel=INFO
)
createEcalBarrelPositionedCells2.positionsTool = cellPositionEcalBarrelTool2
createEcalBarrelPositionedCells2.hits.Path = ecalBarrelCellsName2
createEcalBarrelPositionedCells2.positionedHits.Path = "ECalBarrelPositionedCells2"


# Create cells in HCal
# 1 - merge hits into cells with the default readout
hcalBarrelCellsName = "HCalBarrelCells"
createHcalBarrelCells = CreateCaloCells("CreateHCalBarrelCells",
                                        doCellCalibration=True,
                                        calibTool=calibHcells,
                                        addCellNoise=False,
                                        filterCellNoise=False,
                                        addPosition=True,
                                        hits=hcalBarrelHitsName,
                                        cells=hcalBarrelCellsName,
                                        OutputLevel=INFO)

# 2 - attach positions to the cells
from Configurables import CellPositionsHCalBarrelPhiThetaSegTool
cellPositionHcalBarrelTool = CellPositionsHCalBarrelPhiThetaSegTool(
    "CellPositionsHCalBarrel",
    readoutName=hcalBarrelReadoutName,
    OutputLevel=INFO
)
hcalBarrelPositionedCellsName = "HCalBarrelPositionedCells"
createHcalBarrelPositionedCells = CreateCaloCellPositionsFCCee(
    "CreateHcalBarrelPositionedCells",
    OutputLevel=INFO
)
createHcalBarrelPositionedCells.positionsTool = cellPositionHcalBarrelTool
createHcalBarrelPositionedCells.hits.Path = hcalBarrelCellsName
createHcalBarrelPositionedCells.positionedHits.Path = hcalBarrelPositionedCellsName

# 3 - compute new cellID of cells based on new readout - removing row information
# We use a RedoSegmentation. Using a RewriteBitField with removeIds=["row"],
# wont work because there are tiles with same layer/theta/phi but different row
# as a consequence there will be multiple cells with same cellID in the output collection
# and this will screw up the SW clustering
hcalBarrelCellsName2 = "HCalBarrelCells2"

# first we create new hits with the readout without the row information
# and then merge them into new cells
rewriteHCalBarrel = RedoSegmentation("ReSegmentationHcal",
                                         # old bitfield (readout)
                                         oldReadoutName=hcalBarrelReadoutName,
                                         # specify which fields are going to be altered (deleted/rewritten)
                                         oldSegmentationIds=["row", "theta", "phi"],
                                         # new bitfield (readout), with new segmentation (merged modules and theta cells)
                                         newReadoutName=hcalBarrelReadoutName2,
                                         OutputLevel=INFO,
                                         debugPrint=200,
                                         inhits=hcalBarrelPositionedCellsName,
                                         outhits="HCalBarrelCellsWithoutRow")

createHcalBarrelCells2 = CreateCaloCells("CreateHCalBarrelCells2",
                                             doCellCalibration=False,
                                             addCellNoise=False,
                                             filterCellNoise=False,
                                             OutputLevel=INFO,
                                             hits=rewriteHCalBarrel.outhits.Path,
                                             cells=hcalBarrelCellsName2)

# 4 - attach positions to the new cells
from Configurables import CellPositionsHCalBarrelPhiThetaSegTool
hcalBarrelPositionedCellsName2 = "HCalBarrelPositionedCells2"
cellPositionHcalBarrelTool2 = CellPositionsHCalBarrelPhiThetaSegTool(
    "CellPositionsHCalBarrel2",
    readoutName=hcalBarrelReadoutName2,
    OutputLevel=INFO
)
createHcalBarrelPositionedCells2 = CreateCaloCellPositionsFCCee(
    "CreateHCalBarrelPositionedCells2",
    OutputLevel=INFO
)
createHcalBarrelPositionedCells2.positionsTool = cellPositionHcalBarrelTool2
createHcalBarrelPositionedCells2.hits.Path = hcalBarrelCellsName2
createHcalBarrelPositionedCells2.positionedHits.Path = hcalBarrelPositionedCellsName2

from Configurables import CalibrateBenchmarkMethod
benchmark_calib = CalibrateBenchmarkMethod("CalibrateBenchmarkMethod",
                                      readoutNames=[ecalBarrelReadoutName2, hcalBarrelReadoutName2],
                                      energy=momentum,
                                      ECalSystemID=4,
                                      HCalSystemID=8,
                                      numLayersECal=12,
                                      firstLayerHCal=0,
                                      fixedParameters=[1,5],
                                      OutputLevel=INFO)
benchmark_calib.ecalBarrelCells.Path = ecalBarrelCellsName2
benchmark_calib.hcalBarrelCells.Path = hcalBarrelCellsName2

THistSvc().Output = ["rec DATAFILE='benchmark_calibration_output_pMin_"+str(momentum*1000)+"_thetaMin_"+str(thetaMin)+".root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().PrintAll=True
THistSvc().AutoSave=True
THistSvc().AutoFlush=False
THistSvc().OutputLevel=INFO

#CPU information
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
benchmark_calib.AuditExecute = True

'''
### not needed 
import uuid
### PODIO algorithm
out = PodioOutput("out", OutputLevel=WARNING)
out.outputCommands = ["keep *", "drop ECalBarrelHits", "drop ECalBarrelCellsStep*", "drop ECalBarrelPositionedHits", "drop HCalBarrelHits", "drop HCalBarrelPositionedHits"]
#out.filename = "fccee_caloBenchmarkCalib_%ideg_%igev_%s.root" % (thetaMin, momentum, uuid.uuid4().hex[0:16])
out.filename = "fccsw_output_pdgID_211_pMin_%i_thetaMin_%i.root" % (momentum*1000, thetaMin)
''' 

#from Configurables import EventCounter
event_counter = EventCounter('event_counter')
event_counter.Frequency = 10

#from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [
              event_counter,
              genAlg,
              hepmc_converter,
              geantsim,
              createEcalBarrelCells,
              createEcalBarrelPositionedCells,
              resegmentEcalBarrel,
              createEcalBarrelCells2,
              createEcalBarrelPositionedCells2,
              createHcalBarrelCells,
              createHcalBarrelPositionedCells,
              rewriteHCalBarrel,
              createHcalBarrelCells2,
              createHcalBarrelPositionedCells2,
              benchmark_calib,
              #out
              ],
                EvtSel = 'NONE',
                EvtMax = Nevts,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [podioevent, geoservice, geantservice, audsvc],
                OutputLevel = WARNING
)
