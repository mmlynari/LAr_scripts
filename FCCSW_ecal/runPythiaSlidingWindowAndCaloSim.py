import os

from GaudiKernel.SystemOfUnits import MeV, GeV, tesla

################## Pythia
## Pythia generator
from Configurables import PythiaInterface
pythia8gentool = PythiaInterface()
pythia8gentool.Filename = "MCGeneration/ee_Z_ee.cmd"


## Write the HepMC::GenEvent to the data service
from Configurables import GenAlg
pythia8gen = GenAlg()
pythia8gen.SignalProvider = pythia8gentool
pythia8gen.hepmc.Path = "hepmc"

# Input for simulations (momentum is expected in GeV!)
# theta from 80 to 100 degrees corresponds to -0.17 < eta < 0.17 
#thetaMin = 90.
#thetaMax = 90.
magneticField = False

from Gaudi.Configuration import *

from Configurables import FCCDataSvc
podioevent  = FCCDataSvc("EventDataSvc")

################## Particle gun setup
_pi = 3.14159

genParticlesOutputName = "genParticles"
genVerticesOutputName = "genVertices"


### Reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter("Converter")
hepmc_converter.hepmc.Path = "hepmc"
hepmc_converter.genparticles.Path = genParticlesOutputName
hepmc_converter.genvertices.Path = genVerticesOutputName

################## Simulation setup
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
#path_to_detector = os.environ.get("FCC_DETECTORS", "")
path_to_detector = os.environ.get("FCCSWBASEDIR", "")
detectors_to_use=[
                    'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster.xml',
                    #'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster_original.xml',
                  ]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, "..", _det) for _det in detectors_to_use]
geoservice.OutputLevel = WARNING
#geoservice.OutputLevel = DEBUG

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions="SimG4FullSimActions")

# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
if magneticField == 1:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool", FieldComponentZ=-2*tesla, FieldOn=True,IntegratorStepper="ClassicalRK4")
else:
    field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool",FieldOn=False)

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits

# Detector readouts
# ECAL
#ecalBarrelReadoutName = "ECalBarrelTheta"
#ecalBarrelReadoutNamePhiTheta = "ECalBarrelPhiTheta"
ecalBarrelReadoutName = "ECalBarrelEta"
ecalBarrelReadoutNamePhiTheta = "ECalBarrelPhiEta"
# HCAL
hcalReadoutName = "HCalBarrelReadout"
extHcalReadoutName = "HCalExtBarrelReadout"

# Configure saving of calorimeter hits
from Configurables import SimG4SaveCalHits
saveecalbarreltool = SimG4SaveCalHits("saveECalBarrelHits", readoutNames = [ecalBarrelReadoutName])
saveecalbarreltool.positionedCaloHits.Path = "ECalBarrelPositionedHits"
saveecalbarreltool.caloHits.Path = "ECalBarrelHits"

savehcaltool = SimG4SaveCalHits("saveHCalBarrelHits", readoutNames = [hcalReadoutName])
savehcaltool.positionedCaloHits.Path = "HCalPositionedHits"
savehcaltool.caloHits.Path = "HCalBarrelHits"

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.genParticles.Path = genParticlesOutputName

from Configurables import SimG4Alg
geantsim = SimG4Alg("SimG4Alg",
                       outputs= ["SimG4SaveCalHits/saveECalBarrelHits",
                                 "SimG4SaveCalHits/saveHCalBarrelHits",
                       ],
                       eventProvider=particle_converter,
                       OutputLevel=INFO)

############## Digitization (Merging hits into cells, EM scale calibration)
# EM scale calibration (sampling fraction)
from Configurables import CalibrateInLayersTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                   # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                   samplingFraction = [0.306224547517] * 1 + [0.111664096145] * 1 + [0.135828948734] * 1 + [0.151690753987] * 1 + [0.163564788854] * 1 + [0.172627758875] * 1 + [0.180023941499] * 1 + [0.186642705553] * 1 + [0.192229484697] * 1 + [0.197347321559] * 1 + [0.202461342545] * 1 + [0.225390187901] * 1,
                                   readoutName = ecalBarrelReadoutName,
                                   layerFieldName = "layer")

from Configurables import CalibrateCaloHitsTool
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="41.66")

# Create cells in ECal barrel
# 1. step - merge hits into cells with default Theta segmentation
# 2. step - rewrite the cellId using the Phi-Theta segmentation
from Configurables import CreateCaloCells
createEcalBarrelCellsStep1 = CreateCaloCells("CreateECalBarrelCellsStep1",
                               doCellCalibration=True,
                               calibTool = calibEcalBarrel,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="ECalBarrelHits",
                               cells="ECalBarrelCellsStep1")

# Ecal barrel cell positions
from Configurables import CreateVolumeCaloPositions
positionsEcalBarrel = CreateVolumeCaloPositions("positionsBarrelEcal", OutputLevel = INFO)
positionsEcalBarrel.hits.Path = "ECalBarrelCellsStep1"
positionsEcalBarrel.positionedHits.Path = "ECalBarrelPositions"
# Use Phi-Theta segmentation in ECal barrel
from Configurables import RedoSegmentation
resegmentEcalBarrel = RedoSegmentation("ReSegmentationEcal",
                             # old bitfield (readout)
                             oldReadoutName = ecalBarrelReadoutName,
                             # specify which fields are going to be altered (deleted/rewritten)
                             oldSegmentationIds = ["module"],
                             # new bitfield (readout), with new segmentation
                             newReadoutName = ecalBarrelReadoutNamePhiTheta,
                             OutputLevel = INFO,
                             inhits = "ECalBarrelPositions",
                             outhits = "ECalBarrelCellsStep2")

EcalBarrelCellsName = "ECalBarrelCells"
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                               doCellCalibration=False,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="ECalBarrelCellsStep2",
                               cells=EcalBarrelCellsName)

# Ecal barrel cell positions (good for physics, all coordinates set properly)
from Configurables import CellPositionsECalBarrelTool
cellPositionEcalBarrelTool = CellPositionsECalBarrelTool("CellPositionsECalBarrel", readoutName = ecalBarrelReadoutNamePhiTheta, OutputLevel = DEBUG)

from Configurables import CreateCaloCellPositions
createEcalBarrelPositionedCells = CreateCaloCellPositions("ECalBarrelPositionedCells", OutputLevel = DEBUG)
createEcalBarrelPositionedCells.positionsECalBarrelTool = cellPositionEcalBarrelTool
createEcalBarrelPositionedCells.hits.Path = "ECalBarrelCells"
createEcalBarrelPositionedCells.positionedHits.Path = "ECalBarrelPositionedCells"

# Create cells in HCal
# 1. step - merge hits into cells with the default readout
createHcalBarrelCells = CreateCaloCells("CreateHCaloCells",
                               doCellCalibration=True,
                               calibTool=calibHcells,
                               addCellNoise = False, filterCellNoise = False,
                               OutputLevel = INFO,
                               hits="HCalBarrelHits",
                               cells="HCalBarrelCells")

# sliding window clustering
#Empty cells for parts of calorimeter not implemented yet
from Configurables import CreateEmptyCaloCellsCollection
createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
createemptycells.cells.Path = "emptyCaloCells"

from Configurables import CaloTowerTool
towers = CaloTowerTool("towers",
                               deltaEtaTower = 0.01, deltaPhiTower = 2*_pi/768.,
                               ecalBarrelReadoutName = ecalBarrelReadoutNamePhiTheta,
                               ecalEndcapReadoutName = "",
                               ecalFwdReadoutName = "",
                               hcalBarrelReadoutName = "",
                               hcalExtBarrelReadoutName = "",
                               hcalEndcapReadoutName = "",
                               hcalFwdReadoutName = "",
                               OutputLevel = INFO)
towers.ecalBarrelCells.Path = EcalBarrelCellsName
towers.ecalEndcapCells.Path = "emptyCaloCells"
towers.ecalFwdCells.Path = "emptyCaloCells"
towers.hcalBarrelCells.Path = "emptyCaloCells"
towers.hcalExtBarrelCells.Path = "emptyCaloCells"
towers.hcalEndcapCells.Path = "emptyCaloCells"
towers.hcalFwdCells.Path = "emptyCaloCells"

# Cluster variables
windE = 9
windP = 17
posE = 5
posP = 11
dupE = 7
dupP = 13
finE = 9
finP = 17
# approx in GeV: changed from default of 12 in FCC-hh
threshold = 0.5

from Configurables import CreateCaloClustersSlidingWindow
createClusters = CreateCaloClustersSlidingWindow("CreateClusters",
                                                 towerTool = towers,
                                                 nEtaWindow = windE, nPhiWindow = windP,
                                                 nEtaPosition = posE, nPhiPosition = posP,
                                                 nEtaDuplicates = dupE, nPhiDuplicates = dupP,
                                                 nEtaFinal = finE, nPhiFinal = finP,
                                                 energyThreshold = threshold,
                                                 attachCells = True,
                                                 OutputLevel = INFO
                                                 )
createClusters.clusters.Path = "CaloClusters"

################ Output
from Configurables import PodioOutput
out = PodioOutput("out",
                  OutputLevel=INFO)
#Save information about generated particles & calorimeter cells, drop G4 hits and the intermediate steps 
#out.outputCommands = ["drop *", "keep ECalBarrelCells", "keep HCalBarrelCells", "keep GenParticles","keep GenVertices"]
#Save all
out.outputCommands = ["keep *"]

import uuid
#out.filename = "output_fullCalo_SimAndDigi_withCluster_noMagneticField_"+str(momentum)+"GeV_"+uuid.uuid4().hex+".root"
out.filename = "output_fullCalo_SimAndDigi_withCluster_noMagneticField_pythia.root"

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
pythia8gen.AuditExecute = True
hepmc_converter.AuditExecute = True
geantsim.AuditExecute = True
createEcalBarrelCellsStep1.AuditExecute = True
positionsEcalBarrel.AuditExecute = True
resegmentEcalBarrel.AuditExecute = True
createEcalBarrelCells.AuditExecute = True
createHcalBarrelCells.AuditExecute = True
out.AuditExecute = True

from Configurables import EventCounter
event_counter = EventCounter('event_counter')
event_counter.Frequency = 10

from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [
              event_counter,
              pythia8gen,
              hepmc_converter,
              geantsim,
              createEcalBarrelCellsStep1,
              positionsEcalBarrel,
              resegmentEcalBarrel,
              createEcalBarrelCells,
              createEcalBarrelPositionedCells,
              #createHcalBarrelCells,
              createemptycells,
              createClusters,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = 10,
    ExtSvc = [geoservice, podioevent, geantservice, audsvc],
    StopOnSignal = True,
 )
