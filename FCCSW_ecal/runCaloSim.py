import os

from GaudiKernel.SystemOfUnits import MeV, GeV, tesla

# Input for simulations (momentum is expected in GeV!)
momentum = 1
# theta from 80 to 100 degrees corresponds to -0.17 < eta < 0.17 
#thetaMin = 45.
thetaMin = 90.
#thetaMax = 135.
thetaMax = 90.
magneticField = True

from Gaudi.Configuration import *

from Configurables import FCCDataSvc
podioevent  = FCCDataSvc("EventDataSvc")

################## Particle gun setup
_pi = 3.14159

from Configurables import  MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("ParticleGun_Electron")
#pgun.PdgCodes = [11] # electron
pgun.PdgCodes = [13] # muon
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin = 0
pgun.PhiMax = 2 * _pi
pgun.ThetaMin = thetaMin * _pi / 180.
pgun.ThetaMax = thetaMax * _pi / 180.

from Configurables import GenAlg
genalg_pgun = GenAlg()
genalg_pgun.SignalProvider = pgun 
genalg_pgun.hepmc.Path = "hepmc"

from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.genparticles.Path="GenParticles"
hepmc_converter.genvertices.Path="GenVertices"

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
ecalBarrelReadoutName = "ECalBarrelTheta"
ecalBarrelReadoutNamePhiTheta = "ECalBarrelPhiTheta"
#ecalBarrelReadoutName = "ECalBarrelEta"
#ecalBarrelReadoutNamePhiTheta = "ECalBarrelPhiEta"
# HCAL
hcalReadoutName = "HCalBarrelReadout"
extHcalReadoutName = "HCalExtBarrelReadout"

# Configure saving of calorimeter hits
from Configurables import SimG4SaveCalHits
saveecalbarreltool = SimG4SaveCalHits("saveECalBarrelHits", readoutNames = [ecalBarrelReadoutName])
saveecalbarreltool.positionedCaloHits.Path = "ECalBarrelPositionedHits"
saveecalbarreltool.caloHits.Path = "ECalBarrelHits"

savehcaltool = SimG4SaveCalHits("saveHCalBarrelHits",readoutNames = [hcalReadoutName])
savehcaltool.positionedCaloHits.Path = "HCalPositionedHits"
savehcaltool.caloHits.Path = "HCalBarrelHits"

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.genParticles.Path = "GenParticles"

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
                                   samplingFraction = [0.303451138049] * 1 + [0.111872504159] * 1 + [0.135806495306] * 1 + [0.151772636618] * 1 + [0.163397436122] * 1 + [0.172566977313] * 1 + [0.179855253903] * 1 + [0.186838417657] * 1 + [0.192865946689] * 1 + [0.197420241611] * 1 + [0.202066552306] * 1 + [0.22646764465] * 1,
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
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                               doCellCalibration=False,
                               addCellNoise=False, filterCellNoise=False,
                               OutputLevel=INFO,
                               hits="ECalBarrelCellsStep2",
                               cells="ECalBarrelCells")
positionsEcalBarrelStep2 = CreateVolumeCaloPositions("ECalBarrelPositionedCells", OutputLevel = INFO)
positionsEcalBarrelStep2.hits.Path = "ECalBarrelCells"
positionsEcalBarrelStep2.positionedHits.Path = "ECalBarrelPositionedCells"

# Create cells in HCal
# 1. step - merge hits into cells with the default readout
createHcalBarrelCells = CreateCaloCells("CreateHCaloCells",
                               doCellCalibration=True,
                               calibTool=calibHcells,
                               addCellNoise = False, filterCellNoise = False,
                               OutputLevel = INFO,
                               hits="HCalBarrelHits",
                               cells="HCalBarrelCells")


################ Output
from Configurables import PodioOutput
out = PodioOutput("out",
                  OutputLevel=INFO)
#Save information about generated particles & calorimeter cells, drop G4 hits and the intermediate steps 
#out.outputCommands = ["drop *", "keep ECalBarrelCells", "keep HCalBarrelCells", "keep GenParticles","keep GenVertices"]
#Save all
out.outputCommands = ["keep *"]

import uuid
out.filename = "output_fullCalo_SimAndDigi_"+str(momentum)+"GeV_"+uuid.uuid4().hex+".root"

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
genalg_pgun.AuditExecute = True
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
              genalg_pgun,
              hepmc_converter,
              geantsim,
              createEcalBarrelCellsStep1,
              positionsEcalBarrel,
              resegmentEcalBarrel,
              createEcalBarrelCells,
              #positionsEcalBarrelStep2,
              createHcalBarrelCells,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = 100,
    ExtSvc = [geoservice, podioevent, geantservice, audsvc],
    StopOnSignal = True,
 )
