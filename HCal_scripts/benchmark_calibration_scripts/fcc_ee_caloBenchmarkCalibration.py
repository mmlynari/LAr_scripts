# k4run fcc_ee_caloBenchmarkCalibration.py
# steering file to obtain ALLEGRO ECal+HCal barrel benchmark parameters 
# goal is to obtain correct energy calibration for simulation of charged pions 

#
# COMMON IMPORTS
#

# Logger
from Gaudi.Configuration import INFO, DEBUG  # , VERBOSE
# units and physical constants
from GaudiKernel.PhysicalConstants import pi
#from Configurables import GenAlg

#
# SETTINGS
#

# - default settings, that can be overridden via CLI
inputfile = "ALLEGRO_sim.root"            # input file produced with ddsim - can be overridden with IOSvc.Input
outputfile = "ALLEGRO_sim_digi_reco.root" # output file produced by this steering file - can be overridden with IOSvc.Output
Nevts = -1                                # -1 means all events in input file (can be overridden with -n or --num-events option of k4run
## energy of generated particles, needs to be changed for every energy point 
momentum=50

# HCal barrel/endcap parameters for digitisation 
## HCal should be calibrated at HAD scale
hcalBarrelLayers = 13
hcalBarrelInvSamplingFraction = 27.42 # 27.42 is at 60 degrees at EM scale 
hcalEndcapInvSamplingFraction = 27.42 # FIXME: to be updated for ddsim

# ECAL barrel parameters for digitisation
ecalBarrelLayers = 11
ecalBarrelSamplingFraction = [0.3800493723322256] * 1 + [0.13494147915064658] * 1 + [0.142866851721152] * 1 + [0.14839315921940666] * 1 + [0.15298362570665006] * 1 + [0.15709704561942747] * 1 + [0.16063717490147533] * 1 + [0.1641723795419055] * 1 + [0.16845490287689746] * 1 + [0.17111520115997653] * 1 + [0.1730605163148862] * 1
ecalBarrelUpstreamParameters = [[0.028158491043365624, -1.564259408365951, -76.52312805346982, 0.7442903558010191, -34.894692961350195, -74.19340877431723]]
ecalBarrelDownstreamParameters = [[0.00010587711361028165, 0.0052371999097777355, 0.69906696456064, -0.9348243433360095, -0.0364714212117143, 8.360401126995626]]
if ecalBarrelSamplingFraction and len(ecalBarrelSamplingFraction)>0:
    assert(ecalBarrelLayers == len(ecalBarrelSamplingFraction))

resegmentECalBarrel = False

#
# ALGORITHMS AND SERVICES SETUP
#
TopAlg = []  # alg sequence
ExtSvc = []  # list of external services

# Event counter
from Configurables import EventCounter
eventCounter = EventCounter("EventCounter",
                            OutputLevel=INFO,
                            Frequency=10)
TopAlg += [eventCounter]
# add a message sink service if you want a summary table at the end (not needed..)
# ExtSvc += ["Gaudi::Monitoring::MessageSvcSink"]

# CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
ExtSvc += [audsvc]

# Detector geometry
# prefix all xmls with path_to_detector
# if K4GEO is empty, this should use relative path to working directory
from Configurables import GeoSvc
import os
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = [
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml'
]
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
geoservice.OutputLevel = INFO
ExtSvc += [geoservice]

# Input/Output handling
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = inputfile
io_svc.Output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]

# Digitisation (merging hits into cells, EM scale calibration via sampling fractions)

# - ECAL readouts
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"      # barrel, original segmentation (baseline)
ecalBarrelReadoutName2 = "ECalBarrelModuleThetaMerged2"    # barrel, after re-segmentation (for optimisation studies)
ecalEndcapReadoutName = "ECalEndcapTurbine"                # endcap, turbine-like (baseline)
# - HCAL readouts
hcalBarrelReadoutName = "HCalBarrelReadout"            # barrel, original segmentation (theta-phi)
hcalEndcapReadoutName = "HCalEndcapReadout"            # endcap, original segmentation


# - EM scale calibration (sampling fraction)
from Configurables import CalibrateInLayersTool
#   * ECAL barrel
calibEcalBarrel = CalibrateInLayersTool("CalibrateECalBarrel",
                                        samplingFraction=ecalBarrelSamplingFraction,
                                        readoutName=ecalBarrelReadoutName,
                                        layerFieldName="layer")
#   * ECAL endcap
calibEcalEndcap = CalibrateInLayersTool("CalibrateECalEndcap",
                                        samplingFraction=[0.16419] * 1 + [0.192898] * 1 + [0.18783] * 1 + [0.193203] * 1 + [0.193928] * 1 + [0.192286] * 1 + [0.199959] * 1 + [0.200153] * 1 + [0.212635] * 1 + [0.180345] * 1 + [0.18488] * 1 + [0.194762] * 1 + [0.197775] * 1 + [0.200504] * 1 + [0.205555] * 1 + [0.203601] * 1 + [0.210877] * 1 + [0.208376] * 1 + [0.216345] * 1 + [0.201452] * 1 + [0.202134] * 1 + [0.207566] * 1 + [0.208152] * 1 + [0.209889] * 1 + [0.211743] * 1 + [0.213188] * 1 + [0.215864] * 1 + [0.22972] * 1 + [0.192515] * 1 + [0.0103233] * 1,
                                        readoutName=ecalEndcapReadoutName,
                                        layerFieldName="layer")

from Configurables import CalibrateCaloHitsTool
# HCAL barrel
calibHCalBarrel = CalibrateCaloHitsTool(
    "CalibrateHCalBarrel", invSamplingFraction=hcalBarrelInvSamplingFraction)  
calibHCalEndcap = CalibrateCaloHitsTool(
    "CalibrateHCalEndcap", invSamplingFraction=hcalEndcapInvSamplingFraction) 

# - cell positioning tools
from Configurables import CellPositionsECalBarrelModuleThetaSegTool
cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
# the noise tool needs the positioning tool, but if I reuse the previous one the code crashes..
cellPositionEcalBarrelToolForNoise = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrelForNoise",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)
if resegmentECalBarrel:
    cellPositionEcalBarrelTool2 = CellPositionsECalBarrelModuleThetaSegTool(
        "CellPositionsECalBarrel2",
        readoutName=ecalBarrelReadoutName2,
        OutputLevel=INFO
    )

from Configurables import CellPositionsECalEndcapTurbineSegTool
cellPositionEcalEndcapTool = CellPositionsECalEndcapTurbineSegTool(
    "CellPositionsECalEndcap",
    readoutName=ecalEndcapReadoutName,
    OutputLevel=INFO
)

from Configurables import CellPositionsHCalPhiThetaSegTool
cellPositionHCalBarrelTool = CellPositionsHCalPhiThetaSegTool(
        "CellPositionsHCalBarrel",
        readoutName=hcalBarrelReadoutName,
        OutputLevel=INFO
    )
cellPositionHCalEndcapTool = CellPositionsHCalPhiThetaSegTool(
        "CellPositionsHCalEndcap",
        readoutName=hcalEndcapReadoutName,
        OutputLevel=INFO
    )


# Create cells in ECal barrel (calibrated and positioned - optionally with xtalk and noise added)
# from uncalibrated cells (+cellID info) from ddsim
ecalBarrelPositionedCellsName = ecalBarrelReadoutName + "Positioned"
from Configurables import CreatePositionedCaloCells
createEcalBarrelCells = CreatePositionedCaloCells("CreatePositionedECalBarrelCells",
                                                  doCellCalibration=True,
                                                  calibTool=calibEcalBarrel,
                                                  positionsTool=cellPositionEcalBarrelTool,
                                                  addCrosstalk=False,
                                                  addCellNoise=False,
                                                  filterCellNoise=False,
                                                  noiseTool=None,
                                                  geometryTool=None,
                                                  OutputLevel=INFO,
                                                  hits=ecalBarrelReadoutName,
                                                  cells=ecalBarrelPositionedCellsName
                                                  )
TopAlg += [createEcalBarrelCells]

# -  now, if we want to also save cells with coarser granularity:
if resegmentECalBarrel:
    # rewrite the cellId using the merged theta-module segmentation
    # (merging several modules and severla theta readout cells).
    # Add noise at this step if you derived the noise already assuming merged cells
    # Step a: compute new cellID of cells based on new readout
    # (merged module-theta segmentation with variable merging vs layer)
    from Configurables import RedoSegmentation
    resegmentEcalBarrelTool = RedoSegmentation("ReSegmentationEcal",
                                               # old bitfield (readout)
                                               oldReadoutName=ecalBarrelReadoutName,
                                               # specify which fields are going to be altered (deleted/rewritten)
                                               oldSegmentationIds=["module", "theta"],
                                               # new bitfield (readout), with new segmentation (merged modules and theta cells)
                                               newReadoutName=ecalBarrelReadoutName2,
                                               OutputLevel=INFO,
                                               debugPrint=200,
                                               inhits=ecalBarrelPositionedCellsName,
                                               outhits="ECalBarrelCellsMerged")

    # Step b: merge new cells with same cellID together
    # do not apply cell calibration again since cells were already
    # calibrated in Step 1
    # noise and xtalk off assuming they were applied earlier
    ecalBarrelPositionedCellsName2 = ecalBarrelReadoutName2 + "Positioned"
    createEcalBarrelCells2 = CreatePositionedCaloCells("CreatePositionedECalBarrelCells2",
                                                       doCellCalibration=False,
                                                       positionsTool=cellPositionEcalBarrelTool2,
                                                       calibTool=None,
                                                       crosstalkTool=None,
                                                       addCrosstalk=False,
                                                       addCellNoise=False,
                                                       filterCellNoise=False,
                                                       OutputLevel=INFO,
                                                       hits="ECalBarrelCellsMerged",
                                                       cells=ecalBarrelPositionedCellsName2)
    TopAlg += [
        resegmentEcalBarrelTool,
        createEcalBarrelCells2,
    ]

# Create cells in ECal endcap (needed if one wants to apply cell calibration,
# which is not performed by ddsim)
ecalEndcapPositionedCellsName = ecalEndcapReadoutName + "Positioned"
createEcalEndcapCells = CreatePositionedCaloCells("CreatePositionedECalEndcapCells",
                                                  doCellCalibration=True,
                                                  positionsTool=cellPositionEcalEndcapTool,
                                                  calibTool=calibEcalEndcap,
                                                  crosstalkTool=None,
                                                  addCrosstalk=False,
                                                  addCellNoise=False,
                                                  filterCellNoise=False,
                                                  OutputLevel=INFO,
                                                  hits=ecalEndcapReadoutName,
                                                  cells=ecalEndcapPositionedCellsName)
TopAlg += [createEcalEndcapCells]

# Apply calibration and positioning to cells in HCal barrel
hcalBarrelPositionedCellsName = hcalBarrelReadoutName + "Positioned"
createHCalBarrelCells = CreatePositionedCaloCells("CreatePositionedHCalBarrelCells",
                                                   doCellCalibration=True,
                                                   calibTool=calibHCalBarrel,
                                                   positionsTool=cellPositionHCalBarrelTool,
                                                   addCellNoise=False,
                                                   filterCellNoise=False,
                                                   hits=hcalBarrelReadoutName,
                                                   cells=hcalBarrelPositionedCellsName,
                                                   OutputLevel=INFO)
TopAlg += [createHCalBarrelCells]

# Create cells in HCal endcap
hcalEndcapPositionedCellsName = hcalEndcapReadoutName + "Positioned"
createHCalEndcapCells = CreatePositionedCaloCells("CreatePositionedHCalEndcapCells",
                                                   doCellCalibration=True,
                                                   calibTool=calibHCalEndcap,
                                                   addCellNoise=False,
                                                   filterCellNoise=False,
                                                   positionsTool=cellPositionHCalEndcapTool,
                                                   OutputLevel=INFO,
                                                   hits=hcalEndcapReadoutName,
                                                   cells=hcalEndcapPositionedCellsName)
TopAlg += [createHCalEndcapCells]


## Calibrate ECal and HCal barrel 
from Configurables import CalibrateBenchmarkMethod
benchmarkCalibration = CalibrateBenchmarkMethod("CalibrateBenchmarkMethod",
                                      readoutNames=[ecalBarrelReadoutName2, hcalBarrelReadoutName],
                                      energy=momentum,
                                      ECalSystemID=4,
                                      HCalSystemID=8,
                                      numLayersECal=11,
                                      firstLayerHCal=0,
                                      fixedParameters=[1,5],
                                      OutputLevel=INFO)
benchmarkCalibration.ecalBarrelCells.Path = ecalBarrelPositionedCellsName
benchmarkCalibration.hcalBarrelCells.Path = hcalBarrelPositionedCellsName

TopAlg += [benchmarkCalibration]

from Configurables import THistSvc
THistSvc().Output = ["rec DATAFILE='benchmark_calibration_output_pMin_" + str(momentum * 1000) + ".root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().PrintAll=True
THistSvc().AutoSave=True
THistSvc().AutoFlush=False
THistSvc().OutputLevel=INFO

#from Configurables import EventCounter
event_counter = EventCounter('event_counter')
event_counter.Frequency = 10

# configure the application
print(TopAlg)
print(ExtSvc)
from k4FWCore import ApplicationMgr
applicationMgr = ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel='NONE',
    EvtMax=Nevts,
    ExtSvc=ExtSvc,
    StopOnSignal=True,
)

for algo in applicationMgr.TopAlg:
    algo.AuditExecute = True
