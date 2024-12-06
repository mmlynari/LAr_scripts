# run_reco_HCal.py
# steering file for the ALLEGRO reconstruction

#
# COMMON IMPORTS
#

# Logger
from Gaudi.Configuration import INFO, DEBUG  # , VERBOSE
# units and physical constants
from GaudiKernel.PhysicalConstants import pi

#
# SETTINGS
#

# - general settings
#
inputfile = "ALLEGRO_sim_pi_endcap.root"            # input file produced with ddsim
outputfile = "ALLEGRO_sim_digi_reco_pi_endcap.root" # output file produced by this steering file
Nevts = -1                                # -1 means all events
doSWClustering = False
doTopoClustering = False

#
# ALGORITHMS AND SERVICES SETUP
#
TopAlg = []  # alg sequence
ExtSvc = []  # list of external services


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
io_svc.input = inputfile
io_svc.output = outputfile
ExtSvc += [EventDataSvc("EventDataSvc")]


# - cell positioning tools

# - ECAL readouts
ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"      # barrel, original segmentation (baseline)

from Configurables import CellPositionsECalBarrelModuleThetaSegTool
cellPositionEcalBarrelTool = CellPositionsECalBarrelModuleThetaSegTool(
    "CellPositionsECalBarrel",
    readoutName=ecalBarrelReadoutName,
    OutputLevel=INFO
)

# - HCAL readouts
hcalBarrelReadoutName = "HCalBarrelReadout"            # barrel, original segmentation (row-phi)
hcalEndcapReadoutName = "HCalEndcapReadout"            # endcap, original segmentation

from Configurables import CalibrateCaloHitsTool
# HCAL barrel
calibHCalBarrel = CalibrateCaloHitsTool(
        "CalibrateHCalBarrel", invSamplingFraction="1.")
# HCAL endcap
calibHCalEndcap = CalibrateCaloHitsTool(
        "CalibrateHCalEndcap", invSamplingFraction="1.")  # FIXME: to be updated for ddsim

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

from Configurables import CreatePositionedCaloCells
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


if doSWClustering or doTopoClustering:
    from Configurables import CreateEmptyCaloCellsCollection
    createemptycells = CreateEmptyCaloCellsCollection("CreateEmptyCaloCells")
    createemptycells.cells.Path = "emptyCaloCells"
    TopAlg += [createemptycells]

# Function that sets up the sequence for producing SW clusters given an input cell collection
def setupSWClusters(inputCells,
                    inputReadouts,
                    outputClusters,
                    clusteringThreshold,
                    applyUpDownstreamCorrections,
                    applyMVAClusterEnergyCalibration,
                    addShapeParameters,
                    runPhotonIDTool):

    global TopAlg

    from Configurables import CaloTowerToolFCCee
    from Configurables import CreateCaloClustersSlidingWindowFCCee

    # Clustering parameters
    # - phi-theta window sizes
    windT = 9
    windP = 17
    posT = 5
    posP = 11
    dupT = 7
    dupP = 13
    finT = 9
    finP = 17
    # - minimal energy to create a cluster in GeV (FCC-ee detectors have to reconstruct low energy particles)
    threshold = clusteringThreshold

    from Configurables import CaloTowerToolFCCee
    from Configurables import CreateCaloClustersSlidingWindowFCCee

    towerTool = CaloTowerToolFCCee(outputClusters + "TowerTool",
                                   deltaThetaTower=4 * 0.009817477 / 4, deltaPhiTower=2 * 2 * pi / 1536.,
                                   ecalBarrelReadoutName=inputReadouts.get("ecalBarrel", ""),
                                   ecalEndcapReadoutName=inputReadouts.get("ecalEndcap", ""),
                                   ecalFwdReadoutName=inputReadouts.get("ecalFwd", ""),
                                   hcalBarrelReadoutName=inputReadouts.get("hcalBarrel", ""),
                                   hcalExtBarrelReadoutName=inputReadouts.get("hcalExtBarrel", ""),
                                   hcalEndcapReadoutName=inputReadouts.get("hcalEndcap", ""),
                                   hcalFwdReadoutName=inputReadouts.get("hcalFwd", ""),
                                   OutputLevel=INFO)
    towerTool.ecalBarrelCells.Path = inputCells.get("ecalBarrel", "emptyCaloCells")
    towerTool.ecalEndcapCells.Path = inputCells.get("ecalEndcap", "emptyCaloCells")
    towerTool.ecalFwdCells.Path = inputCells.get("ecalFwd", "emptyCaloCells")
    towerTool.hcalBarrelCells.Path = inputCells.get("hcalBarrel", "emptyCaloCells")
    towerTool.hcalExtBarrelCells.Path = inputCells.get("hcalExtBarrel", "emptyCaloCells")
    towerTool.hcalEndcapCells.Path = inputCells.get("hcalEndcap", "emptyCaloCells")
    towerTool.hcalFwdCells.Path = inputCells.get("hcalFwd", "emptyCaloCells")

    clusterAlg = CreateCaloClustersSlidingWindowFCCee("Create" + outputClusters,
                                                      towerTool=towerTool,
                                                      nThetaWindow=windT, nPhiWindow=windP,
                                                      nThetaPosition=posT, nPhiPosition=posP,
                                                      nThetaDuplicates=dupT, nPhiDuplicates=dupP,
                                                      nThetaFinal=finT, nPhiFinal=finP,
                                                      energyThreshold=threshold,
                                                      energySharingCorrection=False,
                                                      attachCells=True,
                                                      OutputLevel=INFO
                                                      )
    clusterAlg.clusters.Path = outputClusters
    clusterAlg.clusterCells.Path = outputClusters.replace("Clusters", "Cluster") + "Cells"
    TopAlg += [clusterAlg]
    clusterAlg.AuditExecute = True

if doSWClustering:
    # HCAL clusters
    if runHCal:
        CaloClusterInputs = {
            "hcalBarrel": hcalBarrelPositionedCellsName,
            "hcalEndcap": hcalEndcapPositionedCellsName,
        }
        CaloClusterReadouts = {
            "hcalBarrel": hcalBarrelReadoutName,
            "hcalEndcap": hcalEndcapReadoutName,
        }
        setupSWClusters(CaloClusterInputs,
                        CaloClusterReadouts,
                        "CaloClusters",
                        0.04,
                        False,
                        False,
                        False,
                        False)

# Function that sets up the sequence for producing Topo clusters given an input cell collection
def setupTopoClusters(inputCells,
                      inputReadouts,
                      inputPositioningTools,  # TODO: check if we still need these since the cells are positioned..
                      outputClusters,
                      neighboursMap,
                      noiseMap,
                      applyUpDownstreamCorrections,
                      applyMVAClusterEnergyCalibration,
                      addShapeParameters,
                      runPhotonIDTool):

    global TopAlg

    from Configurables import CaloTopoClusterInputTool
    from Configurables import TopoCaloNeighbours
    from Configurables import TopoCaloNoisyCells
    from Configurables import CaloTopoClusterFCCee


    # Clustering parameters
    seedSigma = 4
    neighbourSigma = 2
    lastNeighbourSigma = 0

    # tool collecting the input cells
    topoClusterInputTool = CaloTopoClusterInputTool(outputClusters + "InputTool",
                                                    ecalBarrelReadoutName=inputReadouts.get("ecalBarrel", ""),
                                                    ecalEndcapReadoutName=inputReadouts.get("ecalEndcap", ""),
                                                    ecalFwdReadoutName=inputReadouts.get("ecalFwd", ""),
                                                    hcalBarrelReadoutName=inputReadouts.get("hcalBarrel", ""),
                                                    hcalExtBarrelReadoutName=inputReadouts.get("hcalExtBarrel", ""),
                                                    hcalEndcapReadoutName=inputReadouts.get("hcalEndcap", ""),
                                                    hcalFwdReadoutName=inputReadouts.get("hcalFwd", ""),
                                                    OutputLevel=INFO)
    topoClusterInputTool.ecalBarrelCells.Path = inputCells.get("ecalBarrel", "emptyCaloCells")
    topoClusterInputTool.ecalEndcapCells.Path = inputCells.get("ecalEndcap", "emptyCaloCells")
    topoClusterInputTool.ecalFwdCells.Path = inputCells.get("ecalFwd", "emptyCaloCells")
    topoClusterInputTool.hcalBarrelCells.Path = inputCells.get("hcalBarrel", "emptyCaloCells")
    topoClusterInputTool.hcalExtBarrelCells.Path = inputCells.get("hcalExtBarrel", "emptyCaloCells")
    topoClusterInputTool.hcalEndcapCells.Path = inputCells.get("hcalEndcap", "emptyCaloCells")
    topoClusterInputTool.hcalFwdCells.Path = inputCells.get("hcalFwd", "emptyCaloCells")

    # tool providing the map of cell neighbours
    neighboursTool = TopoCaloNeighbours(outputClusters + "NeighboursMap",
                                        fileName=neighboursMap,
                                        OutputLevel=INFO)

    # tool providing expected noise levels per cell
    noiseTool = TopoCaloNoisyCells(outputClusters + "NoiseMap",
                                   fileName=noiseMap,
                                   OutputLevel=INFO)

    # algorithm creating the topoclusters
    clusterAlg = CaloTopoClusterFCCee("Create" + outputClusters,
                                      TopoClusterInput=topoClusterInputTool,
                                      # expects neighbours map from cellid->vec < neighbourIds >
                                      neigboursTool=neighboursTool,
                                      # tool to get noise level per cellid
                                      noiseTool=noiseTool,
                                      # cell positions tools for all sub - systems
                                      positionsECalBarrelTool=inputPositioningTools.get('ecalBarrel', None),
                                      # positionsEMECTool=inputPositioningTools.get('ecalEndcap', None),
                                      # positionsEMFwdTool=inputPositioningTools.get('ecalFwd', None),
                                      positionsHCalBarrelTool=inputPositioningTools.get('hcalBarrel', None),
                                      positionsHCalBarrelNoSegTool=None,
                                      positionsHCalExtBarrelTool=inputPositioningTools.get('hcalEndcap', None),
                                      # positionsHECTool=inputPositioningTools.get('hcalEndcap', None),
                                      # positionsHFwdTool=inputPositioningTools.get('hcalFwd', None),
                                      noSegmentationHCal=False,
                                      # algorithm parameters
                                      seedSigma=seedSigma,
                                      neighbourSigma=neighbourSigma,
                                      lastNeighbourSigma=lastNeighbourSigma,
                                      OutputLevel=INFO)
    clusterAlg.clusters.Path = outputClusters
    clusterAlg.clusterCells.Path = outputClusters.replace("Clusters", "Cluster") + "Cells"
    TopAlg += [clusterAlg]



if doTopoClustering:
    # HCAL clusters
    CaloTopoClusterInputs = {
            "hcalBarrel": hcalBarrelPositionedCellsName,
            "hcalEndcap": hcalEndcapPositionedCellsName
        }
    CaloTopoClusterReadouts = {
            "hcalBarrel": hcalBarrelReadoutName,
            "hcalEndcap": hcalEndcapReadoutName
        }
    CaloTopoClusterPositioningTools = {
            "ecalBarrel" : cellPositionEcalBarrelTool,
            "hcalBarrel": cellPositionHCalBarrelTool,
            "hcalEndcap": cellPositionHCalEndcapTool,
        }
    setupTopoClusters(CaloTopoClusterInputs,
                          CaloTopoClusterReadouts,
                          CaloTopoClusterPositioningTools,
                          "CaloTopoClusters",
                          "neighbours_map.root",
                          "cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged_hcalB_thetaphi.root",
                          False,
                          False,
                          False,
                          False)


# drop lumi, vertex, DCH, Muons (unless want to keep for event display)
io_svc.outputCommands.append("drop Lumi*")
# io_svc.outputCommands.append("drop Vertex*")
# io_svc.outputCommands.append("drop DriftChamber_simHits*")
io_svc.outputCommands.append("drop MuonTagger*")

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
