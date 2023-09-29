from Gaudi.Configuration import *
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
import os
path_to_detector = os.environ.get("FCCDETECTORS", "")
print(path_to_detector)
detectors_to_use=[
                    'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster_thetamodulemerged.xml'
                  ]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

ecalBarrelReadoutName = "ECalBarrelModuleThetaMerged"
#hcalBarrelReadoutName = "ECalBarrelPhiEta"
hcalBarrelReadoutName = "HCalBarrelReadout"
BarrelNoisePath = os.environ['FCCBASEDIR']+"/LAr_scripts/data/elecNoise_ecalBarrelFCCee_theta.root"
ecalBarrelNoiseHistName = "h_elecNoise_fcc_"

from Configurables import CellPositionsECalBarrelModuleThetaSegTool
ECalBcells = CellPositionsECalBarrelModuleThetaSegTool("CellPositionsECalBarrel",
                                                       readoutName = ecalBarrelReadoutName)
#                                                       OutputLevel = DEBUG)
#print(ECalBcells)

from Configurables import CreateFCCeeCaloNoiseLevelMap, ReadNoiseFromFileTool, ReadNoiseVsThetaFromFileTool
ECalNoiseTool = ReadNoiseVsThetaFromFileTool("ReadNoiseFromFileToolECal",
                                             useSegmentation = False,
                                             cellPositionsTool = ECalBcells,
                                             readoutName = ecalBarrelReadoutName,
                                             noiseFileName = BarrelNoisePath,
                                             elecNoiseHistoName = ecalBarrelNoiseHistName,
                                             setNoiseOffset = False,
                                             activeFieldName = "layer",
                                             addPileup = False,
                                             numRadialLayers = 12,
                                             scaleFactor = 1/1000., #MeV to GeV
                                             OutputLevel = INFO)

HCalNoiseTool = ReadNoiseFromFileTool("ReadNoiseFromFileToolHCal",
                                      useSegmentation = True,
                                      readoutName = hcalBarrelReadoutName,
                                      noiseFileName = BarrelNoisePath,
                                      elecNoiseHistoName = ecalBarrelNoiseHistName,
                                      setNoiseOffset = False,
                                      activeFieldName = "layer",
                                      addPileup = False,
                                      numRadialLayers = 12,
                                      scaleFactor = 1/1000., #MeV to GeV
                                      OutputLevel = INFO)

noisePerCell = CreateFCCeeCaloNoiseLevelMap("noisePerCell", 
                                            ECalBarrelNoiseTool = ECalNoiseTool, 
                                            ecalBarrelSysId = 4,
                                            HCalBarrelNoiseTool = HCalNoiseTool,
                                            readoutNamesModuleTheta=[ecalBarrelReadoutName],
                                            systemNamesModuleTheta=["system"],
                                            systemValuesModuleTheta=[4],
                                            activeFieldNamesModuleTheta=["layer"],
                                            activeVolumesNumbers = [12],
                                            #activeVolumesEta = [1.2524, 1.2234, 1.1956, 1.1561, 1.1189, 1.0839, 1.0509, 0.9999, 0.9534, 0.91072],
                                            readoutNamesVolumes = [],
                                            outputFileName = "cellNoise_map_electronicsNoiseLevel_thetamodulemerged.root",
                                            OutputLevel = DEBUG)

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [geoservice, noisePerCell],
                OutputLevel=INFO
               )
