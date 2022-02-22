import os
from Gaudi.Configuration import *

# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
path_to_detector = os.environ.get("FCCDETECTORS", "")
print(path_to_detector)
detectors_to_use=[
                    'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster.xml',
                  ]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import CreateFCChhCaloNeighbours
neighbours = CreateFCChhCaloNeighbours("neighbours", 
                                       outputFileName = "neighbours_map_barrel.root",
                                       readoutNamesPhiEta = ["ECalBarrelPhiEta"], 
                                       systemNamesPhiEta = ["system"],
                                       systemValuesPhiEta = [4],
                                       activeFieldNamesPhiEta = ["layer"],
                                       activeVolumesNumbers = [12],
                                       #activeVolumesEta = [1.2524, 1.2234, 1.1956, 1.1561, 1.1189, 1.0839, 1.0509, 0.9999, 0.9534, 0.91072],
                                       readoutNamesVolumes = [],
                                       connectBarrels = False, 
                                       OutputLevel = INFO)

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [geoservice, neighbours],
                OutputLevel=INFO
)
