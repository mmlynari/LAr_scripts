import os
from Gaudi.Configuration import *

# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if K4GEO, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use=[
#    'Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectMaster_thetamodulemerged.xml',
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml'    
]

# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import CreateFCCeeCaloNeighbours
neighbours = CreateFCCeeCaloNeighbours("neighbours", 
                                       outputFileName = "neighbours_map_barrel_thetamodulemerged.root",
                                       readoutNamesModuleTheta = ["ECalBarrelModuleThetaMerged"],
#                                       readoutNamesModuleTheta = ["ECalBarrelModuleThetaMerged2"],
                                       systemNamesModuleTheta = ["system"],
                                       systemValuesModuleTheta = [4],
                                       activeFieldNamesModuleTheta = ["layer"],
                                       activeVolumesNumbers = [12],
                                       #activeVolumesTheta = [1.2524, 1.2234, 1.1956, 1.1561, 1.1189, 1.0839, 1.0509, 0.9999, 0.9534, 0.91072],
                                       includeDiagonalCells = False,
                                       readoutNamesVolumes = [],
                                       connectBarrels = False, 
                                       OutputLevel = DEBUG)

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [geoservice, neighbours],
                OutputLevel=INFO
)
