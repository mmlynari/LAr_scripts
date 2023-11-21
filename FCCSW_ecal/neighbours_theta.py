from Configurables import ApplicationMgr
from Configurables import CreateFCCeeCaloNeighbours
import os
from Gaudi.Configuration import *

# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if K4GEO is empty, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use = [
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ALLEGRO_o1_v02.xml'
]

# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(
    path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
neighbours = CreateFCCeeCaloNeighbours("neighbours",
                                       readoutNames=[
                                       "ECalBarrelModuleThetaMerged"],
                                       systemNames=["system"],
                                       systemValues=[4],
                                       activeFieldNames=["layer"],
                                       activeVolumesNumbers=[12],
                                       includeDiagonalCells=False,
                                       connectBarrels=False,
                                       OutputLevel=DEBUG)

# ApplicationMgr
ApplicationMgr(TopAlg=[],
               EvtSel='NONE',
               EvtMax=1,
               # order is important, as GeoSvc is needed by G4SimSvc
               ExtSvc=[geoservice, neighbours],
               OutputLevel=INFO
               )
