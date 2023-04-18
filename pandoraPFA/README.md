# Tools for running PandoraPFA on the noble liquid detector concept

Particle Flow Algorithms (PFA) need the full detector information (or at least tracks and calorimeter hits). 
Since we currently do not have track in the drift chamber, one workaround to play with PFA is to take the CLD detector and replace the calorimeter by the Noble Liquid one.

The `CLD_LAr` folder hosts a detector configuration corresponding to CLD from `k4geo` v00-18 where the SiW ECAL Barrel was replaced by the Noble Liquid Barrel one from `FCCDetectors` v0.1pre09 (it is not yet good for physics, just use it to start playing with PFA, once a full Noble Liquid based detector full sim is available we should move to this)

## Recipe

1) Make sure you have installed a local version of FCCDetectors and that you have set the environment variables properly:

To be done only once:
```
git clone https://github.com/HEP-FCC/FCCDetectors
cd FCCDetectors
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ../../
```

To be done each time (in the folder containing FCCDetectors)
```
export FCCDETECTORS=$PWD/FCCDetectors/
PATH=$PWD/FCCDetectors/install/bin/:$PATH
CMAKE_PREFIX_PATH=$PWD/FCCDetectors/install/:$CMAKE_PREFIX_PATH
LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/FCCDetectors/install/python:$PYTHONPATH
LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib64:$LD_LIBRARY_PATH
```

2) Produce events with Simhits from CLD_LAr:
`./run_ddsim_cld_lar.sh`

3) Run the reconstruction on it to get tracks and LAr ECAL deposits:
`fccrun track_gaudi_produce_CLD_events_from_ddsim.py`

4) TO BE IMPLEMENTED: run MarlinPandoraPFA on the output root file of the previous step

## Ongoing attempt

Currently trying to produce CLD_LAr events fully through Gaudi with

`produce_CLD_LAr_events.py`

but it is not working at the moment due to difference in the way Geant4 steps are saved in ddsim and in the FCC framework.

