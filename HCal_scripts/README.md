# HCal and ECal+HCal scripts 

This directory includes various scripts for standalone HCal and combined ECal+HCal simulation and calibration.
This directory reflects the status in the summer 2024 and this version will be used for FCC feasibility study report (FSR). 

## Important notes 
Due to the k4hep software migration, some of these scripts might not work in the newer releases, so here are a few notes 
- the scripts for running the simulation were tested in release source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-12-11 
- BRT training was synchronized with the nightly release on 6 September 2024 and likely it will not work out of the box in the latest release
- plotting scripts might require renaming some variable branches and containers following the latest updates 
- for running the simulation in newer releases, please use [this how-to](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ALLEGRO/ALLEGRO_o1_v03) and 
here is an example code for running the digi+reco step with HCal Barrel and Endcap: [link](https://github.com/HEP-FCC/k4RecCalorimeter/blob/main/RecFCCeeCalorimeter/tests/options/ALLEGRO_o1_v03_digi_reco.py)
 
## General workflow 
- calibrate standalone HCal to EMscale (by changing the sampling fraction) and check performance with single electrons and charged pions (xml file needs to be changed to remove all other subdetectors)
- look at the combined performance of ECal+HCal, one can use BDT calibration or benchmark method, see description below (for both of them updated parameters need to be derived)  
- one can look at the performance of calorimeter cells, sliding window clusters and topological clusters (neighbours and noise maps are required)

## HCal sampling fraction (SF) calculation
For this, one needs to first remove ECal and other subdetectors in front of HCal in the geometry [xml file](https://github.com/mmlynari/LAr_scripts/blob/main/HCal_scripts/k4geo_Tilestandalone_xml/ALLEGRO_o1_v03_tileStandalone.xml) and run standalone HCal simulation for a FIXED THETA (e.g. for the barrel 68-69 degrees). 
Then you will evaluate the performance using the same theta angle as the one for which the SF was obtained. If you intend to evaluate the performance for a different theta angle, you need to recalculate the SF.   
First, set invSamplingFraction=1 in calibHCalBarrel in the [run_reco_HCal.py](https://github.com/mmlynari/LAr_scripts/blob/main/HCal_scripts/run_simulation/run_reco_HCal.py) script and run the simulation with 100 GeV electrons (for EM scale) or charged pions (for HAD scale). 
Basic scripts to obtain the SF and invSF are in the HCal_SF_calibration directory. The script expects as an input a flat ntuple that can be produced with
[FCCAnalyses caloNtuplizer](https://github.com/HEP-FCC/FCCAnalyses/blob/master/examples/FCCee/fullSim/caloNtupleizer/analysis.py) and an adapted HCal version can be found in 
[FCCAnalyses_updated_scripts directory](FCCAnalyses_updated_scripts/analysis_HCal_new.py).
For the record, the HCal Barrel SFs in Allegro_o1_v03 calculated with 5000 events at theta=68deg (28 Aug 2024): 
EM invSF = 30.3953
HAD invSF = 35.2556

### Setup the environment
Use the nightly release
``` 
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```
or stable stack (several releases available)
``` 
source /cvmfs/sw.hsf.org/key4hep/setup.sh
```
### run simulation for a fixed theta
```
ddsim --enableGun --gun.distribution uniform --gun.energy "100*GeV" --gun.thetaMin "68*deg" --gun.thetaMax "68*deg" --gun.particle e- --numberOfEvents 5000 --outputFile ALLEGRO_sim_e.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03_tileStandalone.xml
```
### run reconstruction with invSamplingFraction=1
```
k4run run_reco_HCal.py --IOSvc.input ALLEGRO_sim_e.root --IOSvc.output ALLEGRO_reco_pMin_100000_e.root
```
### create a flat ntuple
```
python analysis_HCal_new.py
```
### calculate the sampling fraction 
```
python SF_calibration.py 
```
### rerun the reconstruction with invSamplingFraction set to the correct value
```
k4run run_reco_HCal.py --IOSvc.input ALLEGRO_sim_e.root --IOSvc.output ALLEGRO_reco_pMin_100000_e.root
```

## Make performance plots
Simulate 10k events for each energy point, you can launch the simulations jobs on lxbatch [launch_condor_ddsim.sh](https://github.com/mmlynari/LAr_scripts/blob/main/HCal_scripts/run_simulation/launch_condor_ddsim.sh). 
An example of a script to make resolution plots using cells information can be found here: [perfPlots_HCal_cells_only.py](https://github.com/mmlynari/LAr_scripts/blob/main/HCal_scripts/plotting_scripts/perfPlots_HCal_cells_only.py). 
One can also compare cells and clusters performance, the script it [perfPlots_HCal_clusters_cells.py](https://github.com/mmlynari/LAr_scripts/blob/main/HCal_scripts/plotting_scripts/perfPlots_HCal_clusters_cells.py)
Note that the clustering parameters might need some tuning (e.g. size of the SW).
More scripts for performance plots can be found in (LAr_scripts)[https://github.com/BrieucF/LAr_scripts], but these scripts are not supported anymore.  

## Running the simulation 
Scripts to run the simulation on lxbatch can be found in the [simulation_scripts](simulation_scripts) directory, for standalone HCal (change the xml file if you want to run the combined ECal+HCal simulation).

Below is the how-to for running the simulation locally from Archil (Nov 2024). NOTE: As of December 2024, all updates were pushed to main k4geo and k4Rec directories, so these can be cloned instead of Archil's github and nightly release can be used. 

### Setup the environment (I am working with the stable release 2024-10-03):
```
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-10-03
```
### Ceckout packages:
```
git clone https://github.com/Archil-AD/k4RecCalorimeter.git
git clone https://github.com/Archil-AD/k4geo.git
```

### Compile packages:
```
cd k4geo/
mkdir build install
cd build/
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ..
k4_local_repo
export K4GEO=$PWD/
cd ../

cd k4RecCalorimeter/
mkdir build install
cd build/
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ..
k4_local_repo
cd ../
```

### Run the simulation
```
ddsim --enableGun --gun.distribution uniform --gun.energy "50*GeV" --gun.thetaMin "10*deg" --gun.thetaMax "170*deg" --gun.particle pi- --numberOfEvents 5000 --outputFile ALLEGRO_sim_pi.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
```

### Run the reconstruction
```
cp /afs/cern.ch/work/a/adurglis/public/ALLEGRO/run_reco_HCal.py .
k4run run_reco_HCal.py --IOSvc.input ALLEGRO_sim_pi.root --IOSvc.output ALLEGRO_reco_pi.root
```

### Check the content of output file
```
podio-dump ALLEGRO_reco_pi.root
```
 -- collection of HCal Barrel and Endcap cells are HCalBarrelReadoutPositioned and HCalEndcapReadoutPositioned, respectively.


### An example of drawing HCal barrel energy response:
```
root -l    ALLEGRO_reco_pi.root
events->Draw("Sum$(HCalBarrelReadoutPositioned.energy)/50.","")
```

### ALLEGRO pandora development:
An effort started towards implementing the Pandora particle flow in ALLEGRO detector benchmark, including TileCal HCal. More information and the first very promising results can be found in presentations below  
- [Presentation by Archil in the FCC Full sim meeting](https://indico.cern.ch/event/1481286/#55-pandorapfa-on-allegro)
- [ALLEGRO Pandora how-to](https://github.com/Archil-AD/ALLEGRO_PandoraPFA/tree/main)

One of the first tasks is to tune EM shower parameters so they are correctly identified. This can be done via changing the setup in the settings xml file to play with EMshower identification algorithm parameters:
```    
<LCEmShowerId>
<RmsLowECut>80</RmsLowECut>
</LCEmShowerId>
```
and this is the place in the code where you can find all other parameters:
- [https://github.com/PandoraPFA/LCContent/blob/master/src/LCPlugins/LCParticleIdPlugins.cc](https://github.com/PandoraPFA/LCContent/blob/master/src/LCPlugins/LCParticleIdPlugins.cc#L229)

## Benchmark calibration for ECal+HCal combined simulation of charged hadrons
HCal should be calibrated to HAD scale via invSF, ECal on EM scale.
The benchmark calibration as it is implemented now has 6 parameters, and two of them are fixed during the minimization (one can try to include them)
- p[0] brings ECal to HAD scale,
- p[1] scales the HCal energy, but as it is already on HAD scale, it is fixed to 1
- p[2] corrects for the energy lost between ECal and HCal
- p[3] corrects ECal energy  
- p[4] corrects for energy loses in the upstream material (e.g. ECal cryostat)
- p[5] is for the residual corrections, for now it is fixed to 0 as it did not improve the results 

After running [fcc_ee_caloBenchmarkCalibration.py](benchmark_calibration_scripts/fcc_ee_caloBenchmarkCalibration.py) to obtain the benchmark parameters, the output is stored in a histogram in a root file>
For FCC-ee, it was observed that in the energy range of interest, the benchmark parameters depend on the energy (unlike in FCC-hh case).
Therefore, one needs to obtain benchmark parameters for various energies (in the range cca1GeV-150GeV), plot these distributions to find out the dependency of each parameter on the energy
and fit with appropriate formula - example script is [draw_benchmark_parameters.py](benchmark_calibration_scripts/draw_benchmark_parameters.py).
Then these formulas can be used to correct cluster energy by applying CorrectCaloClusters in the simulation code. Note, that one also needs to provide the approximate benchmark parameters
for the initial energy calculation (currently correspond to 50 GeV pion) - this approximate energy is then used to obtain final benchmark parameters.


## Energy calibration with BRT
Original scripts were prepared for standalone ECal simulation and its energy calibration, so here is a modified version for HCal and ECal+HCal. 
The scripts take information from FCCAnalyses caloNtupleizer, so this part needs to be modified, to reflect the correct number of radial layers in each simulation setup and to read the corresponding 
geometry - see [FCCAnalyses_updated_scripts directory](FCCAnalyses_updated_scripts).
For the creation of training and testing samples, one needs to be careful and ensure that these are different events. 
In addition, for the training, you want to have a continuous energy distribution (possible with ddsim using momentumMin and momentumMax instead of the energy?), while for the evaluation, we usually derive the energy resolution
from several energy points. Outdated scripts used to prepare training and validation samples can be found in [outdated_scripts](run_simulation/outdated_scripts). The code for training and plotting is executed from run_all_chain.sh.  

