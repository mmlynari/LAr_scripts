# HCal and ECal+HCal scripts 

This directory includes various scripts for standalone HCal and combined ECal+HCal simulation and calibration.
This directory reflects the status in the summer 2024 and this version will be used for FCC feasibility study report (FSR). 

## Important notes 
Due to the k4hep software migration, some of these scripts might not work in the newer releases, so here are a few notes 
- the scripts for running the simulation were tested in release /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh -r 2024-08-01 
- BRT training was synchronized with the nightly release on 6 September so it should work
- plotting scripts should also work in the newer releases as these are written in python
- for running the simulation in newer releases, please use [this how-to](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ALLEGRO/ALLEGRO_o1_v03) and 
here is an example code for running the digi+reco step with HCal Barrel and Endcap: [link](https://github.com/HEP-FCC/k4RecCalorimeter/blob/main/RecFCCeeCalorimeter/tests/options/ALLEGRO_o1_v03_digi_reco.py)
 

## HCal sampling fraction (SF) calculation
For this, one needs to first remove ECal and other subdetectors in front of HCal in the geometry xml file and run standalone HCal simulation for a FIXED THETA (e.g. for the barrel 69 degrees). Then you will evaluate the performance using the same theta angle as the one for which the SF was obtained.   
Set invSF=1 in the run_reco_HCal.py script and run the simulation with 100 GeV electrons (for EM scale) or charged pions (for HAD scale). 
Basic scripts to obtain the SF and invSF are in the HCal_SF_calibration directory. The script expects as an input a flat ntuple that can be produced with
[FCCAnalyses caloNtuplizer](https://github.com/HEP-FCC/FCCAnalyses/blob/master/examples/FCCee/fullSim/caloNtupleizer/analysis.py) and an adapted HCal version can be found in 
[FCCAnalyses_updated_scripts directory](FCCAnalyses_updated_scripts/analysis_HCal_new.py).
For the record, the HCal Barrel SFs in Allegro_o1_v03 calculated with 5000 events (28 Aug 2024): 
EM invSF = 30.3953
HAD invSF = 35.2556
 

## Benchmark calibration for ECal+HCal combined simulation of charged hadrons
HCal should be calibrated to HAD scale via invSF, ECal on EM scale. 
The benchmark calibration as it is implemented now has 6 parameters, and two of them are fixed during the minimization (one can try to include them)
- p[0] brings ECal to HAD scale, 
- p[1] scales the HCal energy, but as it is already on HAD scale, it is fixed to 1 
- p[2] corrects for the energy lost between ECal and HCal 
- p[3] corrects ECal energy  
- p[4] corrects for energy loses in the upstream material (e.g. ECal cryostat)
- p[5] is for the residual corrections, it is fixed to 0
 
After running [fcc_ee_caloBenchmarkCalibration.py](benchmark_calibration_scripts/fcc_ee_caloBenchmarkCalibration.py) to obtain the benchmark parameters, the output is stored in a histogram in a root file, where each bin corresponds to one parameter.
For FCC-ee, it was observed that in the energy range of interest, the benchmark parameters depend on the energy (unlike in FCC-hh case). 
Therefore, one needs to obtain benchmark parameters for various energies (in the range cca1GeV-150GeV), plot these distributions to find out the dependency of each parameter on the energy 
and fit with appropriate formula - example script is [draw_benchmark_parameters.py](benchmark_calibration_scripts/draw_benchmark_parameters.py). 
Then these formulas can be used to correct cluster energy by applying CorrectCaloClusters in the simulation code. Note, that one also needs to provide the approximate benchmark parameters 
for the initial energy calculation (currently correspond to 50 GeV pion) - this approximate energy is then used to obtain final benchmark parameters.

## Running the simulation 
Outdated scripts to run the simulation locally or on lxbatch can be found in the simulation_scripts directory, both for standalone HCal and combined ECal+HCal simulation.

New how-to to run the simulation below from Archil (Nov 2024)

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

 

## Energy calibration with BRT
Original scripts were prepared for standalone ECal simulation and its energy calibration, so here is a modified version for HCal and ECal+HCal. 
The scripts communicate with FCCAnalyses caloNtupleizer, so this part needs to be modified, to reflect the correct number of radial layers in each simulation setup and to read the corresponding 
geometry - see [FCCAnalyses_updated_scripts directory](FCCAnalyses_updated_scripts/analysis_HCal_new.py).
For the creation of training and testing samples, one needs to be careful and ensure that these are different events (need to modify the randomSeed for each submitted job) and this is why we have 
two separate lxbatch submission scripts for this. In addition, for the training, you want to have a continuous energy distribution, while for the evaluation, we usually derive the energy resolution
from several energy points. Use BRT_training_launch_several_energies.sh and BRT_validation_launch_several_energies.sh scripts (stored in run_simulation directory) to get training and testing datasets.  
The code for training and plotting is executed from run_all_chain.sh.  

