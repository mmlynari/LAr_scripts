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
For this, one needs to first remove ECal and other subdetectors in front of HCal in the geometry xml file and run standalone HCal simulation.  
Set invSF=1 in the [run_thetamerged_tileStandalone.py](run_simulation/run_thetamerged_tileStandalone.py) script and run the simulation with 100 GeV electrons with pdgID=11 (for EM scale) or charged pions with pdgID=211 (for HAD scale). 
Basic scripts to obtain the SF and invSF are in the HCal_SF_calibration directory. The script expects as an input a flat ntuple that can be produced with
[FCCAnalyses caloNtuplizer](https://github.com/HEP-FCC/FCCAnalyses/blob/master/examples/FCCee/fullSim/caloNtupleizer/analysis.py) and an adapted version can be found in 
[FCCAnalyses_updated_scripts directory](FCCAnalyses_updated_scripts/analysis_HCal.py).
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
The scripts to run the simulation locally or on lxbatch can be found in the simulation_scripts directory, both for standalone HCal and combined ECal+HCal simulation. 

## Energy calibration with BRT
Original scripts were prepared for standalone ECal simulation and its energy calibration, so here is a modified version for HCal and ECal+HCal. 
The scripts communicate with FCCAnalyses caloNtupleizer, so this part needs to be modified, to reflect the correct number of radial layers in each simulation setup and to read the corresponding 
geometry - see [FCCAnalyses_updated_scripts directory](FCCAnalyses_updated_scripts/analysis_HCal.py).
For the creation of training and testing samples, one needs to be careful and ensure that these are different events (need to modify the randomSeed for each submitted job) and this is why we have 
two separate lxbatch submission scripts for this. In addition, for the training, you want to have a continuous energy distribution, while for the evaluation, we usually derive the energy resolution
from several energy points. Use BRT_training_launch_several_energies.sh and BRT_validation_launch_several_energies.sh scripts (stored in run_simulation directory) to get training and testing datasets.  
The code for training and plotting is executed from run_all_chain.sh.  

