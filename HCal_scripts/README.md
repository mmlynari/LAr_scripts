# HCal and ECal+HCal scripts 

Various scripts for HCal studies and combined ECal+HCal simulation

## HCal sampling fraction (SF) calculation
For this, one needs to first remove ECal which is in front of HCal from the geometry xml file and it in the simulation.  
Set invSF=1 in the run_thetamerged_tileStandalone.py script and run the simulation with 100 GeV electrons (for EM scale) or pions (for HAD scale). 
Basic scripts to obtain the SF and invSF are in the calibration directory. 

## Benchmark calibration for ECal+HCal simulation of charged hadrons
HCal should be calibrated to HAD scale via invSF, ECal on EM scale. 
The benchmark calibration as it is implemented now has 5 parameters
- p[0] brings ECal to HAD scale, 
- p[1] scales the HCal energy, but as it is already on HAD scale, it is fixed to 1 
- p[2] corrects for the energy lost between ECal and HCal 
- p[3] corrects for energy loses in the upstream material (e.g. ECal cryostat)
- p[4] is for the residual corrections (currently set to 0) 
For FCC-ee, it was observed that in the energy range of interest, the benchmark parameters depend on the energy (unlike in FCC-hh case). 
Therefore, one needs to obtain benchmark parameters for various energies, plot these distributions to find out the dependency of each parameter on the energy. 
Then these formulas can be used to correct cluster energy by applying CorrectCaloClusters. Note, that one also needs to provide the approximate benchmark parameters 
for the initial energy calculation (currently correspond to 50 GeV pion) - this approximate energy is then used to obtain final benchmark parameters. 
 


