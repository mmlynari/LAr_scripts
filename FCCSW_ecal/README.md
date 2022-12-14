# Explanation of the scripts

- `condor_submit_fccsw.py` is used to submit full sim jobs on the lxplus condor batch system. Check the help message with `python condor_submit_fccsw.py -h`. It will also automatically write the subsequent scripts that you have to run like `hadd.sh`, `rm.sh`, `fcc_analysis_NUM.sh` and `perfPlots.sh` if you use it in the `caloReco` mode. You can also use this to submit sampling fraction jobs (option `-jobType samplingFraction`) which will automatically write a `sf.sh` script to derive the sampling fraction and update it in the config automnatically (see also `FCC_calo_analysis_cpp/plot_samplingFraction.py`) or dead materail correction jobs.
- `deal_with_failed_job.py` is used to diagnose failed job, tells you which one could potentially succeed if you resubmit them and prints a list of rootfile associated to the failed job (remove them before to run the hadd if a small loss of statistics is ok for your purpose). NB: the script parses the .err and check if they contain the expected number of line (which is not 0 even for succesful jobs), you thus have to set this value in line 7 to 13 in `deal_with_failed_job.py`.
- `fcc_ee_samplingFraction_inclinedEcal.py` Gaudi steering file producing samples used to derive the sampling fraction
- `fcc_ee_upstream_inclinedEcal.py` Gaudi steering file producing samples used to derive the dead material corrections
- `fcc_ee_upstream_with_clusters.py` similar to the above but it also keeps cluster collections
- `launch_several_energies.sh` script to submit full sim particle gun jobs for different energies. Useful to e.g. produce energy resolution cruve
- `neighbours.py` Gaudi steering file used to produce the cell neighbor map necessary for Topo clustering
- `noise_launch_several_energies.sh` script to submit full sim particle gun jobs for different energies, including noise (the job splitting is different than in `launch_several_energies.sh` to account for the longer run time)
- `noise_map.py` Gaudi steering file used to produce the noise map file necessary for Topo clustering
- `pi0_condor_submit_fccsw.py` adaptation of the condor_submit script to prioduce pi0/photon samples. #FIXME: `condor_submit_fccsw.py` should be adapted to include this possibility as well instead of having to separated files
- `read_upstream_json.py` script to read the .json produced in dead material correction derivation. It is called in the script automatically produced by `condor_submit_fccsw.py` when running in 'upstream' mode
- `run*.py` (except `runParallel.py`) Gaudi steering files to run full sim in different configurations
- `run_calibration.sh` bash script to launch the derivation of sampling fraction. It first makes sure the `xml` geometries are consistent with `write_calibration_xml.py` (copy the central geometry `xml` and set the absorbers/readout as sensitive) and then call condor_submit in `samplingFraction` mode)
- `run_pythia.sh` bash script to run full sim with Pythia physics events instead of particle gun
- `run_sampling_factions.sh` bash script example to derive sampling fraction once the appropriate samples have been produced. It is called by the scipts automatically produced in `condor_submit_fccsw.py` when running in `samplingFraction` mode
- `run_upstream.sh` bash script to launch full sim and produce samples needed to derive dead material correction
- `write_calibration_xml.py` python script that manipulates the geometry `xml`, writes the one needed for sampling fraction or dead material correction based on the 'regular' one and update the different parameters in the Gaudi steering file with `sed` (`totalNumLayers`, `lastLayerIDs`, etc) 
