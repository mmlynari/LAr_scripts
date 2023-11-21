# Some notes

Upon detector geometry change, the noise values, noise maps and neighbor maps (used for topo clustering) should be updated. This can be done with the following scripts:
- `create_capacitance_file.py` creates a TFile with detector capacitance values taking into account various combination of capa source (shields, transmission lines, ...)
- `create_noise_file_chargePreAmp.py` uses the output of `create_capacitance_file.py` to create a TFile with noise vs eta for the different layers. The noise is derived based on [this](https://indico.cern.ch/event/1066234/contributions/4708987/attachments/2387716/4080914/20220209_Brieuc_Francois_Noble_Liquid_Calorimetry_forFCCee_FCCworkshop2022.pdf#page=7).
- you then have to go to `../FCCSW_ecal` to run `neighbours.py` and `noise_map.py` (make sure you modify `BarrelNoisePath` in the latter script to point to your new noise values)

If using the theta-module segmentation (with inclined modules and merging of cells along theta and/or module directions) then the following scripts should be used instead:
- `create_capacitance_file_theta.py`
- `create_noise_file_chargePreAmp_theta.py`
- you then have to go to `../FCCSW_ecal` and execute `fccrun neighbours_theta.py` and `fccrun noise_map_theta.py`

NB: At the moment those files are provided from `/eos/user/b/brfranco/rootfile_storage/` and `/eos/user/g/gmarchio/rootfile_storage/` (which have restricted writing permissions) as explained in the main README. This should probably be updated to take them from e.g. a public cern box link but meanwhile you can just point `TopoCaloNeighbours` and `TopoCaloNoisyCells` to the local updated files or contact Brieuc/Giovanni to update the files directly there.
