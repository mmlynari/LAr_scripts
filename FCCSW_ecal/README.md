# Full chain

## Setup

```
export FCCDETECTORS=$PWD/FCCDetectors/
cd LAr_scripts/FCCSW_ecal
```

Alternatively: 

```
setupFCCCalo
```

## Changing geometry

Modify XML

Then make sure the geometry changes are propagated to the calibration xml
`python write_calibration_xml.py ../../FCCDetectors/Detector/DetFCCeeECalInclined/compact/FCCee_ECalBarrel.xml`

Recompute sampling fraction
`./run_sampling.sh`
includes `--sed` to propagate to various jOs

=> `sed` to be replaced by command line arguments

Recompute upstream/downstream corrections
`python runParallel.py`

Propagate corrections to jOs
`python read_upstream_json.py input.json`

=> to be replaced by command line arguments

## Run and analyse electron clusters

`python runParallel.py`

runs `runTopoAndSlidingWindowAndCaloSum.py`

## Analyse root files and store into csv files

`python compute_resolutions.py --inputDir Blah --outFile Blah/blih.csv --clusters CaloClusters
CorrectedCaloClusters CaloTopoClusters CorrectedCaloTopoClusters`

## Make plots

`python plot_resolutions.py -h`

Has `plot` and `compare` subcommands
