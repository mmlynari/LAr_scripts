# Full chain

## Setup

```
export FCCDETECTORS=$PWD/FCCDetectors/
cd LAr_scripts/FCCSW_ecal
```

What else ?

## Changing geometry

Modify XML

Then make sure the geometry changes are propagated to the calibration xml
`python write_calibration_xml.py ../../FCCDetectors/Detector/DetFCCeeECalInclined/compact/FCCee_ECalBarrel.xml`

Recompute sampling fraction
`../../MyScripts/run_sampling.sh`
includes `--sed` to propagate to various jOs

Recompute upstream/downstream corrections
`python ../../MyScripts/runParallel.py`

Propagate corrections to jOs
`python read_upstream_json.py input.json`

## Run and analyse electron clusters

