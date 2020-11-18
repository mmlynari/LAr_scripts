python write_calibration_xml.py ../../Detector/DetFCCeeECalInclined/compact/FCCee_ECalBarrel.xml # make sure the geometry changes are propagated to the calibration xml
python condor_submit_calibration.py /eos/user/b/brfranco/rootfile_storage/ $(date +"%y%d%m")_condor_calib_5kEvt  #launch batch job for different energy points

