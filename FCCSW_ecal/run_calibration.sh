python write_calibration_xml.py ../../Detector/DetFCCeeECalInclined/compact/FCCee_ECalBarrel.xml # make sure the geometry changes are propagated to the calibration xml
python condor_submit_fccsw.py -outputFolder /eos/user/b/brfranco/rootfile_storage/ -campaignName $(date +"%y%d%m")_condor_calib_10kEvt \
 -gaudiConfig $PWD/fcc_ee_samplingFraction_inclinedEcal.py -jobType samplingFraction \
 -energies 10000 -polarAngles 90 -energiesForDifferentPolarAngles 10000 \
 -nEvt 10000 -originalNjobs 50 --submit
 #-energies 300 1000 10000 50000 100000 -polarAngles 90 70 50 -energiesForDifferentPolarAngles 10000 50000 \
# once the job is done, do the hadd and run 
 # cd FCC_calo_analysis_cpp
 # python FCC_calo_analysis_cpp/plot_samplingFraction.py /eos/user/b/brfranco/rootfile_storage/200712_condor_calib_10kEvt/calibration_output_pdgID_22_pMin_?_pMax_?_thetaMin_90_thetaMax_90.root 10 -r 10000 --preview -outputfolder plots_sampling_fraction_$(date +"%y%d%m") --sed
# After you have to change the sampling fraction in runCaloSim.py and in  upstream material config
