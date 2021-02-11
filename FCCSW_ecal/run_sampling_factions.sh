# Careful, you can not provide the rootfile name without '?', it crashes
python FCC_calo_analysis_cpp/plot_samplingFraction.py /eos/user/b/brfranco/rootfile_storage/210122_condor_calib_10kEvt/calibration_output_pdgID_22_pMin_?_pMax_?_thetaMin_90_thetaMax_90.root 10 -r 10000 --preview -outputfolder FCC_calo_analysis_cpp/plots_sampling_fraction_$(date +"%y%m%d") --sed

