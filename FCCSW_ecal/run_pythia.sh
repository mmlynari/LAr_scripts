python condor_submit_fccsw.py -outputFolder /eos/user/b/brfranco/rootfile_storage/ -campaignName $(date +"%y%d%m")_pythia -gaudiConfig $PWD/runSlidingWindowAndCaloSim.py -jobType caloReco -pythia True \
    -nEvt 100000 -originalNjobs 500 --submit
