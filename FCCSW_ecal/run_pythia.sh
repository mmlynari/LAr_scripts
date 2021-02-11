python condor_submit_fccsw.py -outputFolder /eos/user/b/brfranco/rootfile_storage/ -campaignName $(date +"%y%m%d")_pythia -gaudiConfig $PWD/runSlidingWindowAndCaloSim.py -jobType caloReco -pythia True -pythiaCfg $PWD/MCGeneration/ee_Z_ee.cmd  \
    -nEvt 10000 -originalNjobs 200 --submit
