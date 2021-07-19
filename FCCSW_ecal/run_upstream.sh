python condor_submit_fccsw.py -outputFolder /eos/user/b/brfranco/rootfile_storage/ -campaignName $(date +"%y%m%d")_condor_upstream \
 -gaudiConfig $PWD/fcc_ee_upstream_inclinedEcal.py -jobType upstream \
 -nEvt 10000 -originalNjobs 2 \
 -energies 1000 5000 10000 15000 20000 30000 50000 75000 100000 125000 150000 200000 #--submit
 #-nEvt 10000 -originalNjobs 5 #--submit
 #-energies 1000 5000 10000 15000 20000 25000 30000 50000 75000 100000 125000 150000 200000 -polarAngles 90 -energiesForDifferentPolarAngles 10000 \
# once the job is done, do the hadd and run 
 #-energies 1000 30000 50000 #--submit
