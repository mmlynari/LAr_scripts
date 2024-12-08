# eospath=rootfile_storage
# script=runTopoAndSlidingWindowAndCaloSim.py

eospath=FCC_fellow/FCC_rootfile_storage/MVA_training_v28Aug24_FSR
script=run_thetamodulemerged_tileStandaloneBarrel.py
# before submitting, consider turning off saving hits/cells/hcal in run_thetamodulemerged.py to save disk space

user=`whoami`
initial=$(printf %.1s "$user")
eosdir=/eos/user/$initial/$user/${eospath}/
echo $eosdir
python BRT_validation_condor_submit_fccsw.py -outputFolder $eosdir -campaignName $(date +"%y%m%d")_energies_10kevt_cells_SW_noNoise_HCal_EMscale_BRT_validation \
 -gaudiConfig $PWD/$script -jobType caloReco \
 -energies 1000 2000 5000 10000 20000 30000 45000 100000 150000 200000 \
 -nEvt 10000 -originalNjobs 2 -energies 1000 2000 3000 4000 5000 10000 15000 20000 30000 45000 75000 100000 120000 150000 -condorQueue longlunch --submit
