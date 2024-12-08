# eospath=rootfile_storage
# script=runTopoAndSlidingWindowAndCaloSim.py

eospath=FCC_fellow/FCC_rootfile_storage/MVA_training_v28Aug24_FSR
script=run_thetamodulemerged_tileStandaloneBarrel.py
# before submitting, consider turning off saving hits/cells/hcal in run_thetamodulemerged.py to save disk space

user=`whoami`
initial=$(printf %.1s "$user")
eosdir=/eos/user/$initial/$user/${eospath}/
echo $eosdir
python BRT_training_condor_submit_fccsw.py -outputFolder $eosdir -campaignName $(date +"%y%m%d")_energies_1mil_cells_SW_noNoise_HCal_EMscale_BRT_training_trial2 \
 -gaudiConfig $PWD/$script -jobType caloReco \
 -nEvt 1000000 -originalNjobs 1000 -condorQueue longlunch --submit

