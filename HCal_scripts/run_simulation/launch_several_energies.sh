# eospath=rootfile_storage
# script=runTopoAndSlidingWindowAndCaloSim.py

eospath=FCC_fellow/FCC_rootfile_storage/HCal_v28Aug24_FSR
script=run_thetamodulemerged_tileStandaloneBarrel.py
# before submitting, consider turning off saving hits/cells/hcal in run_thetamodulemerged.py to save disk space

user=`whoami`
initial=$(printf %.1s "$user")
eosdir=/eos/user/$initial/$user/${eospath}/
echo $eosdir
python condor_submit_fccsw.py -outputFolder $eosdir -campaignName $(date +"%y%m%d")_energies_10kevt_SW_noNoise_HCal \
 -gaudiConfig $PWD/$script -jobType caloReco \
 -energies 1000 2000 5000 10000 20000 30000 45000 100000 150000 200000 \
 -nEvt 10000 -originalNjobs 2 -condorQueue longlunch --submit
