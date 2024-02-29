# eospath=rootfile_storage
# script=runTopoAndSlidingWindowAndCaloSim.py

eospath=fcc/lar/root
script=run_thetamodulemerged.py
# before submitting, consider turning off saving hits/cells/hcal in run_thetamodulemerged.py to save disk space

user=`whoami`
initial=$(printf %.1s "$user")
eosdir=/eos/user/$initial/$user/${eospath}/
echo $eosdir
python condor_submit_fccsw.py -outputFolder $eosdir -campaignName $(date +"%y%m%d")_energies_10kevt_topoAndSW_noNoise \
 -gaudiConfig $PWD/$script -jobType caloReco \
 -energies 300 650 1000 2000 5000 10000 20000 30000 45000 100000 150000 200000 -polarAngles 90 -energiesForDifferentPolarAngles 10000 \
 -nEvt 10000 -originalNjobs 2 -condorQueue microcentury --submit
