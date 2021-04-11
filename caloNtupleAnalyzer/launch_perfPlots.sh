basePath=/afs/cern.ch/user/b/brfranco/work/public/Fellow/FCCSW/FCCAnalysesRepos/FCCAnalyses/210217_caloReco/
for file in $basePath/output_fullCalo_SimAndDigi_withCluster_MagneticField_True_*GeV_ThetaMinMax_90.0_90.0_pdgId_*11_pythiaFalse/*
do
    python perfPlots.py -input $file
done
