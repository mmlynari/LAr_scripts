#!/usr/bin/bash

#runname="trapezoid_LKr_1mm_13factor_XGBoost/"
runname="baseline_LAr_noNoise_1_bugfix/"

# Compute sampling fractions
#python runParallel.py --outDir $runname/sampling --nEvt 1000 --energies 20000 --sampling

# Compute upstream corrections
#python runParallel.py --outDir $runname/upstream --nEvt 1000 --energies 1000 5000 10000 15000 20000 30000 50000 75000 100000 --upstream --SF $runname/sampling/SF.json

# Generate clusters for upstream studies
#python runParallel.py --outDir $runname/upstreamProd --nEvt 1000000 --upstreamProd --SF $runname/sampling/SF.json

# Generate clusters for MVA training
#python runParallel.py --outDir $runname/production --nEvt 300000 --production --SF $runname/sampling/SF.json --corrections $runname/upstream/corr_params_1d.json

# Run clustering algs
#python runParallel.py --outDir $runname/clusters2 --nEvt 1000 --energies 500 1000 5000 10000 15000 20000 30000 50000 75000 100000 --clusters --SF $runname/sampling/SF.json --corrections $runname/upstream/corr_params_1d.json

# Compute resolutions and responses
python compute_resolutions.py --inputDir $runname/clusters2 --outFile $runname/results.csv --clusters CaloClusters CorrectedCaloClusters CaloTopoClusters CorrectedCaloTopoClusters

# And make plots
python plot_resolutions.py --outDir $runname --doFits plot $runname/results.csv --all
python plot_resolutions.py --outDir $runname --doFits compare clusters CaloClusters CorrectedCaloClusters CaloTopoClusters CorrectedCaloTopoClusters $runname/results.csv --all
