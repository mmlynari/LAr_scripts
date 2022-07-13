#!/usr/bin/bash

runname="trapezoid_LKr_angleScaling_1mm_13factor_noNoise_1/"

# Compute sampling fractions
python runParallel.py --outDir $runname/sampling --nEvt 1000 --energies 10000 --sampling

# Compute upstream corrections
python runParallel.py --outDir $runname/upstream --nEvt 1000 --energies 1000 5000 10000 15000 20000 30000 50000 75000 100000 --upstream --SF $runname/sampling/SF.json

# Run clustering algs
python runParallel.py --outDir $runname/clusters --nEvt 1000 --energies 500 1000 5000 10000 15000 20000 30000 50000 75000 100000 --clusters --SF $runname/sampling/SF.json --corrections $runname/upstream/corr_params_1d.json

# Compute resolutions and responses
python compute_resolutions.py --inputDir $runname/clusters --outFile $runname/results.csv --clusters CaloClusters CorrectedCaloClusters CaloTopoClusters CorrectedCaloTopoClusters

# And make plots
python plot_resolutions.py --outDir $runname --doFits plot $runname/results.csv --all
python plot_resolutions.py --outDir $runname --doFits compare clusters CaloClusters CorrectedCaloClusters CaloTopoClusters CorrectedCaloTopoClusters $runname/results.csv --all
