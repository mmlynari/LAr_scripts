#!/usr/bin/bash

runname="thetamodulemerged_topoclusters"
xmlbasedir=../../k4geo
xmldir=FCCee/ALLEGRO/compact/ALLEGRO_o1_v02
xmlfileFullDet=ALLEGRO_o1_v02
xmlfileECal=ECalBarrel_thetamodulemerged
today=`date +%y%m%d`


# Remake calibration xml files from the main xml file
#
python write_calibration_xml.py $xmlbasedir/$xmldir/$xmlfileECal.xml


# Compute the X0 plot (material upstream and ECAL separately, and then full detector)
#
cd ../geometry

# tracker only
# prepare steering file
cp -f material_scan.py material_scan_tracker.py
sed -i 's%#detcard%'$xmldir/$xmlfileFullDet'_trackeronly.xml%' material_scan_tracker.py
sed -i 's%#suffix%tracker%' material_scan_tracker.py
sed -i 's%#etamax%2.7%' material_scan_tracker.py
sed -i 's%#etabinning%0.1%' material_scan_tracker.py
# scan
fccrun material_scan_tracker.py
# plot vs costheta
python material_plot_vs_theta.py --f out_material_scan_tracker.root --s _tracker -c 1.0
# plot vs theta
python material_plot_vs_theta.py --f out_material_scan_tracker.root --s _tracker -t 0.0

# ecal only
# prepare steering file
cp -f material_scan.py material_scan_ecal.py
sed -i 's%#detcard%'$xmldir/$xmlfileFullDet'_ecalonly.xml%' material_scan_ecal.py
sed -i 's%#suffix%ecal%' material_scan_ecal.py
sed -i 's%#etamax%2.9%' material_scan_ecal.py
sed -i 's%#etabinning%0.1%' material_scan_ecal.py
# scan
fccrun material_scan_ecal.py
# plot vs costheta
python material_plot_vs_theta.py --f out_material_scan_ecal.root --s _ecal -c 1.0
# plot vs theta
python material_plot_vs_theta.py --f out_material_scan_ecal.root --s _ecal -t 0.0

# full detector
# prepare steering file
cp -f material_scan.py material_scan_all.py
sed -i 's%#detcard%'$xmldir/$xmlfileFullDet'.xml%' material_scan_all.py
sed -i 's%#suffix%all%' material_scan_all.py
sed -i 's%#etamax%2.7%' material_scan_all.py
sed -i 's%#etabinning%0.1%' material_scan_all.py
# scan
fccrun material_scan_all.py
# plot vs costheta
python material_plot_vs_theta.py --f out_material_scan_all.root --s _all -c 1.0
# plot vs theta
python material_plot_vs_theta.py --f out_material_scan_all.root --s _all -t 0.0

cd ../FCCSW_ecal/


# Archive the files
#
mkdir -vp $runname/geometry
cp $xmlbasedir/$xmldir/*.xml $runname/geometry
cp ../geometry/plots/*png $runname/geometry
cp ../geometry/plots/*pdf $runname/geometry


# Compute sampling fractions and update scripts
#

# one energy is enough as they are independent of energy and direction
python runParallel.py --outDir $runname/sampling --nEvt 1000 --energies 20000 --sampling

# otherwise, to plot sampling fractions vs energy or direction, one can do
# python runParallel.py --outDir $runname/sampling --nEvt 1000 --energies 1000 10000 50000 100000 --sampling
# cd FCC_calo_analysis_cpp 
# python plot_samplingFraction.py ../$runname/sampling/calibration_sampling_output_energy_?_theta_90.root 1 10 20 50 100 -r 1000 10000 20000 50000 100000 --totalNumLayers 12 --preview -outputfolder plots_sampling_fraction_$today --plotSFvsEnergy
# cd ..
# python runParallel.py --outDir $runname/sampling --nEvt 1000 --energies 10000 --thetas 90 80 70 60 50 --sampling
# cd FCC_calo_analysis_cpp 
# python plot_samplingFraction.py ../$runname/sampling/calibration_sampling_output_energy_10000_theta_?.root 50 60 70 80 90 -r 50 60 70 80 90 --totalNumLayers 12 --preview -outputfolder plots_sampling_fraction_$today --plotSFvsEnergy --theta
# cd ..


# Compute upstream and downstream corrections and update scripts
#
# The script generates samples of particles of various energies, that are used to 
# calculate the profiles of E(upstream) vs E(layer 0) and of E(downstream) vs E(layer -1)
# which are then fitted with some parametric functions
# The values of the parameters vs particle energy are then fitted to obtain
# a parameterisation of the corrections vs energy
#
python runParallel.py --outDir $runname/upstream --nEvt 1000 --energies 1000 5000 10000 15000 20000 30000 50000 75000 100000 --upstream --SF $runname/sampling/SF.json


# Generate clusters for upstream studies (??)
#
# python runParallel.py --outDir $runname/upstreamProd --nEvt 1000000 --upstreamProd --SF $runname/sampling/SF.json

# Generate clusters for MVA training
# Only 300k events here. Move up to 3M if needed
python runParallel.py --outDir $runname/production --nEvt 3000000 --production --SF $runname/sampling/SF.json --corrections $runname/upstream/corr_params_1d.json

# Then train the MVA on CaloClusters and CaloTopoClusters
# python training.py CaloClusters -i $runname/production/ --json $runname/upstream/corr_params_1d.json -o $runname/training_calo.json
# python training.py CaloTopoClusters -i $runname/production/ --json $runname/upstream/corr_params_1d.json -o $runname/training_topo.json

# Run clustering algs to compute resolutions
# python runParallel.py --outDir $runname/clusters --nEvt 1000 --energies 500 1000 5000 10000 15000 20000 30000 50000 75000 100000 --clusters --SF $runname/sampling/SF.json --corrections $runname/upstream/corr_params_1d.json
#python runParallel.py --outDir $runname/clusters --nEvt 100 --energies 5000 --clusters --SF $runname/sampling/SF.json --corrections $runname/upstream/corr_params_1d.json


# Compute resolutions and responses
#python compute_resolutions.py --inputDir $runname/clusters --outFile $runname/results.csv --clusters CaloClusters CorrectedCaloClusters CaloTopoClusters CorrectedCaloTopoClusters --MVAcalibCalo $runname/training_calo.json --MVAcalibTopo $runname/training_topo.json

# And make plots
#python plot_resolutions.py --outDir $runname --doFits plot $runname/results.csv --all
#python plot_resolutions.py --outDir $runname --doFits compare clusters CaloClusters CorrectedCaloClusters CaloTopoClusters CorrectedCaloTopoClusters CalibratedCaloClusters CalibratedCaloTopoClusters $runname/results.csv --all
#