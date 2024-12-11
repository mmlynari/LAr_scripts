#!/bin/bash

# Define the environment setup script
setup_script="/cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh -r 2024-12-11"

# Define the energies and thetas to use
energies=("2" "3" "5" "10" "20" "30" "50" "100" "120")
thetas=("68" "40")
events="10000"

# Define output folder
outputfolder="/eos/user/m/mmlynari/FCC_fellow/FCC_rootfile_storage/HCal_v11Dec24/pionSIM"

# Function to create a Condor submit file
function create_condor_submit_file {
    local energy=$1
    local theta=$2
    local submit_file=$3
    local executable_name="run_simulation_${energy}GeV_${theta}deg.sh"
    local basename=$(basename $executable_name .sh)  # Extract the base name without the extension

    # Create the Condor submit file
    cat > $submit_file <<EOF
executable     = $executable_name
Log            = ${basename}.log
Output         = ${basename}.out
Error          = ${basename}.err
requirements    = ( (OpSysAndVer =?= "AlmaLinux9") && (Machine =!= LastRemoteHost) && (TARGET.has_avx2 =?= True) )
max_retries    = 3
+JobFlavour    = "tomorrow"
RequestCpus    = 1
queue
EOF
}

# Function to create a simulation script
function create_simulation_script {
    local energy=$1
    local theta=$2
    local events=$3
    local script_name="run_simulation_${energy}GeV_${theta}deg.sh"

    # Create the simulation script
    cat > $script_name <<EOF
#!/bin/sh

# Setup the environment
source $setup_script

# Run the ddsim command
ddsim --enableGun --gun.distribution uniform --gun.energy "${energy}*GeV" \
--gun.particle pi- --gun.thetaMin ${theta}*degree --gun.thetaMax ${theta}*degree \
--numberOfEvents $events --outputFile ${outputfolder}/pi_gun_pMin_${energy}_theta_${theta}_events${events}_ALLEGRO_SIM.root \
--random.enableEventSeed --random.seed 42 \
--compactFile /afs/cern.ch/user/m/mmlynari/workspace//FCC_11Dec24/k4geo/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03_tileStandalone.xml
EOF

    # Make the script executable
    chmod +x $script_name
}

# Submit jobs for each energy and theta combination
for energy in "${energies[@]}"; do
    for theta in "${thetas[@]}"; do
        # Create the simulation script
        create_simulation_script $energy $theta $events

        # Create the Condor submit file
        submit_file="submit_${energy}GeV_${theta}deg.submit"
        create_condor_submit_file $energy $theta $submit_file

        # Submit the job
        condor_submit $submit_file
    done
done
