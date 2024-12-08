#!/bin/bash

# Define the environment setup script
setup_script="/cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh"
# Define the energies to use
energies=("10" "15")

# Define the job submission function
function create_condor_submit_file {
    local energy=$1
    local submit_file=$2
    local executable_name="run_simulation_${energy}GeV.sh"

    # Create the Condor submit file
    cat > $submit_file <<EOF
executable     = $executable_name
Log            = \$(filename)/condor.log
Output         = \$(filename)/condor.out
Error          = \$(filename)/condor.err
requirements    = ( (OpSysAndVer =?= "AlmaLinux9") && (Machine =!= LastRemoteHost) && (TARGET.has_avx2 =?= True) )
max_retries    = 3
+JobFlavour    = "longlunch"
RequestCpus    = 1
queue filename matching files $executable_name
EOF
}

# Define the job submission function
function create_simulation_script {
    local energy=$1
    local script_name="run_simulation_${energy}GeV.sh"

    # Create the simulation script
    cat > $script_name <<EOF
#!/bin/sh

# Setup the environment
source $setup_script

# Run the ddsim command
ddsim --enableGun --gun.distribution uniform --gun.energy "${energy}*GeV" \
--gun.particle pi- --gun.thetaMin "68*degree" --gun.thetaMax "68*degree" \
--numberOfEvents 10 --outputFile electron_gun_${energy}GeV_ALLEGRO_SIM.root \
--random.enableEventSeed --random.seed 42 \
--compactFile \$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
EOF

    # Make the script executable
    chmod +x $script_name
}

# Submit jobs for each energy
for energy in "${energies[@]}"; do
    # Create the simulation script
    create_simulation_script $energy

    # Create the Condor submit file
    submit_file="submit_${energy}GeV.submit"
    create_condor_submit_file $energy $submit_file

    # Submit the job
    condor_submit $submit_file
done
