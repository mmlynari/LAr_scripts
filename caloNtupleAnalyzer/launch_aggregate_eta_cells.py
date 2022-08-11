import os
from datetime import date

def get_condor_submit_header(executable_regex, jobFlavour = 'longlunch'):
    return """executable     = $(filename)
Log            = $(filename).log
Output         = $(filename).out
Error          = $(filename).err
requirements    = ( (OpSysAndVer =?= "CentOS7") && (Machine =!= LastRemoteHost) && (TARGET.has_avx2 =?= True) )
max_retries    = 3
+JobFlavour    = "{0}"
RequestCpus = 1
queue filename matching files {1}
""".format(jobFlavour, executable_regex)

exec_file_template = """#!/bin/bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh
#source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
python PATHTOAGGREGATECODE/aggregate_eta_cells.py -inputFile INPUTFILE -postfix POSTFIX -startEvt STARTEVT -endEvt ENDEVT
"""


pi0 = False
#condor_dir = "condor_dir_eta_aggregation_" + date.today().strftime("%y%m%d")
if pi0:
    condor_dir = "condor_dir_eta_aggregation_pi0"
    inputFile = "/eos/user/b/brfranco/rootfile_storage/fcc_analysis_ouput/220808_pi0_flat_1_100_noNoise_caloReco/fccsw_output_pdgID_111_pMin_1000_pMax_100000_thetaMin_50_thetaMax_130.root"
    condor_submit_file = "submit_aggregateCell_pi0.cmd"
else:
    condor_dir = "condor_dir_eta_aggregation_gamma"
    inputFile = "/eos/user/b/brfranco/rootfile_storage/fcc_analysis_ouput/220810_gamma_flat_1_100_noNoise_caloReco/fccsw_output_pdgID_22_pMin_1000_pMax_100000_thetaMin_50_thetaMax_130.root"  
    condor_submit_file = "submit_aggregateCell_gamma.cmd"

if not os.path.isdir(condor_dir):
    os.mkdir(condor_dir)

exec_file_name_path_template = "run_aggregate_cells_ID.sh"
evt_btach = 500
evt_tot = 100000
for jobID in range(int(evt_tot/float(evt_btach))):
    exec_file_name_path = os.path.join(condor_dir, exec_file_name_path_template.replace("ID", str(jobID)))
    with open(exec_file_name_path, "w") as f:
        f.write( exec_file_template.replace("PATHTOAGGREGATECODE", os.environ.get("PWD", "")).replace("INPUTFILE", inputFile).replace("POSTFIX", str(jobID)).replace("STARTEVT", str(jobID * evt_btach)).replace("ENDEVT", str((jobID + 1) * evt_btach - 1)) )

with open(condor_submit_file, "w") as f:
    f.write(get_condor_submit_header(os.path.join(os.environ.get("PWD", ""), condor_dir, exec_file_name_path_template.replace("ID", "*"))))
