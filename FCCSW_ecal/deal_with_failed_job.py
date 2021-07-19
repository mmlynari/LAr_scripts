import os, sys, glob, re

condor_dir = sys.argv[1]
if '_calib_' in condor_dir:
    n_line_errfile = 1
elif '_upstream' in condor_dir:
    n_line_errfile = 2
else:
    n_line_errfile = 14
# 201012_pythia/exec_pdgID_22_pMin_0_pMax_0_thetaMin_0_thetaMax_0_evt_200_jobid_369.sh.err
# "fccsw_output_pythia_ee_Z_ee_jobid_313.root"
rootfile_dir = os.path.join("/eos/user/b/brfranco/rootfile_storage/", condor_dir)

fixable_patterns = ["WriteBasketImpl"]
fixable_shs = []

dict_energy_njob_nfailedjob = {}
n_failed_jobs = 0
n_jobs = 0
corrupted_rootfiles = []
error_files = []
for errorfile in glob.glob(os.path.join(condor_dir, "*.err")):
    energy = errorfile.split('_pMin_')[1].split("_")[0]
    if energy in list(dict_energy_njob_nfailedjob.keys()):
        dict_energy_njob_nfailedjob[energy][0] += 1
    else:
        dict_energy_njob_nfailedjob[energy] = [1, 0]
    fixable = False
    n_jobs += 1
    with open(errorfile, 'r') as read_obj:
        nlines = 0
        for line in read_obj:
            nlines += 1
            for fixable_pattern in fixable_patterns:
                if fixable_pattern in line:
                    fixable = True
        if nlines != n_line_errfile:
                dict_energy_njob_nfailedjob[energy][1] += 1
            #if "Aborted" in line or "segmentation violation" in line:
                jobid = errorfile.split("_")[-1].split('.')[0]
                #print "Job %s failed"%errorfile
                n_failed_jobs += 1
                if fixable:
                    fixable_shs.append(errorfile.replace(".sh.err", ".sh"))
                if 'pythia' in condor_dir:
                    corrupted_rootfile = glob.glob(os.path.join(rootfile_dir, "*_%s.root"%jobid))
                elif '_upstream' in condor_dir:
                    pMin = errorfile.split('pMin_')[1].split('_')[0]
                    corrupted_rootfile = glob.glob(os.path.join(rootfile_dir, "fccsw_upstream_output_pMin_%s_pMax_%s_jobid_%s.root"%(pMin, pMin, jobid)))
                else:
                    file_pattern = errorfile.split('/')[-1]
                    file_pattern = re.sub(r'exec.*_pdgID', '*pdgID', file_pattern).replace(".sh.err", ".root")
                    corrupted_rootfile = glob.glob(os.path.join(rootfile_dir, file_pattern))
                if len(corrupted_rootfile) == 0:
                    print("Already removed %s"%(os.path.join(rootfile_dir, file_pattern)))
                else:
                    corrupted_rootfiles.append(corrupted_rootfile[0])
                error_files.append(errorfile)
print("Failed jobs: %d out of %d"%(n_failed_jobs, n_jobs))
print("Rootfile list from failed jobs:")
print(" ".join(map(str, corrupted_rootfiles)))
print("Error logs:")
for errorfile in error_files:
    print(errorfile)
print("Failed job that could be recovered:")
for fixable_sh in fixable_shs:
    print(fixable_sh)
print("Number of job, failed job and percentage per energy:")
for energy in list(dict_energy_njob_nfailedjob.keys()):
    dict_energy_njob_nfailedjob[energy].append(round(100*dict_energy_njob_nfailedjob[energy][1]/float(dict_energy_njob_nfailedjob[energy][0]), 1))
print(dict_energy_njob_nfailedjob)
