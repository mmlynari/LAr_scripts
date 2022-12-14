import os, sys, glob, re
# usage: python deal_with_failed_job.py condorFolderName
# Will print information about failed jobs and the location of the rootfile to be removed before to do the hadd (you can of course also resubmit the failed jobs)

condor_dir = sys.argv[1]

# put here the number of line that go into the .err stream for a succesful job (it is unfortunately not  0)
if '_calib_' in condor_dir:
    n_line_errfile = 1
elif '_upstream' in condor_dir:
    n_line_errfile = 2
else:
    n_line_errfile = 1
# 201012_pythia/exec_pdgID_22_pMin_0_pMax_0_thetaMin_0_thetaMax_0_evt_200_jobid_369.sh.err
# "fccsw_output_pythia_ee_Z_ee_jobid_313.root"
rootfile_dir = os.path.join("/eos/user/b/brfranco/rootfile_storage/", condor_dir)

fixable_patterns = ["WriteBasketImpl", "error flushing file"]
logFileErrorPatterns = ["wall time exceeded"]
fixable_shs = []

dict_energy_njob_nfailedjob = {}
n_failed_jobs = 0
n_failed_jobs_logs = 0
n_jobs = 0
corrupted_rootfiles = []
error_files = []
for jobfile in glob.glob(os.path.join(condor_dir, "*jobid_*.sh")):
    n_jobs += 1
    failedDueToLog = False
    fileNotFound = False
    errorfile = jobfile + ".err"
    logfile = jobfile + ".log"
    try:
        with open(logfile, 'r') as read_obj:
            for line in read_obj:
                for pattern in logFileErrorPatterns:
                    #print(pattern, " ", line)
                    if pattern in line:
                        failedDueToLog = True
                        n_failed_jobs_logs += 1
    except IOError:
        fileNotFound = True
    energy = errorfile.split('_pMin_')[1].split("_")[0]
    if energy in list(dict_energy_njob_nfailedjob.keys()):
        dict_energy_njob_nfailedjob[energy][0] += 1
    else:
        dict_energy_njob_nfailedjob[energy] = [1, 0]
    fixable = False
    nlines = 0
    try:
        with open(errorfile, 'r') as read_obj:
            for line in read_obj:
                nlines += 1
                for fixable_pattern in fixable_patterns:
                    if fixable_pattern in line:
                        fixable = True
    except IOError:
        fileNotFound = True

    if nlines != n_line_errfile or failedDueToLog or fileNotFound:
            dict_energy_njob_nfailedjob[energy][1] += 1
        #if "Aborted" in line or "segmentation violation" in line:
            jobid = errorfile.split("_")[-1].split('.')[0]
            if failedDueToLog:
                print("Job %s failed but only visible in the log file (probably wall time exceeded and condor killed the job --> change the queue 'microcentury', 'longlunch', 'workday', etc)"%(logfile))
            elif fileNotFound:
                print("Warning: one the file (.log or .err) not found for job %s"%jobfile)
            else:
                print("Job %s failed with messages in the .err file, Can be recovered? %s"%(errorfile, str(fixable)))
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
print("Failed jobs due to wall time: %d"%(n_failed_jobs_logs))
print("Rootfile list from failed jobs:")
print(" ".join(map(str, corrupted_rootfiles)))
print("Error logs:")
#for errorfile in error_files:
#    print(errorfile)
print("Failed job that could be recovered:")
for fixable_sh in fixable_shs:
    print(fixable_sh)
print("Number of job, failed job and percentage per energy:")
for energy in list(dict_energy_njob_nfailedjob.keys()):
    dict_energy_njob_nfailedjob[energy].append(round(100*dict_energy_njob_nfailedjob[energy][1]/float(dict_energy_njob_nfailedjob[energy][0]), 1))
print(dict_energy_njob_nfailedjob)
