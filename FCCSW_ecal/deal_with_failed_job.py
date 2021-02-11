import os, sys, glob

condor_dir = sys.argv[1]
# 201012_pythia/exec_pdgID_22_pMin_0_pMax_0_thetaMin_0_thetaMax_0_evt_200_jobid_369.sh.err
# "fccsw_output_pythia_ee_Z_ee_jobid_313.root"
rootfile_dir = os.path.join("/eos/user/b/brfranco/rootfile_storage/", condor_dir)

n_failed_jobs = 0
n_jobs = 0
corrupted_rootfiles = []
error_files = []
for errorfile in glob.glob(os.path.join(condor_dir, "*.err")):
    n_jobs += 1
    with open(errorfile, 'r') as read_obj:
        nlines = 0
        for line in read_obj:
            nlines += 1
        if nlines != 15:
            #if "Aborted" in line or "segmentation violation" in line:
                jobid = errorfile.split("_")[-1].split('.')[0]
                print "Job id %s failed"%jobid
                n_failed_jobs += 1
                corrupted_rootfile = glob.glob(os.path.join(rootfile_dir, "*_%s.root"%jobid))
                corrupted_rootfiles.append(corrupted_rootfile[0])
                error_files.append(errorfile)
print "Failed jobs: %d out of %d"%(n_failed_jobs, n_jobs)
print "Rootfile list from failed jobs:"
print " ".join(map(str, corrupted_rootfiles))
for errorfile in error_files:
    print errorfile
