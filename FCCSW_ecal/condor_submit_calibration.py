import subprocess
import commands
import sys, os, stat, time, math
import glob

#__________________________________________________________
def getCommandOutput(command):
    p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE,universal_newlines=True)
    (stdout,stderr) = p.communicate()
    return {"stdout":stdout, "stderr":stderr, "returncode":p.returncode}

#__________________________________________________________
def SubmitToCondor(cmd,nbtrials):
    submissionStatus=0
    cmd=cmd.replace('//','/') # -> dav : is it needed?
    for i in range(nbtrials):
        outputCMD = getCommandOutput(cmd)
        stderr=outputCMD["stderr"].split('\n')
        stdout=outputCMD["stdout"].split('\n') # -> dav : is it needed?

        if len(stderr)==1 and stderr[0]=='' :
            print ("------------GOOD SUB ")
            submissionStatus=1
        else:
            print ("++++++++++++ERROR submitting, will retry")
            print ("Trial : "+str(i)+" / "+str(nbtrials))
            print ("stderr : ",len(stderr))
            print (stderr)

            time.sleep(10)


        if submissionStatus==1:
            return 1

        if i==nbtrials-1:
            print ("failed sumbmitting after: "+str(nbtrials)+" trials, will exit")
            return 0

def get_condor_submit_header(executable_regex):
    return """executable     = $(filename)
Log            = $(filename).log
Output         = $(filename).out
Error          = $(filename).err
#getenv         = True
requirements    = ( (OpSysAndVer =?= "CentOS7") && (Machine =!= LastRemoteHost) && (TARGET.has_avx2 =?= True) )
#on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
max_retries    = 3
+JobFlavour    = "longlunch"
RequestCpus = 1
queue filename matching files {0}
""".format(executable_regex)

def get_exec_file_header():# assumes you installed FCCSW locally with the 'install' folder at the root of your FCCSW repository
    return """#!/bin/bash
source %s/../init.sh
source %s/setup.sh
"""%(os.environ.get("FCCSWBASEDIR", ""), os.environ.get("FCCSWBASEDIR", ""))




if __name__ == "__main__":
    storage_path = sys.argv[1]
    calibration_campaign_name = sys.argv[2] #'20201118_condor_calib_5kEvt' #sys.argv[1]
    if not os.path.isdir(calibration_campaign_name):
        os.mkdir(calibration_campaign_name)
    outfile_storage = os.path.join(storage_path, calibration_campaign_name)
    if not os.path.isdir(outfile_storage):
        os.mkdir(outfile_storage)

    command_template = """fccrun %s/fcc_ee_samplingFraction_inclinedEcal.py -n EVT --MomentumMin PMIN --MomentumMax PMAX --ThetaMin THETAMIN --ThetaMax THETAMAX --PdgCodes PDGID --Output "rec DATAFILE='OUTPUTDIR/calibration_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_pihalved_thetaMax_pihalved.root' TYP='ROOT' OPT='RECREATE'" --filename OUTPUTDIR/fccsw_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_pihalved_thetaMax_pihalved.root"""%(os.environ.get("PWD", ""))
    exec_filename_template = os.path.join(calibration_campaign_name, "exec_pdgID_PDGID_pMin_PMIN_pMax_PMAX_evt_EVT_jobid_JOBID.sh")

    # make sure you put the energies in ascending order to have an optimal job splitting
    #energies = [300, 1000, 10000, 50000, 100000] # in MeV
    energies = [300, 10000] # in MeV
    theta = 1.57079
    pdgid = 22
    original_n_jobs = 1
    total_evt_to_generate = 100
    #total_evt_to_generate = 5000

    total_n_job = 0
    for index in range(len(energies)):
        energy = energies[index]
        job_idx = 0
        if index != 0:
            n_jobs = int(n_jobs * math.floor(energy/energies[index-1]))
        else:
            n_jobs = original_n_jobs
        evt_per_job = int(round(total_evt_to_generate/n_jobs))

        evt_already_launched = 0
        while evt_already_launched < total_evt_to_generate - evt_per_job:
            command = command_template.replace('EVT', str(evt_per_job)).replace('PMIN', str(energy)).replace('PMAX', str(energy)).replace('THETAMIN', str(theta)).replace('THETAMAX', str(theta)).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid))
            exec_filename = exec_filename_template.replace('EVT', str(evt_per_job)).replace('PMIN', str(energy)).replace('PMAX', str(energy)).replace('JOBID', str(job_idx)).replace('PDGID', str(pdgid))
            with open(exec_filename, "w") as f:
                f.write(get_exec_file_header())
                f.write(command)
            st = os.stat(exec_filename)
            os.chmod(exec_filename, st.st_mode | stat.S_IEXEC)
            evt_already_launched += evt_per_job
            job_idx += 1
            total_n_job += 1

        evt_last_job = total_evt_to_generate - evt_already_launched
        if evt_last_job >= 0:
            command = command_template.replace('EVT', str(evt_last_job)).replace('PMIN', str(energy)).replace('PMAX', str(energy)).replace('THETAMIN', str(theta)).replace('THETAMAX', str(theta)).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid))
            exec_filename = exec_filename_template.replace('EVT', str(evt_per_job)).replace('PMIN', str(energy)).replace('PMAX', str(energy)).replace('JOBID', str(job_idx)).replace('PDGID', str(pdgid))
            with open(exec_filename, "w") as f:
                f.write(get_exec_file_header())
                f.write(command)
            st = os.stat(exec_filename)
            os.chmod(exec_filename, st.st_mode | stat.S_IEXEC)
            job_idx += 1
            total_n_job += 1

    # write the condor submit file
    condor_submit_path = calibration_campaign_name + ".sub"
    with open(condor_submit_path, "w") as f:
        f.write(get_condor_submit_header(exec_filename_template.replace('EVT', '*').replace('PMIN', '*').replace('PMAX', '*').replace('JOBID', '*').replace('PDGID', '*')))
    submit_cmd = "condor_submit %s"%condor_submit_path
    print "%d jobs prepared in %s"%(total_n_job, calibration_campaign_name)
    print (submit_cmd)
    job=SubmitToCondor(submit_cmd, 10)
    print "Will write the output in %s"%outfile_storage







