import subprocess
import commands
import sys, os, stat, time, math
import glob
import argparse
from datetime import date

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

def get_condor_submit_header(executable_regex, jobFlavour = 'microcentury'):
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

def get_exec_file_header():# assumes you installed FCCSW locally with the 'install' folder at the root of your FCCSW repository
    return """#!/bin/bash
source %s/../init.sh
source %s/setup.sh
"""%(os.environ.get("FCCSWBASEDIR", ""), os.environ.get("FCCSWBASEDIR", ""))




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-outputFolder", default = "/eos/user/b/brfranco/rootfile_storage/", help = "Output folder absolute path for the rootfile", type = str)
    parser.add_argument("-campaignName", default = date.today().strftime("%y%m%d"), help = "Folder name used to store the submission script, logs, etc, as well as the output rootfile in the outputFolder", type = str)
    parser.add_argument("-gaudiConfig", default = "%s/runCaloSim.py"%os.environ.get("PWD", ""), help = "Absolute path to the gaudi config to use", type = str)
    parser.add_argument("-jobType", default = "caloReco", help = "Tell the type of job we launch. Can be samplingFraction, caloReco, upstreamCorrection", type = str)
    parser.add_argument("-inputFiles", help = "Regex used to get all the input files, if any is needed - not implemented yet.", type = str)
    parser.add_argument("-energies", default = [10], help = "Energies in MeV for the process to generate, behavior depends on fixedEnergy. Make sure you put the energies in ascending order to have an optimal job splitting.", type = int, nargs = '+')
    parser.add_argument("-fixedEnergies", default = True, help = "Do we launch several fixed energies gun or an energy range? If range, energies should have two entries: min and max - not implemented yet", type = bool)
    parser.add_argument("-polarAngles", default = [90], help="Polar angles in degrees for the process to generate, behavior depends on fixedPolarAngles", type = int, nargs = '+')
    parser.add_argument("-fixedPolarAngles", default = True, help = "Do we launch several fixed polar angles gun or a polar angle range? If range, polarAngles should have two entries: min and max", type = bool)
    parser.add_argument("-energiesForDifferentPolarAngles", help = "Use this if you want to generate the additional polar angles only for some energies", type = int, nargs = '+')
    parser.add_argument("-pdgId", default = 22, help = "PDG ID of the particle to shoot", type = int)
    parser.add_argument("-originalNjobs", default = 1, help = "If more than one energy point, the script will submit with increasing number of jobs for increasing energies, tells how many jobs should be used for the first point.", type = int)
    parser.add_argument("-nEvt", default = 1000, help = "How many events to generate.", type = int)
    parser.add_argument("--submit", help="Won't actually submit the jobs unless this is provided.", action = 'store_true')

    args = parser.parse_args()

    # make sure you put the energies in ascending order to have an optimal job splitting
    energies = args.energies  # in MeV
    thetas = args.polarAngles # degrees, will transform to radians later
    if args.energiesForDifferentPolarAngles:
        energies_using_other_thetas = args.energiesForDifferentPolarAngles # may not be interested in having all the theta points for all the energies
    else:
        energies_using_other_thetas = energies
    pdgid = args.pdgId
    original_n_jobs = args.originalNjobs
    total_evt_to_generate = args.nEvt
    gaudi_config_path = args.gaudiConfig
    
    storage_path = args.outputFolder
    campaign_name = args.campaignName
    if not os.path.isdir(campaign_name):
        os.mkdir(campaign_name)
    outfile_storage = os.path.join(storage_path, campaign_name)
    if not os.path.isdir(outfile_storage):
        os.mkdir(outfile_storage)

    if args.jobType == 'samplingFraction':
        command_template = """fccrun %s -n EVT --MomentumMin PMIN --MomentumMax PMAX --ThetaMin THETAMINRADIAN --ThetaMax THETAMAXRADIAN --PdgCodes PDGID --Output "rec DATAFILE='OUTPUTDIR/calibration_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_jobid_JOBID.root' TYP='ROOT' OPT='RECREATE'" --filename OUTPUTDIR/fccsw_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_jobid_JOBID.root"""%(gaudi_config_path)
    elif args.jobType == 'caloReco':
        command_template = """fccrun %s -n EVT --MomentumMin PMIN --MomentumMax PMAX --ThetaMin THETAMINRADIAN --ThetaMax THETAMAXRADIAN --PdgCodes PDGID --filename OUTPUTDIR/fccsw_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_jobid_JOBID.root"""%(gaudi_config_path)
    else:
        print "Wrong jobType provided, read the help to see what is available."
        sys.exit()
    exec_filename_template = os.path.join(campaign_name, "exec_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_evt_EVT_jobid_JOBID.sh")


    total_n_job = 0
    hadd_commands = ""
    for index in range(len(energies)):
        energy = energies[index]
        energy_min = energy
        energy_max = energy
        if index != 0:
            n_jobs = int(n_jobs * math.floor(energy/energies[index-1]))
        else:
            n_jobs = original_n_jobs
        evt_per_job = int(round(total_evt_to_generate/n_jobs))
        if evt_per_job == 0:
            evt_per_job = 1

        if energy in energies_using_other_thetas:
            thetas_for_loop = thetas
        else:
            thetas_for_loop = [thetas[0]]
        for theta in thetas_for_loop:
            print "Treating energy {0}, theta {1}".format(energy, theta)
            theta_min = theta
            theta_max = theta
            job_idx = 0
            evt_already_launched = 0
            if total_evt_to_generate - evt_per_job < 0:
                print "Careful, total_evt_to_generate is smaler than evt_per_job"
                evt_per_job = total_evt_to_generate
            while evt_already_launched <= total_evt_to_generate - evt_per_job:
                command = command_template.replace('EVT', str(evt_per_job)).replace('PMIN', str(energy_min)).replace('PMAX', str(energy_max)).replace('THETAMINRADIAN', str(math.radians(theta_min))).replace('THETAMAXRADIAN', str(math.radians(theta_max))).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid)).replace('THETAMIN', str(theta_min)).replace('THETAMAX', str(theta_max)).replace('JOBID', str(job_idx))
                exec_filename = exec_filename_template.replace('EVT', str(evt_per_job)).replace('PMIN', str(energy)).replace('PMAX', str(energy)).replace('THETAMIN', str(theta)).replace('THETAMAX', str(theta)).replace('JOBID', str(job_idx)).replace('PDGID', str(pdgid))
                with open(exec_filename, "w") as f:
                    f.write(get_exec_file_header())
                    f.write(command)
                st = os.stat(exec_filename)
                os.chmod(exec_filename, st.st_mode | stat.S_IEXEC)
                evt_already_launched += evt_per_job
                job_idx += 1
                total_n_job += 1

            evt_last_job = total_evt_to_generate - evt_already_launched
            if evt_last_job > 0:
                command = command_template.replace('EVT', str(evt_last_job)).replace('PMIN', str(energy_min)).replace('PMAX', str(energy_max)).replace('THETAMINRADIAN', str(math.radians(theta_min))).replace('THETAMAXRADIAN', str(math.radians(theta_max))).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid)).replace('THETAMIN', str(theta_min)).replace('THETAMAX', str(theta_max)).replace('JOBID', str(job_idx))
                exec_filename = exec_filename_template.replace('EVT', str(evt_last_job)).replace('PMIN', str(energy)).replace('PMAX', str(energy)).replace('THETAMIN', str(theta)).replace('THETAMAX', str(theta)).replace('JOBID', str(job_idx)).replace('PDGID', str(pdgid))
                with open(exec_filename, "w") as f:
                    f.write(get_exec_file_header())
                    f.write(command)
                st = os.stat(exec_filename)
                os.chmod(exec_filename, st.st_mode | stat.S_IEXEC)
                job_idx += 1
                total_n_job += 1

            if args.jobType == 'samplingFraction':
                hadd_commands += "rm OUTPUTDIR/fccsw_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_jobid_*.root\n".replace('PMIN', str(energy_min)).replace('PMAX', str(energy_max)).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid)).replace('THETAMIN', str(theta_min)).replace('THETAMAX', str(theta_max))
                hadd_commands += "hadd OUTPUTDIR/calibration_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX.root OUTPUTDIR/calibration_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_jobid_*.root\n".replace('PMIN', str(energy_min)).replace('PMAX', str(energy_max)).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid)).replace('THETAMIN', str(theta_min)).replace('THETAMAX', str(theta_max))
                hadd_commands += "#rm OUTPUTDIR/calibration_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_jobid_*.root\n".replace('PMIN', str(energy_min)).replace('PMAX', str(energy_max)).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid)).replace('THETAMIN', str(theta_min)).replace('THETAMAX', str(theta_max))
            else:
                hadd_commands += "hadd  OUTPUTDIR/fccsw_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX.root OUTPUTDIR/fccsw_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_jobid_*.root\n".replace('PMIN', str(energy_min)).replace('PMAX', str(energy_max)).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid)).replace('THETAMIN', str(theta_min)).replace('THETAMAX', str(theta_max))
                hadd_commands += "#rm OUTPUTDIR/fccsw_output_pdgID_PDGID_pMin_PMIN_pMax_PMAX_thetaMin_THETAMIN_thetaMax_THETAMAX_jobid_*.root\n".replace('PMIN', str(energy_min)).replace('PMAX', str(energy_max)).replace('OUTPUTDIR', outfile_storage).replace('PDGID', str(pdgid)).replace('THETAMIN', str(theta_min)).replace('THETAMAX', str(theta_max))

    # write the hadd script
    hadd_script_path = os.path.join(campaign_name, "hadd.sh")
    with open(hadd_script_path, "w") as f:
        f.write(hadd_commands)
    st = os.stat(hadd_script_path)
    os.chmod(hadd_script_path, st.st_mode | stat.S_IEXEC)

    # write the condor submit file
    condor_submit_path = campaign_name + ".sub"
    with open(condor_submit_path, "w") as f:
        f.write(get_condor_submit_header(exec_filename_template.replace('EVT', '*').replace('PMIN', '*').replace('PMAX', '*').replace('THETAMIN', '*').replace('THETAMAX', '*').replace('JOBID', '*').replace('PDGID', '*')))
    submit_cmd = "condor_submit %s"%condor_submit_path
    print "%d jobs prepared in %s"%(total_n_job, campaign_name)
    print (submit_cmd)
    if args.submit:
        job=SubmitToCondor(submit_cmd, 10)
    print "Will write the output in %s"%outfile_storage


