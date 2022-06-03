#!/usr/bin/env python

import os
import os.path
import subprocess
import multiprocessing as mp
import shutil
import sys

class JobProcessor:
    def __init__(self, script, outdir, output_tag):
        self.script=script
        self.outdir=outdir
        self.output_tag = output_tag
        os.makedirs(outdir, exist_ok=True)

    def process(self, nevt, energy, jobId):
        return subprocess.run([f"fccrun {self.script} -n {nevt} --MomentumMin {energy} --MomentumMax {energy} \
            --filename {self.outdir}/{self.output_tag}_energy_{energy}_jobid_{jobId}.root --seedValue {jobId}"], shell=True)

    def hadd(self, energy):
        filestub = f"{self.outdir}/{self.output_tag}_energy_{energy}"
        return subprocess.run([f"hadd {filestub}.root {filestub}_jobid_*.root"], shell=True)

    def rm(self, energy):
        return subprocess.run([f"rm -f {self.outdir}/{self.output_tag}_energy_{energy}_jobid_*.root"], shell=True)


class UpstreamJobProcessor(JobProcessor):
    def __init__(self, outdir, postprocess_scripts_dir=None):
        script = "fcc_ee_upstream_inclinedEcal.py"
        output_tag = "upstream_output"
        super().__init__(script, outdir, output_tag)
        self.postprocess_dir = os.path.join(os.getenv("FCCBASEDIR"), 'k4SimGeant4/Detector/DetStudies/scripts/')

    def postprocess(self, energy):
        subprocess.run([f"{self.postprocess_dir}/cec_process_events -i {self.outdir}/{self.output_tag}_energy_{energy}.root -t upstream --plot-file-format png \
        -o {self.outdir}/fit_results.json --plot-directory {self.outdir} --func-from 0.005 --func-to 0.06"], shell=True)

        subprocess.run([f"{self.postprocess_dir}/cec_derive1 -i {self.outdir}/fit_results.json -t upstream --plot-file-format png \
        --plot-directory {self.outdir} --functions '[0]+[1]/(x-[2])' '[0]+[1]/(x-[2])'"], shell=True)

        subprocess.run([f"{self.postprocess_dir}/cec_process_events -i {self.outdir}/{self.output_tag}_energy_{energy}.root -t downstream --function pol2  \
        --plot-file-format png -o {self.outdir}/fit_results.json --plot-directory {self.outdir} --func-from 0.005 --func-to 10"], shell=True)

        subprocess.run([f"{self.postprocess_dir}/cec_derive1 -i {self.outdir}/fit_results.json -t downstream --plot-file-format png \
        --plot-directory {self.outdir} --functions '[0]+[1]*x' '[0]+[1]/sqrt(x)' '[0]+[1]/x'"], shell=True)
        return 0


class ClusterJobProcessor(JobProcessor):
    def __init__(self, outdir, postprocess_scripts_dir=None):
        script = "runTopoAndSlidingWindowAndCaloSim.py"
        output_tag = "upstream_output"
        super().__init__(script, outdir, output_tag)


def run_the_jobs(jobProcessor, energies, nEvt, do_process, do_postprocess):
    with mp.Pool() as p:
        if do_process:
            nEvtMaxPerJob = [int(3000*10000/e) for e in energies] # 3000 evts for 10GeV, with linear scaling
            nEvtPerJob = [min(nEvt, nMax) for nMax in nEvtMaxPerJob]
            args = []
            jobId = 1
            for e, nEvtPerJobForE in zip(energies, nEvtMaxPerJob):
                nEvtToLaunch = nEvt
                while nEvtToLaunch>0:
                    nLaunched = min(nEvtToLaunch, nEvtPerJobForE)
                    args.append((nLaunched, e, jobId))
                    jobId += 1
                    nEvtToLaunch -= nLaunched

            print("About to send jobs with parameters:")
            for a in args:
                print(a)
            res = p.starmap_async(jobProcessor.process, args)
            print("Retcodes of the jobs:")
            print (res.get())
            print()

            print("Hadd'ing results:")
            res = p.map_async(jobProcessor.hadd, energies)
            print("Retcodes of the jobs:")
            print (res.get())
            print()

            print("Removing intermediate files:")
            res = p.map_async(jobProcessor.rm, energies)
            print("Retcodes of the jobs:")
            print (res.get())
            print()

        if do_postprocess:
            # Do it sequentially
            res = p.map(jobProcessor.postprocess, energies, chunksize=50)



def main():
    energies = [500, 1000, 5000, 10000, 15000, 20000, 30000, 50000, 75000, 100000]
    #energies = [1000]
    nEvt = 10000
    outdir = "baseline_LAr_testResol_1"
    #upJobPr = UpstreamJobProcessor(outdir)
    #run_the_jobs(upJobPr, energies, nEvt, True, True)
    clJobPr = ClusterJobProcessor(outdir)
    run_the_jobs(clJobPr, energies, nEvt, True, False)

    #return


if __name__ == "__main__":
    main()
