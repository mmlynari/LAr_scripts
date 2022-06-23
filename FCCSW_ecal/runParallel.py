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
        self.extra_args = ""
        os.makedirs(outdir, exist_ok=True)

    def process(self, nevt, energy, jobId):
        return subprocess.run([f"fccrun {self.script} -n {nevt} --MomentumMin {energy} --MomentumMax {energy} \
            --filename {self.outdir}/{self.output_tag}_energy_{energy}_jobid_{jobId}.root --seedValue {jobId} {self.extra_args}"], shell=True)

    def hadd(self, energy):
        filestub = f"{self.outdir}/{self.output_tag}_energy_{energy}"
        return subprocess.run([f"hadd -f {filestub}.root {filestub}_jobid_*.root"], shell=True)

    def rm(self, energy):
        return subprocess.run([f"rm -f {self.outdir}/{self.output_tag}_energy_{energy}_jobid_*.root"], shell=True)


class SamplingJobProcessor(JobProcessor):
    def __init__(self, outdir):
        script = "fcc_ee_samplingFraction_inclinedEcal.py"
        output_tag = "sampling_output"
        super().__init__(script, outdir, output_tag)

    def postprocess(self, energy):
        energy_GeV = energy/1.e3
        subprocess.run([f"python FCC_calo_analysis_cpp/plot_samplingFraction.py
            {self.outdir}/{self.output_tag}_energy_{energy}.root {energy_GeV}  --totalNumLayers 12
            --preview -outputfolder {self.outdir} --json {self.outdir}/SF.json"], shell=True)
        return 0

class UpstreamJobProcessor(JobProcessor):
    def __init__(self, outdir, sampling_fracs=None, postprocess_scripts_dir=None):
        script = "fcc_ee_upstream_inclinedEcal.py"
        output_tag = "upstream_output"
        super().__init__(script, outdir, output_tag)
        if sampling_fracs:
            self.extra_args += "--samplingFractions "
            self.extra_args += ' '.join([str(s) for s in sampling_fracs])
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
    def __init__(self, outdir, sampling_fracs=None, corrections=None):
        script = "runTopoAndSlidingWindowAndCaloSim.py"
        output_tag = "upstream_output"
        if sampling_fracs:
            self.extra_args += "--samplingFraction "
            self.extra_args += ' '.join([str(s) for s in sampling_fracs])
        if corrections:
            self.extra_args += "--upstreamParameters "
            self.extra_args += ' '.join([str(s) for s in corrections['up']])
            self.extra_args += "--downstreamParameters "
            self.extra_args += ' '.join([str(s) for s in corrections['do']])
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--outDir', default='./', type=str, help='output directory for plots')
    parser.add_argument('--nEvt', default=1000, type=int, help='number of events to process per point')
    parser.add_argument('--energies', default=[10000], action='extend', nargs='+', type=int,
            help='energies to process')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--sampling', action='store_true', help='compute the sampling fractions')
    group.add_argument('--upstream', action='store_true', help='compute the upstream corrections')
    group.add_argument('--clusters', action='store_true', help='run fixed size and topo clusterings')
    parser.add_argument('--SF', default='', type=str, help='JSON file containing sampling fractions')
    parser.add_argument('--corrections', default='', type=str, help='JSON file containing upstream and downstream corrections')
    args = parser.parse_args()

    sampling_fracs = None
    if args.SF:
        with open(args.SF, 'r') as jsonfile:
            sampling_fracs = json.load(jsonfile)
            print("Applying sampling fractions:", sampling_fracs)

    corrections = None
    if args.corrections:
        with open(args.corrections, 'r') as jsonfile:
            source = json.load(jsonfile)
            dict_up = {e['name'] : e for e in source['corr_params'] if e['type']=='upstream'}
            dict_do = {e['name'] : e for e in source['corr_params'] if e['type']=='downstream'}
            corrections = {}
            corrections['up'] = [dict_up['a'], dict_up['b'], dict_up['c'], dict_up['d'], dict_up['e']]
            corrections['do'] = [dict_do['a'], dict_do['b'], dict_do['c'], dict_do['d'], dict_do['e']]

    #energies = [500, 1000, 5000, 10000, 15000, 20000, 30000, 50000, 75000, 100000]
    if args.sampling:
        samJobPr = SamplingJobProcessor(outDir)
        run_the_jobs(samJobPr, energies, nEvt, True, True)
    if args.upstream:
        upJobPr = UpstreamJobProcessor(outDir, sampling_fracs=sampling_fracs)
        run_the_jobs(upJobPr, energies, nEvt, True, True)
    elif args.clusters:
        clJobPr = ClusterJobProcessor(outDir, sampling_fracs=sampling_fracs, corrections=corrections)
        run_the_jobs(clJobPr, energies, nEvt, True, False)


if __name__ == "__main__":
    main()
