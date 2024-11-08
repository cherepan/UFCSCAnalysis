#!/usr/bin/env python3


import subprocess

def submit_condor_job(root_file, python_script, job_name, config):
    submit_file_content = """
universe     = vanilla
rank         = memory
executable   = run.py
output       = Run-Condor_{job_name}.o
error        = Run-Condor_{job_name}.e
log          = Run-Condor_{job_name}.log
getenv       = True
+MaxRuntime  = 9600
notification = Error

# Specify arguments
arguments    = -f {root_file} -p {python_script} -j {job_name} -c {config}

queue 1
""".format(root_file=root_file, python_script=python_script, job_name=job_name, config=config)

    # Write the submit file
    submit_file = "submit_{job_name}.condor".format(job_name=job_name)
    with open(submit_file, 'w') as f:
        f.write(submit_file_content)

    # Submit the job to Condor
    subprocess.call(["condor_submit", submit_file])

    print("Submitted Condor job with job_name: {}".format(job_name))


def main():
    # Example job parameters
    job_params = [

        ##  config 0 - clean; 1 - noise was used in the past
        {"root_file": "../UFCSCRootMaker/CSC_UF_Ntuple_SegmentAlgoUF.root",      "python_script": "segments_efficiency.py", "job_name": "UFSegmentsEfficiency", "config": "0"},
        {"root_file": "../UFCSCRootMaker/CSC_UF_Ntuple_SegmentAlgoDefault.root", "python_script": "segments_efficiency.py", "job_name": "DFSegmentsEfficiency", "config": "0"},
        
        
        # Add more job parameters as needed
    ]

    # Submit Condor jobs for each set of parameters
    for params in job_params:
        submit_condor_job(params["root_file"], params["python_script"], params["job_name"], params["config"])

if __name__ == "__main__":
    main()




