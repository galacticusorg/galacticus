# Provides functionality to manage jobs in queue managers (e.g. SLURM).
# Andrew Benson (25-February-2025)
import os
import sys
import time
import re
import json
import lxml.etree as ET
import subprocess

class QueueManager:
    """Base class for queue managers"""
    def __init__(self,name):
        self.name = name

def factory(args):
    """A simple factory to return a queue manager appropriate for the current system."""
    if not os.path.exists(os.environ['GALACTICUS_EXEC_PATH']+"/galacticusConfig.xml"):
        raise Exception("unable to locate the `galacticusConfig.xml` file")
    config      = ET.parse(os.environ['GALACTICUS_EXEC_PATH']+"/galacticusConfig.xml")
    hosts       = config.findall('.//queueManager/host')
    managerType = None
    for host in hosts:
        name = host.find('name')
        match = re.search(name.text,os.environ['HOSTNAME'])
        if match:
            manager = host.find('manager')
            managerType = manager.text
    if managerType is None:
        raise Exception(f"No manager defined for this system")
    hosts = config.findall('.//'+managerType+'/host')
    hostConfig = None
    for host in hosts:
        name = host.find('name')
        match = re.search(name.text,os.environ['HOSTNAME'])
        if match:
            hostConfig = host
    if hostConfig is None:
        raise Exception(f"No config found for manager on this system")
    if managerType == "slurm":
        return SLURMManager(hostConfig,args)
    else:
        raise Exception(f"Unknown queue manager type '{managerType}'")
        
class SLURMManager(QueueManager):
    """A queue manager that interfaces with SLURM."""
    def __init__(self,config,args):
        super().__init__("SLURM")
        self.options = {}
        for option in 'partition',:
            self.options[option]  = config.find(option).text
        for option in 'jobMaximum', 'waitOnSubmit', 'waitOnActive':
            self.options[option]  = int(config.find(option).text)
        if args.partition    is not None:
            self.options['partition'   ] = args.partition
        if args.jobMaximum   is not None:
            self.options['jobMaximum'  ] = args.jobMaximum
        if args.waitOnSubmit is not None:
            self.options['waitOnSubmit'] = args.waitOnSubmit
        if args.waitOnActive is not None:
            self.options['waitOnActive'] = args.waitOnActive
            
    def submitJobs(self,jobs):
        """Submit jobs to a SLURM queue manager and wait for completion."""
        # Define components of SLURM interface,
        activeStates = [ "RUNNING", "PENDING" ]
        optionMap    = {
            "label":           "job-name"       ,
            "partition":       "partition"      , 
            "nodes":           "nodes"          ,
            "tasksPerNode":    "ntasks-per-node",
            "cpusPerTask":     "cpus-per-task"  ,
            "memoryPerThread": "mem-per-cpu"    ,
            "memory":          "mem"            ,
            "walltime":        "time"           ,
            "logOutput":       "output"         ,
            "logError":        "error"
            }
        # Iterate until all jobs are finished.
        activeJobs = {}
        while len(activeJobs) > 0 or len(jobs) > 0:
            # Find all running jobs for this user.
            runningJobs = {}
            squeue      = subprocess.run(['squeue', '--me', '--json'], capture_output=True, text=True)
            if squeue.returncode == 0:
                # Parse the JSON output from squeue
                jobData = json.loads(squeue.stdout)
                for job in jobData['jobs']:
                    if job['job_state'][0] not in activeStates:
                        continue
                    jobID              = str(job['job_id'])
                    jobState           = job['job_state'][0]
                    runningJobs[jobID] = jobState
                    if jobID not in activeJobs:
                        continue
                    if jobState == "RUNNING" and activeJobs[jobID]['state'] != "RUNNING":
                        print(f'Job "{activeJobs[jobID]["label"]}" has started')
                    activeJobs[jobID]['state'] = jobState
                # Remove jobs that are no longer active.
                jobIDsToRemove = []
                for jobID in activeJobs:
                    if jobID not in runningJobs:
                        jobIDsToRemove.append(jobID)
                for jobID in jobIDsToRemove:
                    print(f'Job "{activeJobs[jobID]["label"]}" has finished')
                    if 'onCompletion' in activeJobs[jobID]:
                        activeJobs[jobID]['onCompletion'](activeJobs[jobID])
                    del activeJobs[jobID]
                # Decide if we should submit a new job.
                if len(runningJobs) < self.options['jobMaximum'] and len(jobs) > 0:
                    # Get the next job.
                    job = jobs.pop()
                    # Set job defaults.
                    if "partition" not in job:
                        job['partition'] = self.options['partition']
                    # Determine number of tasks such that we do not exceed the available memory.
                    if "memoryPerThread" in job and not "tasksPerNode" in job:
                        # Determine if we are optimizing for nodes or for threads.
                        optimizeFor = job['optimizeFor'] if "optimizeFor" in job else "threads"
                        # Deterine if we should use the maximum or minimum memory across nodes.
                        useMemory = None
                        if optimizeFor == "threads":
                            useMemory = "maximum"
                        if optimizeFor == "nodes":
                            useMemory = "minimum"
                        # Specify the fraction of memory on a node we will use (i.e. allow some buffer so as not to use all memory).
                        memoryFraction = 0.8
                        # Get info on nodes in this partition.
                        sinfo = subprocess.run(['sinfo', '--partition', job['partition'], '--json'], capture_output=True, text=True)
                        if sinfo.returncode != 0:
                            raise Exception("`sinfo` failed")
                        partitionData = json.loads(sinfo.stdout)
                        job['tasksPerNode'] = 0 if useMemory == "maximum" else 1000000
                        job['cpusPerTask' ] = 1
                        for node in partitionData['sinfo']:
                            countCPUsLimit = min(node['cpus']['maximum'],int(memoryFraction*node['memory'][useMemory]/job['memoryPerThread']))
                            if useMemory == "maximum":
                                if countCPUsLimit > job['tasksPerNode']:
                                    job['tasksPerNode'] = countCPUsLimit
                            else:
                                if countCPUsLimit < job['tasksPerNode']:
                                    job['tasksPerNode'] = countCPUsLimit
                    # Create the job submit file.
                    fileBatch = open(job['launchFile'],"w")
                    fileBatch.write('#!/bin/bash\n')
                    for option in job:
                        if option in optionMap:
                            suffix = "M" if re.match("^mem\-",optionMap[option],) else ""
                            fileBatch.write(f'#SBATCH --{optionMap[option]}={job[option]}{suffix}\n')
                    fileBatch.write(f'ulimit -t unlimited\n')
                    fileBatch.write(f'ulimit -c unlimited\n')
                    if "countOpenMPThreads" in job:
                        fileBatch.write(f'export OMP_NUM_THREADS={job["countOpenMPThreads"]}\n')
                    fileBatch.write(f'{job["command"]}\n')
                    fileBatch.write(f'exit\n')
                    fileBatch.close()
                    print(f'Submitting job "{job["label"]}"')
                    sbatch = subprocess.run(['sbatch', job['launchFile']], capture_output=True, text=True)
                    if sbatch.returncode == 0:
                        # Extract the job ID.
                        match = re.search("\d+", sbatch.stdout)
                        if match:
                            job['jobID'] = match.group(0)
                            job['state'] = "PENDING"
                            activeJobs[job['jobID']] = job
                        else:
                            raise Exception(f"`sbatch` command did not return job ID")
                    else:
                        print(sbatch.stderr)
                        raise Exception(f"`sbatch` command failed with return code: {sbatch.returncode}")
                    time.sleep(self.options['waitOnSubmit'])
                elif len(activeJobs) > 0:
                    time.sleep(self.options['waitOnActive'])
            else:
                # `squeue` failed - report this but do not exit - could be a temporary failure.
                print(f"`squeue` command failed with return code: {squeue.returncode}")
                print(squeue.stderr)
                time.sleep(self.options['waitOnActive'])
