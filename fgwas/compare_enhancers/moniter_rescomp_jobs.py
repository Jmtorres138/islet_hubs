#!/usr/bin/env python

import re
import sys
import time

from optparse import OptionParser
from subprocess import Popen, PIPE
from subprocess import check_call
#import subprocess as sp

usage = '%prog [options] <job id> [<job id>...]'
desc = 'Poll qstat and wait until all <job id>s are finished'
parser = OptionParser(usage=usage,description=desc)

array_job_match = '^(\d+)\[\]\.(.*)'
array_job_regex = '^%s\[[0-9]\+\]'

def get_job_ids(jobid_regex):
    # Jason Matthew Torres
    #qstat_p = Popen('qstat -t | grep "%s" | cut -f 4 -d " "'%jobid_regex,shell=True,stdout=PIPE)
    qstat_p = Popen('qstat -t | grep "%s" | cut -f 3 -d " "'%jobid_regex,shell=True,stdout=PIPE)
    stdout, stderr = qstat_p.communicate()
    job_list = stdout.split("\n")
    job_list = [x for x in job_list if len(x)>0]
    return(job_list)


def kill_jobs(job_list):
    print job_list
    for j in job_list:
        command = "qdel " + j
        check_call(command,shell=True)

def is_job_done(jobid) :

    done = False

    # have to handle array jobs differently than standalone
    array_match = re.search(array_job_match,jobid)
    if array_match is not None :
        idnum, rest = array_match.groups()
        jobid_regex = array_job_regex%idnum
        qstat_p = Popen('qstat -t | grep "%s" | cut -f 4 -d " "'%jobid_regex,shell=True,stdout=PIPE)
        stdout, stderr = qstat_p.communicate()
        done = len(stdout) == 0

    else :
        # -j is only for SGE
        qstat_p = Popen('qstat -j %s'%jobid,shell=True,stdout=PIPE,stderr=PIPE)
        qstat_p.wait()
        if qstat_p.returncode == 0 :
            pass
        # assume any != 0 return code means job is done
        else :
            done = True

    return done

def wait_for_jobs(job_list):
    sys.stderr.write('Waiting for jobs to complete\n')
    if len(job_list) > 0:
        jobs_done = [False]*len(job_list)
        try :
            while not all(jobs_done) :
                jobs_not_done = filter(lambda x: not x[1], enumerate(jobs_done))
                for i, jid in jobs_not_done :
                    jobs_done[i] = is_job_done(job_list[i])
                #sys.stderr.write('Jobs done: %d/%d\r'%(sum(jobs_done),len(jobs_done)))
                time.sleep(2)
                #sys.stderr.flush()
        except KeyboardInterrupt :
            sys.stderr.write('\n')
            resp = raw_input('Caught keyboard interrupt, kill all jobs? [y/N] ')
            if resp.lower() == 'y' :
                kill_jobs(job_list)
    else:
        print ("Provided job list is empty")
    sys.stderr.write('done\n')
