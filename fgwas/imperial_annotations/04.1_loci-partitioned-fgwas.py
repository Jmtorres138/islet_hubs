#!/usr/bin/python -O
# Jason Matthew Torres
'''
Run fGWAS across annotations
Usage: python 04.1_loci-partitioned-fgwas.py
'''
# libraries
import sys,os,gzip
import subprocess as sp
import operator
import time
from select import select
from math import ceil
from math import floor
import moniter_rescomp_jobs

# globals
fgwas = "LD_LIBRARY_PATH=/apps/well/gsl/2.2.1-gcc4.9.3/lib /users/mccarthy/jmtorres/software/fgwas-0.3.6/bin/fgwas"
eur_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/"
home_dir = "/well/mccarthy/users/jason/projects/islet_hubs/fgwas/imperial_annotations/"
in_dir=home_dir+"fgwas_input/"
out_dir = home_dir + "fgwas_output/"
#input_file=in_dir+"ukbb_diamante-euro.fgwas.gz"
# Optional: Run 01.1 script and use this for abbreviated annotations: in_dir+"diagram_hrc.renamed.fgwas.gz"
input_file=in_dir+"ukbb_diamante-euro.fgwas.gz" # explicitly removing transcript and intron from analysis
job_dir=home_dir+"jobs/"
log_dir=home_dir+"logs/"
if os.path.isdir(out_dir)==False:
    os.mkdir(out_dir)
if os.path.isdir(job_dir)==False:
    os.mkdir(job_dir)
if os.path.isdir(log_dir)==False:
    os.mkdir(log_dir)
start_index = 9 #0-based index of column in fgwas input file where annotations start
job_prefix = "imp_"

bed_file = eur_dir + "loci-partition.bed"
best_param_file = out_dir + "best-joint-model.params"


def get_annot_list(best_param_file):
    out_list = []
    fin = open(best_param_file,'r')
    fin.readline()
    fin.readline()
    for line in fin:
        l = line.strip().split()
        annot = l[0]
        annot = annot.split("_ln")[0].split("_0_5000")[0]
        out_list.append(annot)
    fin.close()
    return(out_list)

def fgwas_run(annot_list):
    job_file = job_dir+"fgwas_run_loci-partition.sh"
    fout=open(job_file,'w')
    if "distance_tss" in annot_list:
        a_list = list(annot_list)
        a_list.remove("distance_tss")
        command_list = [fgwas, "-i", input_file, "-cc",  "-bed", bed_file,
                        "-dists", "distance_tss:"+home_dir+"dist_model",
                        "-w", "+".join(a_list), "-print","-o", out_dir+"fgwas_run_loci-partition"]
    else:
        command_list = [fgwas, "-i", input_file, "-cc", "-bed", bed_file,
                        "-w", "+".join(annot_list), "-print","-o", out_dir+"fgwas_run_loci-partition"]
    command = " ".join(command_list)
    script='''
#$ -N %s%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
echo "start time" `date`
%s
echo "end time" `date`
    ''' % (job_prefix,"loci-partition", log_dir,"loci-partition",
    log_dir,"loci-partition", command)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call)


def main():
    annot_list = get_annot_list(best_param_file)
    fgwas_run(annot_list)

if (__name__=="__main__"): main()
