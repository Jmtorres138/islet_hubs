#!/usr/bin/python -O
# Jason Matthew Torres
'''
Run fGWAS across annotations
Usage: python 04.1_conditional-fgwas-wrapper.py
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
python = "/apps/well/python/2.7.11/bin/python"
rscript = "/apps/well/R/3.3.1/bin/Rscript" #"/apps/well/R/3.3.1/bin/Rscript"
fgwas = "LD_LIBRARY_PATH=/apps/well/gsl/2.2.1-gcc4.9.3/lib /users/mccarthy/jmtorres/software/fgwas-0.3.6/bin/fgwas"
eur_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/"
home_dir = "/well/mccarthy/users/jason/projects/islet_hubs/fgwas/oxford_and_imperial/"
part_dir = eur_dir + "partition_bedfiles/"
in_dir=home_dir+"conditional/fgwas_input_files/"
out_dir = home_dir + "conditional/fgwas_output_files/"
job_dir=home_dir+"jobs/"
log_dir=home_dir+"logs/"

manual_list = ["210_1","210_2","86_1","87_1", "132_1", "132_2", "132_3", "132_4", "132_5",
                "133_1","133_2","133_3","133_4","133_5","133_6","133_7","133_8","133_9","133_10"]
missing_list = ["20_1","86_2","87_2","163_2","164_2"]

gwas_bed_dir = "/well/got2d/jason/reference/gwas/diamante-ukbb_hrc/conditioned/gwas_bedfiles/"
best_annot_bed = in_dir+"best_anno_input.bed"

best_param_file = home_dir + "fgwas_output/best-joint-model.params"

loc_ref_file = "/well/got2d/jason/reference/gwas/diamante-ukbb_hrc/conditioned/list.for.credible.sets.ALL.txt"

key_file = home_dir + "fgwas_input/annotation_key-file.txt"

# functions

def run_fgwas_input_job(gwas_bed):
    ref_name = gwas_bed.split("ukbb_diamante-euro.")[1].split(".bed")[0]
    out_file = in_dir + "ukbb_diamante-euro." + ref_name + ".fgwas"
    job_file = job_dir+"job."+ref_name+".sh"
    fout=open(job_file,'w')
    command_list = [python,home_dir+"05.0_build_fgwas_input.py",gwas_bed,out_file]
    command = " ".join(command_list)
    script='''
#$ -N cond2_%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
echo "start time" `date`
%s
echo "end time" `date`
        ''' % (ref_name, log_dir,ref_name,log_dir,ref_name, command)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call)

def build_all_inputs():
    gwas_bed_files = os.listdir(gwas_bed_dir)
    for f in gwas_bed_files:
        print f
        run_fgwas_input_job(gwas_bed_dir+f)

def build_key_dic():
    fin = open(key_file,'r')
    dic = {}
    fin.readline() # header
    for line in fin:
        l = line.strip().split()
        a,k = l[0],l[1]
        dic[a]=k
    fin.close()
    return(dic)

def get_annot_list(best_param_file,key=True):
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
    if key==True:
        dic = build_key_dic()
        for i in range(0,len(out_list)):
            mykey = out_list[i]
            ann = dic[mykey]
            print [mykey, ann]
            out_list[i] = ann
    return(out_list)

def run_loc_job(loc_id):
    # Note: loc_id can be locus id (e.g. "188_2") or condition ref name (e.g. "cond2")
    if loc_id in manual_list:
        fgwas_input_file = in_dir + "ukbb_diamante-euro." + loc_id + ".fgwas.gz"
    elif loc_id in missing_list:
        cond = loc_id.split("_")[1]
        fgwas_input_file = in_dir + "ukbb_diamante-euro.cond" + str(cond) + ".fgwas.gz"
    elif "_" in loc_id:
        cond = loc_id.split("_")[1]
        fgwas_input_file = in_dir + "ukbb_diamante-euro.cond" + str(cond) + ".fgwas.gz"
    else:
        cond = loc_id.split("cond")[1]
        fgwas_input_file = in_dir + "ukbb_diamante-euro.cond" + str(cond) + ".fgwas.gz"
    bed_part_file = part_dir + "loci-partition-"+loc_id + ".bed"

    loc_out_dir = out_dir + loc_id + "/"
    if os.path.isdir(loc_out_dir)==False:
        os.mkdir(loc_out_dir)

    job_file = job_dir+"job."+loc_id+".sh"
    fout=open(job_file,'w')

    command_list1 = [rscript,"--vanilla",home_dir+" 05.0_conditional-create-bed.R",loc_id]
    command1 = " ".join(command_list1)

    annot_list = get_annot_list(best_param_file)
    if "distance_tss" in annot_list:
        a_list = list(annot_list)
        a_list.remove("distance_tss")
        command_list2 = [fgwas, "-i", fgwas_input_file, "-cc",  "-bed", bed_part_file,
                        "-dists", "distance_tss:"+home_dir+"dist_model",
                        "-w", "+".join(a_list), "-print","-o", loc_out_dir+"fgwas_run_loci-partition"]
    else:
        command_list2 = [fgwas, "-i", fgwas_input_file, "-cc", "-bed", bed_part_file,
                        "-w", "+".join(annot_list), "-print","-o", loc_out_dir+"fgwas_run_loci-partition"]
    command2 = " ".join(command_list2)

    command_list3 = [rscript, "--vanilla",home_dir+"05.0_get-cond-block.R",loc_id]
    command3 = " ".join(command_list3)

    command_list4 = [python,home_dir+"05.0_subset_to_seg.py",loc_id]
    command4 = " ".join(command_list4)

    command_list5 = [rscript, "--vanilla",home_dir+"05.0_functional_credible_sets.R",loc_id]
    command5 = " ".join(command_list5)

    command_list6 = ["rm", loc_out_dir+"fgwas_run_loci-partition.bfs.gz"]
    command6 = " ".join(command_list6)

    script='''
#$ -N cond2_%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V
echo "start time" `date`
%s
%s
%s
%s
#%s
#%s
echo "end time" `date`
        ''' % (loc_id, log_dir,loc_id,log_dir,loc_id,
        command1, command2, command3, command4, command5, command6)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call)

def run_all_loci():
    loc_list = []
    for i in [1,2,3,4,5,6,7,8,9,10]:
        loc_list.append("cond"+str(i))
    loc_list = loc_list + manual_list
    loc_list = loc_list + missing_list
    #loc_list = manual_list # TEMPORARY REMOVE WHEN FINISHED
    for loc in loc_list:
        run_loc_job(loc)
    print loc_list

def run_cred_sets():
    loc_list = []
    for i in [1,2,3,4,5,6,7,8]:#,9,10]:
        loc_list.append("cond"+str(i))
    loc_list = loc_list + manual_list
    loc_list = loc_list + missing_list
    #loc_list = manual_list # TEMPORARY REMOVE WHEN FINISHED
    #loc_list = missing_list # TEMPORARY REMOVE WHEN FINISHED
    for loc in loc_list:
        print(loc)
        command_list5 = [rscript, "--vanilla",home_dir+"05.0_functional_credible_sets.R",loc]
        sp.check_call(command_list5)
    print loc_list

def main():

    #build_all_inputs()
    #run_all_loci()
    run_cred_sets()


if (__name__=="__main__"): main()
