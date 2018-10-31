
#!/usr/bin/python -O
# Jason Matthew Torres
'''
Input a gwas bed file and an annotation bed file
Usage: python JTbuild_fgwas_bedfile.py
'''
# libraries
import sys,os,gzip
import subprocess as sp

# globals
cur_dir = "/well/mccarthy/users/jason/projects/islet_hubs/fgwas/oxford_and_imperial/"
input_dir = cur_dir + "fgwas_input/"

infile = input_dir+"ukbb_diamante-euro.fgwas.gz"
outfile = input_dir+"ukbb_diamante-euro.renamed.fgwas.gz"

start_index = 9 #0-based index of column in fgwas input file where annotations start

# functions

def rename_annots(input_file,output_file):
    fin = gzip.open(input_file,'rb')
    fout = gzip.open(output_file,'wb')
    head_list = fin.readline().strip().split()
    print "Writing key file..."
    keyfile = open(input_dir+"annotation_key-file.txt",'w')
    keyfile.write("\t".join(["Key","Annotation"])+"\n")
    num = 0
    annot_list = []
    for a in head_list[start_index:]:
        num += 1
        if a=="distance_tss":
            key = a
        else:
            key = "A"+str(num)
        annot_list.append(key)
        keyfile.write("\t".join([key,a])+"\n")
    keyfile.close()
    new_head_list = head_list[0:(start_index)] + annot_list
    print "Writing renamed fgwas input file...."
    fout.write("\t".join(new_head_list)+"\n")
    count=0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        line = line.strip()
        fout.write(line+"\n")
    fin.close()
    fout.close()
    print "\nProcess Complete"


def main():
    rename_annots(infile,outfile)


if (__name__=="__main__"): main()
