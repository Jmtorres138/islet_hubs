#!/usr/bin/python -O
# Jason Matthew Torres
'''
Usage: python 06.0_subset_to_seg.py loc_id
'''
# libraries
import sys,os,gzip


# globals
server_dir = "/well/mccarthy/users/jason/"
work_dir = server_dir + "projects/islet_hubs/fgwas/oxford_and_imperial/"

loc_id = sys.argv[1]
pre = work_dir + "conditional/fgwas_output_files/" + loc_id + "/" + "fgwas_run_loci-partition"
sig_block_file = work_dir + "conditional/fgwas_output_files/" + loc_id + "/" + "fgwas_blk.txt"
output_file = work_dir+"conditional/fgwas_output_files/" + loc_id + "/" +"loci_block_snps.bfs.txt.gz"

def build_dic():
    fin = open(sig_block_file,'r')
    header = fin.readline()
    dic = {}
    for line in fin:
        l = line.strip().split()
        seg, chrom, st, sp, ppa = l[0],l[1],l[2],l[3],l[4]
        try:
            dic[chrom].append([int(st),int(sp),seg,ppa])
        except:
            dic[chrom] = [[int(st),int(sp),seg,ppa]]
    fin.close()
    return(dic)

def subset_to_seg():
    print("Building dictionary...")
    ref_dic =  build_dic()
    fin = gzip.open(pre+".bfs.gz",'rb')
    head_list = fin.readline().strip().split()
    head_list = ["SEGNUMBER"] + head_list
    fout = gzip.open(output_file,'wb')
    print("Writing output file...")
    fout.write(" ".join(head_list)+"\n")
    count=0
    for line in fin:
        count+=1
        #sys.stdout.write("\r%d" % count)
        #sys.stdout.flush()
        l = line.strip().split()
        snpid, chrom, pos = l[0],l[1],int(l[2])
        try:
            for li in ref_dic[chrom]:
                start, stop, seg = li[0],li[1],li[2]
                if pos >= start and pos <= stop:
                    print ("\n" + snpid + "\n")
                    write_list = [seg] + l
                    fout.write(" ".join(write_list)+"\n")
        except:
            pass
    fout.close()
    fin.close()



def main():
    subset_to_seg()

if (__name__=="__main__"): main()
