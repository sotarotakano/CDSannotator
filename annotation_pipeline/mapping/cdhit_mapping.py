import sys
import os
import gc
import argparse 
from os.path import join

from multiprocessing import Pool
from contextlib import closing
from concurrent.futures import ThreadPoolExecutor
import pathlib
import glob
from multiprocessing import Pool
import pandas as pd
import numpy as np  

def search_rep_CDS(args):
    locustag = args[0]
    depthdata = args[1]
    clstrdata = args[2]
    max_clstr_num = args[3]
    outputfile = args[4]

    target_sequence = [x for x in depthdata if locustag in x]
    if len(target_sequence) > 0:
        target_sequence = target_sequence[0]
        length = target_sequence.split("\t")[1]
        reads = target_sequence.split("\t")[2]
        #length = len(target_sequence)
        #reads = sum([int(x.strip("\n").split('\t')[2]) for x in target_sequence])
    else:
        sys.stderr.write("The locus tag of %s was not found ... " + '\n'%locustag)
        return
    
   
    #searching corresponding cluster id
    row_num = [clstrdata.index(x) for x in clstrdata if locustag+"..." in x]
    if len(row_num) > 1:
        print(locustag)
        sys.stderr.write("The locus tag was found in multiple clusters. something wrong... " + '\n')
        return
    elif len(row_num) == 0:
        print(locustag)
        sys.stderr.write("The locus tag was not found in any clusters. something wrong... " + '\n')
        return
    else:
        row_num = row_num[0]
    
    if "*" in clstrdata[row_num]:
        rep_id = clstrdata[row_num].split(">")[1].split("...")[0]
    else:
        search_start = max(0,row_num-max_clstr_num)
        clstr_num = [clstrdata.index(x) for x in clstrdata[search_start:row_num] if ">Cluster" in x][-1]
        search_end = min(clstr_num+max_clstr_num+1,len(clstrdata))
        try:
            rep_cds = [x for x in clstrdata[clstr_num:search_end] if "*" in x][0]
        except:
            print([locustag,search_start,search_end,clstr_num])
        rep_id = rep_cds.split(">")[1].split("...")[0]
    
    summarydata = [rep_id,str(length),str(reads)]

    with open(outputfile, mode='a') as nf:
        nf.write('\t'.join(summarydata) + '\n')

    gc.collect()



def cdhit_mapping(inputfile,clstrfile,outputfile,overwrite):
    
    print("[CDHIT-mapping] Loading %s..."%inputfile)
    
    # load bastadata, diamondata
    with open(inputfile,mode="r") as f:
        lines = f.readlines()
        depthdata= [x for x in lines if "*" not in x]
    locustaglist = sorted({x.split('\t')[0] for x in depthdata})


    with open(clstrfile,mode="r") as f:
        lines = f.readlines()
        clstrdata= [x for x in lines]
    
    
    '''Create depth summary file'''
    if os.path.exists(outputfile) and not overwrite:
        with open(outputfile, mode='r') as f:
            existdata = [x for x in f.readlines() if not "locustag" in x]
            nrows_summary = len(existdata)
        if nrows_summary == len(locustaglist):
            print("[CDHIT-mapping] output files exist. %s \n Skip this process."%outputfile)
            return
        else:
            print("[CDHIT-mapping] summary file was not complete for %s, running the process again..."%inputfile)
            

    columns_name = '\t'.join(['locustag','length','reads']) + '\n'

    with open(outputfile, mode='wt') as nf:
        nf.write(columns_name)
    
    max_clstr_num = max({int(x.split("\t")[0]) for x in clstrdata if ">Cluster" not in x}) + 1

    with Pool(processes=18) as pool:
        pool.map(search_rep_CDS,[(x,depthdata,clstrdata,max_clstr_num,outputfile) for x in locustaglist])
    



if __name__ == "__main__":
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--inputdir", type=str, default="None", help="input_dir")
    p.add_argument("-c","--clstr_dat", type=str, default="None", help="cluster file of cd-hit")
    p.add_argument("-ow","--over_write", action='store_true')
    args = p.parse_args()

    if args.over_write:
        ow = args.over_write 
    else:
        ow = False
    
    if args.inputdir== "None":
        print("[ERROR] input directory must be set by -i option.")
        sys.exit()
    else:
        inputdir = args.inputdir

    if args.clstr_dat == "None":
        print("[ERROR] the location of clstr file in cd-hit must be set by -c option.")
        sys.exit()
    else:
        clstrfile = args.clstr_dat

    #filelist = glob.glob(join(inputdir,"bowtie2_results","*","*_depth.txt"))
    filelist = glob.glob(join(inputdir,"bowtie2_results","*","*_bowtie2stats.txt"))

    #args = [(x,clstrfile,join(pathlib.Path(x).parent,x.split("/")[-1].split("_")[0]+"_depth_summary.txt"),ow) for x in filelist]
    #with Pool(processes=16) as pool:
    #   pool.map(cdhit_mapping_multi,args)
    
    for f in filelist:
        cdhit_mapping(f,clstrfile,join(pathlib.Path(f).parent,f.split("/")[-1].split("_")[0]+"_depth_summary.txt"),ow)