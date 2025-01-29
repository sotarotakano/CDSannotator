import sys
import os
import glob
import argparse
import collections
import numpy as np

from os.path import join
from multiprocessing import Pool
import subprocess
import re
import pathlib
import glob
import pandas as pd    
import configparser
import errno

def get_GTDB_taxonomy(arguments):
    originaldata=arguments[0]
    gtdb_genbankfile=arguments[1]
    diamonddata=arguments[2]
    gtdb_taxonfile=arguments[3]
    outputfile=arguments[4]
    locustag = originaldata["locustag"]
    diamondhit = [x for x in diamonddata if locustag in x]

    if len(diamondhit) > 0:
        max_pident = MAX_PIDENT # minimum threshold for besthit 
        min_length = MIN_LENGTH # minimum threshold for alignment length
        min_evalue = 1e-10 # maximum threshold for e-value
        diamondhit = [x for x in diamondhit if float(x.split("\t")[2]) >= max_pident]
        diamondhit = [x for x in diamondhit if float(x.split("\t")[3]) >= min_length]

        diamond_best = ""
        for hit in diamondhit:
            if float(hit.split("\t")[10]) < min_evalue:
                diamond_best = hit
                min_evalue = float(hit.split("\t")[10])

        if len(diamond_best) > 0:     
            genbank_id = diamond_best.split("\t")[1]
            # accession number of assembly file that becomes a query for GTDB taxonomy searching
            returns = subprocess.run(["rg","".join(genbank_id.split("_")[:-1]), gtdb_genbankfile],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            genbank_hits = returns.stdout.decode("utf-8").split("\n")[:-1]
            acc_no_list = [x.split('\t')[0] for x in genbank_hits]
            
            if len(acc_no_list) == 0:
                print("".join(genbank_id.split("_")[:-1]))
                print("this genbank id are not found in GTDB-TK list. something wrong...")
                return

            GTDB_taxonlist = []
            for acc_no in acc_no_list: 
                returns = subprocess.run(["rg",acc_no, gtdb_taxonfile],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                taxon = returns.stdout.decode("utf-8").replace("\n","").split("\t")[-1]
                if len(taxon) > 0:
                    GTDB_taxonlist.append(taxon)

            if len(GTDB_taxonlist) > 1:
                print("multiple candidates of GTDB taxonomy found... for %s"%locustag)
                GTDB_tax = "Unknown"
                GTDB_accno = "Unknown"
            else:
                GTDB_tax = GTDB_taxonlist[0]
                GTDB_accno = acc_no_list[0]

        else:
            GTDB_tax = "Unknown"
            GTDB_accno = "Unknown"
    else:
        GTDB_tax = "Unknown"
        GTDB_accno = "Unknown"
    
    GTDB_summarydata = [str(x) for x in originaldata] + [GTDB_tax, GTDB_accno]

    with open(outputfile, mode='a') as nf:
        nf.write('\t'.join(GTDB_summarydata) + '\n')



def add_GTDB_taxonomy(inputfile,gtdb_taxonfile,gtdb_genbankfile,diamondfile,outputfile,overwrite):
    
    print("Loading %s..."%inputfile)
    
    with open(diamondfile,mode="r") as f:
        lines = f.readlines()
        diamonddata= [x for x in lines]
    
    '''Create depth annotate file'''
    # load depth summary file
    summary = pd.read_csv(inputfile, sep='\t')
    nrows_summary = summary.shape[0]

    if os.path.exists(outputfile) and not overwrite:
        existdata = pd.read_csv(outputfile, sep='\t')
        nrows_output = existdata.shape[0]
        
        if nrows_summary == nrows_output:
            print("output files exist. skip this process.")
            return
        else:
            print("annotation process was not completed for %s, running the process again..."%inputfile)
    
    # get GTDB taxon info
    columns_name = '\t'.join([x for x in summary.columns.values] + ["GTDB_taxonomy","GTDB_acc_no"]) + '\n'

    with open(outputfile, mode='wt') as nf:
        nf.write(columns_name)
    
    with Pool(processes=Nthreads) as pool:
        pool.map(get_GTDB_taxonomy,[(summary.loc[x,:],gtdb_genbankfile,diamonddata,gtdb_taxonfile,outputfile) for x in summary.index])



def create_summary_annotate_GTDB(inputfile,outputfile,overwrite):
    '''Create annotation summary file'''

    def return_majority(item):
        if len(item) < 1:
            return ""
        else:
            c = collections.Counter(item)
            tophit = sorted(c.items(), key=lambda x:x[1], reverse=True)[0][0]
            return tophit

    #if os.path.exists(outputfile) and not overwrite:
    #    print("output files exist. skip this process.")
    #    return

    columns_name = '\t'.join(['taxon','primary_id','strain','isolate','organism','GTDB_accno','reads','length']) + '\n'

    with open(outputfile, mode='wt') as nf:
        nf.write(columns_name)

    # load depth annotation file
    summary = pd.read_csv(inputfile, sep='\t')
    summary = summary.fillna("")
    # with open(inputfile, mode='r') as f:
    #     lines = f.readlines()
    #     dataline = [x for x in lines if not "locustag" in x]
    # idlist = {x.split('\t')[6].strip('\n') for x in dataline}
    
    idlist = {x for x in summary["primary_id"]}
    
    for primary_id in idlist:
        # matches = [x.strip('\n') for x in dataline if primary_id == x.split('\t')[6].strip('\n')]
        matches = summary.loc[summary["primary_id"] == primary_id,:]
        annotated = matches.loc[matches["besthit"] != "Unknown",:]
       
        #annotated = annotated.sort_values("length",ascending = False)
        if annotated.shape[0] > 0:
            taxon = return_majority([str(x) for x in annotated["besthit"] if len(x) > 0])
            strain = return_majority([str(x) for x in annotated["strain"] if type(x) == "str"])
            isolate = return_majority([str(x) for x in annotated["isolate"] if len(x) > 0])
            organism = return_majority([str(x) for x in annotated["organism"] if len(x) > 0])
            GTDB_acc_no = return_majority([str(x) for x in annotated["GTDB_acc_no"] if len(x) > 0])
        else: 
            best_matches = matches.iloc[0,:]
            taxon = best_matches["besthit"]
            strain = best_matches["strain"]
            isolate = best_matches["isolate"]
            organism = best_matches["organism"]
            GTDB_acc_no = best_matches["GTDB_acc_no"]

        # for m in matches:
        #     if m.split('\t')[4] != "Unknown":
        #         taxon = m.split('\t')[4]
        #         strain = m.split('\t')[7]
        #         isolate = m.split('\t')[8]
        #         organism = m.split('\t')[9]
        #         break        
        # if len(taxon) == 0:
        #     taxon = matches[0].split('\t')[4]
        #     strain = matches[0].split('\t')[7]
        #     isolate = matches[0].split('\t')[8]
        #     organism = matches[0].split('\t')[9]
        reads = sum(matches["reads"])
        length = sum(matches["length"])
        
        #reads = sum([int(x.split('\t')[1]) for x in matches])
        #length = sum([int(x.split('\t')[2]) for x in matches])
        
        output = [taxon,str(primary_id),str(strain),str(isolate),str(organism),str(GTDB_acc_no),str(reads),str(length)]
       
        with open(outputfile, mode='a') as nf:
            nf.write('\t'.join(output) + '\n')


def create_summary_GTDB(inputfile,outputfile,overwrite):
    '''Create annotation summary file based on GTDB info'''
    def return_majority(item):
        if len(item) > 0:
            c = collections.Counter(item)
            tophit = sorted(c.items(), key=lambda x:x[1], reverse=True)[0][0]
            return tophit
        else:
            return ""

    print("Loading %s..."%inputfile)

    if os.path.exists(outputfile) and not overwrite:
        print("%s files exist. skip this process."%outputfile)
        return

    columns_name = '\t'.join(['taxon','GTDB_accno','reads','length']) + '\n'

    with open(outputfile, mode='wt') as nf:
        nf.write(columns_name)

    # load depth annotation file
    summary = pd.read_csv(inputfile, sep='\t')
    summary = summary.fillna("")

    
    if "GTDB_acc_no" not in summary.columns.values:
        print("no GTDB acc no. found. stop the process...")
        return
       
    accnolist = {x for x in summary["GTDB_acc_no"]}
    
    for accno in accnolist:
        matches = summary.loc[summary["GTDB_acc_no"] == accno,:]
       
        if matches.shape[0] > 1:
            taxon = return_majority([str(x) for x in matches["GTDB_taxonomy"] if len(str(x)) > 0])
            GTDB_acc_no = return_majority([str(x) for x in matches["GTDB_acc_no"] if len(str(x)) > 0])
        else: 
            best_matches = matches.iloc[0,:]
            taxon = best_matches["GTDB_taxonomy"]
            GTDB_acc_no = best_matches["GTDB_acc_no"]

        reads = sum(matches["reads"])
        length = sum(matches["length"])
        
        output = [taxon,str(GTDB_acc_no),str(reads),str(length)]
       
        with open(outputfile, mode='a') as nf:
            nf.write('\t'.join(output) + '\n')


if __name__ == "__main__":
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--inputdir", type=str, default="None", help="input_dir")
    p.add_argument("-d","--diamondfile", type=str, default="None", help="diamondfile")
    p.add_argument("-ow","--over_write", action='store_true')
    p.add_argument("-mp","--min_pident", type=int,  default=80, help="minimum_pident")
    p.add_argument("-ml","--min_length", type=int,  default=100, help="minimum_length")
    args = p.parse_args()

    MAX_PIDENT=args.min_pident
    MIN_LENGTH=args.min_length
    
    # Load config info
    config_ini = configparser.ConfigParser()
    configfile=join(str(pathlib.Path(__file__).parent.parent),'python_config.ini')
    config_ini.read(configfile, encoding='utf-8')

    if not os.path.exists(configfile):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), configfile)
    
    GTDB_taxon = config_ini['DEFAULT']['GTDB_TAXON_FILE']
    GTDB_genbank = config_ini['DEFAULT']['GTDB_GENBANK_FILE']
    Nthreads = int(config_ini['DEFAULT']['THREADS'])

    filelist = glob.glob(join(args.inputdir,"bowtie2_results","*","*_depth_annotate.txt"))

    for f in filelist:
        print(f)
        add_GTDB_taxonomy(f,GTDB_taxon,GTDB_genbank,args.diamondfile,join(pathlib.Path(f).parent,f.split("/")[-1].split("_")[0]+"_depth_annotate_GTDB.txt"),args.over_write)
    
    filelist = glob.glob(join(args.inputdir,"bowtie2_results","*","*_depth_annotate_GTDB.txt"))

    for i in filelist:
        print(i)
        create_summary_annotate_GTDB(i,join(pathlib.Path(i).parent,i.split("/")[-1].split("_")[0]+"_annotate_summary.txt"),args.over_write)
        create_summary_GTDB(i,join(pathlib.Path(i).parent,i.split("/")[-1].split("_")[0]+"_annotate_summary_GTDB.txt"),args.over_write)
