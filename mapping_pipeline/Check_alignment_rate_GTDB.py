from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import numpy as np
import time
import subprocess
import glob
import pandas as pd
import os
import sys
from os.path import join
import re
import argparse
import configparser
import errno
import pathlib


def get_assembly_GTDB(parent_dir, sample_id, gtdb_no, organism, assembly_dir, maxret):
    
    print("%s:%s"%(gtdb_no,organism))

    # get assembly based on GTDB id
    try:
        GC_file = re.search("GC.+",gtdb_no).group(0)
        handle = Entrez.esearch(db="assembly", retmax=5,
                                term=GC_file.strip('\n'), report="full")
        record = Entrez.read(handle)
        id_gtdb = record['IdList']
        print("found gtdb ref seq, %s"%gtdb_no)
    except Exception as exc:
        id_gtdb = []
        sys.stderr.write(str(exc) + '\n')

    ASM_data = pd.DataFrame(0, index=id_gtdb,
                            columns=["acc_no", "url", "alignrate", "coverrate"])
    
    # get esummary
    for id_str in id_gtdb:
        try:
            esummary_handle = Entrez.esummary(
                db="assembly", id=id_str, report="full")
            esummary_record = Entrez.read(esummary_handle, validate=False)
        except Exception as exc:
            print("Error in %s"%id_str)
            sys.stderr.write(str(exc) + '\n')
            continue

        # get fna file
        try:
            acc_no = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
        except Exception as exc:
            print("Error in %s"%id_str)
            sys.stderr.write(str(exc) + '\n')
            continue

        if len(acc_no) == 0:
            print("No accession number is found. Skip this assembly file...")
            continue

        #organism_name = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        #organism = str(organism)
        #if len(organism) > 0:
        #    jw_distance = Levenshtein.jaro_winkler(organism_name, organism)
        #else:
        #    jw_distance = None

        if "GCF" in acc_no:
            url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        elif "GCA" in acc_no:
            url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        else:
            print("irregular accession number:%s"%acc_no)


        filename = url.split('/')[-1] + '_genomic.fna.gz'

        ASM_data.loc[id_str, "acc_no"] = acc_no
        ASM_data.loc[id_str, "url"] = url
        ASM_data.loc[id_str,"filename"] = join(assembly_dir,filename)
        ASM_data.loc[id_str,"status"] = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatusSort']
    
    print("%i assembly files are obtained"%len(ASM_data))
    if len(ASM_data) < 1:
        print("No assembly file is obtained...Skip this ID...")
        return

    ASM_data = ASM_data.sort_values("status")
    ASM_data = ASM_data.iloc[:maxret,:]
    
    bowtie2_data = join(parent_dir,sample_id + "_depth.txt")
    
    for idx in ASM_data.index:
        url = ASM_data.loc[idx,"url"].replace("ftp://","https://")
        assembly_file = ASM_data.loc[idx,"filename"]
        filename = os.path.basename(assembly_file)
        print("[ENTREZ] Getting assembly file... %s"%filename)
        # Downloading reference genomic file
        if os.path.exists(assembly_file.replace(".gz","")):
            #print("Found %s in the directory." % assembly_file)
            pass
        else:
            try:
                print("Fetching %s" % join(url, filename))
                wget_cmd = "wget " + "-P " + assembly_dir + " " + \
                    join(url, filename) + "; gzip -d " + \
                    assembly_file
                returns = subprocess.run(
                    ["wget", "-P", assembly_dir, join(url, filename)],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(returns.stderr.decode("utf-8"))
                if re.search("Not Found",returns.stderr.decode("utf-8")):
                    continue
                returns = subprocess.run(
                    ["gzip", "-d", join(assembly_dir, filename)],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                #returns = subprocess.check_output(wget_cmd, shell = True)
            except Exception as exc:
                sys.stderr.write(str(exc) + '\n')
                continue

        # Check file size (if it's larger than 100MB, this is too large for bacteria.)
        if not os.path.exists(assembly_file.replace(".gz","")):
            print("%s download failed..."%assembly_file)
            continue
        else:
            if os.path.getsize(assembly_file.replace(".gz","")) > 1E+08:
                os.remove(assembly_file.replace(".gz",""))
                print("%s is too large, the mapping is canceled..."%assembly_file.replace(".gz",""))
                continue

        # Performing Bowtie2 (build)
        #print("[BOWTIE2] building bowtie2 reference %s"%assembly_file)

        returns = subprocess.run(
            [bt_build_script, "-a", assembly_file.replace(".gz","") ,"--filename"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("[BOWTIE2] Reference file generated " +returns.stdout.decode("utf-8"))
        bowtie2ref = [x for x in returns.stdout.decode("utf-8").split("\n") if "home" in x][0]
        #print(bowtie2ref)

        # Performing Bowtie2 (alignment)
        returns = subprocess.run([bowtie2_script,"-d",parent_dir,"-i",sample_id,"-r",bowtie2ref,
                                  "-o",parent_dir],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(returns.stderr.decode("utf-8"))
        result = [x for x in returns.stderr.decode("utf-8").split("\n") if "overall alignment rate" in x]
        if len(result) < 1: 
            alignment_rate = np.nan
        else:
            alignment_rate = float(result[0].split("%")[0])

        print("%s|%s|%s:%f" % (acc_no, id_str, organism, alignment_rate))
        
        ASM_data.loc[idx,"alignrate"] = alignment_rate

        # Check covered rate
        if os.path.exists(bowtie2_data):
            with open(bowtie2_data, mode='r') as f:
                lines = f.readlines()
                mapdata = np.array([int(x.split('\t')[2]) for x in lines])
        
            mapped = len(mapdata[mapdata > 0]) 
            total = len(mapdata) 
        
            if total > 0:
                ASM_data.loc[idx,"coverrate"] = mapped/total
                print("Cover rate:%f"%(mapped/total))
            else:
                ASM_data.loc[idx,"coverrate"] = 0
        else:
             ASM_data.loc[idx,"coverrate"] = 0


        if idx in id_gtdb:
            ASM_data.loc[idx, "GTDB"] = 1
            print("%s in GTDB taxonomy"%idx)
    
    ASM_data = ASM_data.sort_values("alignrate", ascending=False)
    
    if os.path.exists(bowtie2_data):
        os.remove(bowtie2_data)
    
    return ASM_data



def check_alignment_rate(summary_file,parent_dir,overwrite,max_ret,Ntophits=10,assembly_dir = "/home/soutacano/Ref_genomes"):
    start_time = time.time()
    summary = pd.read_csv(summary_file, sep='\t')
    sample_id = summary_file.split('/')[-1].split('_')[0]

    if os.path.exists(join(parent_dir,"bowtie2_results",sample_id,sample_id + "_besthits_GTDB.txt")) and not overwrite:
        print(join(parent_dir,"bowtie2_results",sample_id,sample_id + "_besthits_GTDB.txt") + "exists. Skip this process.")
        return

    total_reads = sum(summary["reads"])
    total_length = sum(summary["length"])

    summary["frequency"] = 0
    for i in summary.index:
        summary.loc[i,"frequency"] = summary.loc[i,"length"]/total_length


    abundant = summary.sort_values("frequency", ascending = False).head(Ntophits)
    print("Candidates for MFTs.")
    print(abundant.loc[:,["taxon","frequency","GTDB_accno"]].reset_index().drop(columns='index'))
    for i in abundant.index:
        candidate = get_assembly_GTDB(parent_dir, sample_id, abundant.loc[i,"GTDB_accno"],abundant.loc[i,"taxon"].split(";")[-1], assembly_dir,max_ret)

        if candidate is not None:
            candidate['taxon'] = abundant.loc[i,"taxon"]      
            candidate['reads'] = abundant.loc[i,"reads"]              
            candidate['length'] = abundant.loc[i,"length"]
            candidate['frequency'] = abundant.loc[i,"frequency"]
        else:
            continue

        if "besthits" in locals():
            besthits = pd.concat([besthits,candidate.T], axis = 1)
        else:
            besthits = candidate.T


    print("%s finished, took %3f min"%(summary_file,(time.time()-start_time)/60))
    start_time = time.time()
    if "besthits" in locals():     
        besthits = besthits.T
        besthits.to_csv(join(parent_dir,"bowtie2_results",sample_id,sample_id + "_besthits_GTDB.txt"), sep = "\t", index = False)
        return besthits
    else:
        return

if __name__ == "__main__":
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--inputdir", type=str, default="None", help="input_dir")
    p.add_argument("-n","--ntophits",type=int,default=10, help="maximum numbers of candidates")
    p.add_argument("-r","--max_ret", type=int, default=10, help="maximum numbers of retrieves")
    p.add_argument("-ow","--over_write", action='store_true')
    args = p.parse_args()

    # Load config info
    config_ini = configparser.ConfigParser()
    configfile=join(str(pathlib.Path(__file__).parent),'python_config.ini')
    config_ini.read(configfile, encoding='utf-8')

    if not os.path.exists(configfile):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), configfile)


    Nthreads = int(config_ini['DEFAULT']['THREADS'])
    # E-mail address
    # You first visit NCBI registration page and create NCBI account
    # then the fill in the registered e-mail address followings
    Entrez.email = config_ini['DEFAULT']['ENTREZ_EMAIL']

    # NCBI API key: you can get in the account page
    Entrez.api_key =  config_ini['DEFAULT']['ENTREZ_API']
    bowtie2_script = config_ini['DEFAULT']['BOWTIE2_SCRIPT']
    bt_build_script = config_ini['DEFAULT']['BOWTIE2_BUILD_SCRIPT']

    bowtie_data = glob.glob(join(args.inputdir,"bowtie2_results/*/*_annotate_summary_GTDB.txt"))
    print("%i GTDB summary data are found."%len(bowtie_data))

    for filepath in bowtie_data:
        print(filepath)
        besthits = check_alignment_rate(filepath,args.inputdir,args.over_write,args.max_ret,Ntophits=args.ntophits)