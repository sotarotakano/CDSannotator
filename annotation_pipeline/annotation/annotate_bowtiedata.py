import sys
import os
import time
import argparse
import collections
from os.path import join

from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from contextlib import closing
import re
import pathlib
import glob
from multiprocessing import Pool
import pandas as pd
import configparser
import errno


def get_annotation(args):
    
    line = args[0]
    bastadata = args[1]
    diamonddata = args[2]
    outputfile = args[3]
    get_primary_id = args[4]

    line = line.strip("\n")
    
    reads = line.split("\t")[2]
    length = line.split("\t")[1]
    locustag = line.split("\t")[0]
    
    #searching basta hit
    bastahit = [x.strip("\n") for x in bastadata if locustag in x]
    if len(bastahit) > 0:
        bastahit = bastahit[0]
        LCA = bastahit.split("\t")[1]
        besthit = bastahit.split("\t")[2]
    elif len(bastahit) == 0:
        LCA = "Unknown"
        besthit = "Unknown"
    
    diamondhit = [x for x in diamonddata if locustag in x]
    
    if len(diamondhit) > 0:
        max_pident = 80 # minimum threshold for besthit 
        diamond_best = ""
        for hit in diamondhit:
            if float(hit.split("\t")[2]) > max_pident:
                diamond_best = hit
                max_pident = float(hit.split("\t")[2])
        if len(diamond_best) > 0:     
            genbank_id = diamond_best.split("\t")[1]
            if get_primary_id:
                [strain_id, isolate_id, organism_id] = get_entrez_info(genbank_id,"isolate")
            else:
                strain_id = ""
                isolate_id = ""
                organism_id = ""
        else:
            genbank_id = "Unknown"
            strain_id = ""
            isolate_id = ""
            organism_id = ""

    else:
        genbank_id = "Unknown"
        strain_id = ""
        isolate_id = ""
        organism_id = ""
    

    try:
        primary_id = [x for x in [strain_id, isolate_id, organism_id] if len(x) > 0][0]
    except:
        primary_id = "Unknown"
    
    annotation = [locustag,reads,length,LCA,besthit,genbank_id,primary_id,strain_id,isolate_id,organism_id]

    with open(outputfile, mode='a') as nf:
        nf.write('\t'.join(annotation) + '\n')



def annotate_to_depth(inputfile,bastafile,diamondfile,outputfile,overwrite,get_primary_id):
    
    print("Loading %s..."%inputfile)
    
    # load bastadata, diamondata
    with open(bastafile,mode="r") as f:
        lines = f.readlines()
        bastadata= [x for x in lines]

    with open(diamondfile,mode="r") as f:
        lines = f.readlines()
        diamonddata= [x for x in lines]
        
        
    '''Create depth annotate file'''
    # load depth summary file
    with open(inputfile, mode='r') as f:
        lines = f.readlines()
        dataline = [x for x in lines if not "locustag" in x]
        nrows_summary = len(dataline)
    

    if os.path.exists(outputfile) and not overwrite:
        with open(outputfile, mode='r') as f:
            existdata = [x for x in f.readlines() if not "locustag" in x]
            nrows_annotate = len(existdata)
        
        if nrows_summary == nrows_annotate:
            print("output files exist. skip this process.")
            return
        else:
            print("annotation process was not completed for %s, running the process again..."%inputfile)
            

    columns_name = '\t'.join(['locustag','reads','length','LCA','besthit','genbankid','primary_id','strain','isolate','organism']) + '\n'

    with open(outputfile, mode='wt') as nf:
        nf.write(columns_name)
    
    with Pool(processes=Nthreads) as pool:
        pool.map(get_annotation,[(x,bastadata,diamonddata,outputfile,get_primary_id) for x in dataline])
      
    


def create_summary_annotate(inputfile,outputfile,overwrite):
    '''Create annotation summary file'''
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
    
    if "GTDB_acc_no" not in summary.columns.values:
        summary["GTDB_acc_no"] = ""
       
    idlist = {x for x in summary["primary_id"]}
    
    for primary_id in idlist:
        # matches = [x.strip('\n') for x in dataline if primary_id == x.split('\t')[6].strip('\n')]
        matches = summary.loc[summary["primary_id"] == primary_id,:]
        annotated = matches.loc[matches["besthit"] != "Unknown",:]
       
        #annotated = annotated.sort_values("length",ascending = False)
        if annotated.shape[0] > 1:
            taxon = return_majority([str(x) for x in annotated["besthit"] if len(str(x)) > 0])
            strain = return_majority([str(x) for x in annotated["strain"] if len(str(x)) > 0])
            isolate = return_majority([str(x) for x in annotated["isolate"] if len(str(x)) > 0])
            organism = return_majority([str(x) for x in annotated["organism"] if len(str(x)) > 0])
            GTDB_acc_no = return_majority([str(x) for x in annotated["GTDB_acc_no"] if len(str(x)) > 0])
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



def get_entrez_info(genbankid,output="taxon"):
    # E-mail address
    # You first visit NCBI registration page and create NCBI account
    # then the fill in the registered e-mail address followings
    Entrez.email = ENTREZ_EMAIL

    # NCBI API key: you can get in the account page
    Entrez.api_key =ENTREZ_API
    
    error = 1
    for i in range(10):
        try:
            handle = Entrez.efetch(db = 'protein', id = genbankid, rettype = 'gb', retmode = 'text')
            error = 0
            break
        except Exception as exc:
            sys.stderr.write(str(exc) + '\n')
            time.sleep(0.5)

    if error:
        for i in range(10):
            try:
                handle = Entrez.efetch(db = 'nucleotide', id = genbankid, rettype = 'gb', retmode = 'text')
                error = 0
                break
            except Exception as exc:
                sys.stderr.write(str(exc) + '\n')
                print("%s is not found in ncbi database"%genbankid)
                time.sleep(0.5)
    if error:
        return ["", "", ""]
    
    if output == "taxon":
        record = SeqIO.read(handle, "genbank")
        organ = record.annotations['taxonomy']
        return organ
    elif output == "isolate":
        record = SeqIO.read(handle, "genbank")
        for feature in record.features:
            if feature.type == 'source':
                if 'strain' in feature.qualifiers:
                    strain = feature.qualifiers['strain']
                    strain_str = str(strain[0])
                else:
                    strain_str = ""

                if 'isolate' in feature.qualifiers:
                    isolate = feature.qualifiers['isolate']
                    isolate_str = str(isolate[0])
                else:
                    isolate_str = ""

                if 'organism' in feature.qualifiers:
                    organism = feature.qualifiers['organism']
                    organism_str = str(organism[0])
                else:
                    organism_str=""    

        return [strain_str, isolate_str, organism_str]
    else:
        print("\n# target must be set to 'taxon' or 'isolate'.")
        return ["", "", ""]



if __name__ == "__main__":
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--inputdir", type=str, default="None", help="input_dir")
    p.add_argument("-d","--diamond_dat", type=str, default="None", help="output of diamond")
    p.add_argument("-b","--basta_dat", type=str, default="None", help="output of basta")
    p.add_argument("-ow","--over_write", action='store_true')
    p.add_argument("-gp","--get_primary_id", action='store_true')

    args = p.parse_args()

    if args.over_write:
        ow = args.over_write 
    else:
        ow = False
    
    if args.get_primary_id:
        get_primary_id = args.get_primary_id 
    else:
        get_primary_id = False

    if args.inputdir== "None":
        print("[ERROR] input directory must be set by -i option.")
        sys.exit()
    else:
        inputdir = args.inputdir

    if args.diamond_dat == "None":
        print("[ERROR] diamond data must be set by -d option.")
        sys.exit()
    else:
        diamondfile = args.diamond_dat

    if args.basta_dat == "None":
        print("[ERROR] basta data must be set by -b option.")
        sys.exit()
    else:
        bastafile = args.basta_dat


    # Load config info
    config_ini = configparser.ConfigParser()
    configfile=join(str(pathlib.Path(__file__).parent.parent),'python_config.ini')
    config_ini.read(configfile, encoding='utf-8')

    if not os.path.exists(configfile):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), configfile)

    ENTREZ_EMAIL = config_ini['DEFAULT']['ENTREZ_EMAIL']
    ENTREZ_API = config_ini['DEFAULT']['ENTREZ_API']
    Nthreads = int(config_ini['DEFAULT']['THREADS'])

    # annotate GTDB info
    filelist = glob.glob(join(inputdir,"bowtie2_results","*","*_depth_summary.txt"))
    for f in filelist:
        annotate_to_depth(f,bastafile,diamondfile,join(pathlib.Path(f).parent,f.split("/")[-1].split("_")[0]+"_depth_annotate.txt"),ow,get_primary_id)

    filelist = glob.glob(join(inputdir,"bowtie2_results","*","*_depth_annotate.txt"))

    for i in filelist:
        create_summary_annotate(i,join(pathlib.Path(i).parent,i.split("/")[-1].split("_")[0]+"_annotate_summary.txt"),ow)
