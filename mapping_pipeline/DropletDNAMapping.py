# ! usr/bin/env python

__author__ = 'Sotaro Takano'
__version__ = '1.0.0'

'''This scripts mapped the sequence reads of MV droplet sequencing and visualizing.
Also, extracting the enriched genomic region on the chromosome of parental strain.'''

import re
import numpy as np
import sys
import pandas as pd
import os
from os.path import join
import copy
import math
import glob
import gc
import argparse
import subprocess
import configparser
import pathlib
import errno

from scipy.stats import hypergeom
from scipy.stats import binom_test
import matplotlib.pyplot as plt 
from matplotlib import cm
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.colors import LinearSegmentedColormap


def create_bowtie2ref(reffile):
    # Load bowtie2 reference file of user-defined reference genome 
    # The program search for the "genomic.fna" file. This is usually available at NCBI.
    bowtie2ref = list({re.search('.+_bowtie2ref',x).group(0) for x in os.listdir(bowtie2_dir) 
                      if reffile.split("/")[-1].split(".")[0] in x})
    if len(bowtie2ref) > 0:
        bowtie2ref = join(bowtie2_dir,bowtie2ref[0])
    else:
        # if bowtie2 ref does not exist. 
        if os.path.exists(reffile):
            reffile = [reffile]
        else:
            reffile = list({re.search('.+_genomic.fna',x).group(0) for x in os.listdir(Refseq_dir) 
                            if reffile.split("/")[-1].split(".")[0] in x})
        if len(reffile) > 0:
            returns = subprocess.run([bt_build_script,"-a",join(Refseq_dir,reffile[0])], 
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            bowtie2ref = [x for x in returns.stdout.decode("utf-8").split("\n") if "Users" in x][0]
            print(returns.stderr.decode("utf-8"))
        else:
            sys.stderr.write('No bowtie2 reference file and reference GCF file found for %s'%reffile)
    
    return bowtie2ref

def search_positive_droplets(maindir,tophit=True,GTDB=True,params="length",min_params=1000):
    # Check whether more than 10% (default) of reads are mapped to target reference genome
    align_data = glob.glob(join(maindir,"bowtie2_results/*/*_besthits.txt"))

    MV_idlist = []
    for besthit in align_data:
        summary = pd.read_csv(besthit, sep='\t',index_col = 0)
        if tophit:
            if GTDB:
                try:
                    summary_GTDB = summary.loc[summary["GTDB"]==1,:]
                    summary = summary_GTDB.sort_values(params,ascending=False).iloc[0,:]
                except:
                    summary = summary.sort_values(params,ascending=False).iloc[0,:]
            else:
                summary = summary.sort_values(params,ascending=False).iloc[0,:]

            if reffile in summary["filename"] and summary[params] > min_params:
                MV_idlist.append(besthit.split("/")[-1].split("_")[0])
            else:
                continue
        else:
            target = summary.loc[{x for x in summary.index if type(x) == str and 
                                  reffile in summary.loc[x,"filename"]},:]
            if target.shape[0] > 0:
                target_params = np.array([float(x) for x in target[params]])
                max_params = max(target_params)
                if max_params > min_params:
                    MV_idlist.append(besthit.split("/")[-1].split("_")[0])
    
    return MV_idlist


def heatmap_bowtie2(maindir,bowtie2ref,droplets="all",bppercell=20):
    
    '''This script running bowtie2 and create a summary of mapped regions on the genome
    across the droplets. It returns pandas.dataframe format summary object.'''

    # Load all possible droplet lists
    dropletlist = {x.split("/")[-1].split("_")[0] for x in 
                   glob.glob(join(maindir,"done_fastp","*.fq.gz"))}
    if droplets == "all":
        pass
    elif len(droplets) > 0 and type(droplets) == list:
        dropletlist = [x for x in dropletlist if len([y for y in droplets if y in x]) > 0]
    print("%i droplets are going to be analyzed..."%len(dropletlist))
    

    # Running bowtie2 to map sequence reads
    for id in dropletlist:
        bowtie2_data = join(savedir,id + "_depth.txt")
        if not os.path.exists(bowtie2_data):
            returns = subprocess.run([bowtie2_script,"-d",maindir,"-i",id,"-r",bowtie2ref,],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if returns.returncode == 1:
                print(returns.stderr.decode("utf-8"))
                continue

        print("%s finished"%id)
    
    # Create summary of bowtie2 results
    for bowtie2_data in [x for x in os.listdir(savedir) if "_depth.txt" in x]:
        id = bowtie2_data.split("_")[0]

        if len([x for x in dropletlist if id in x]) < 1:
            pass

        with open(join(savedir,bowtie2_data), mode='r') as f:
            lines = f.readlines()
            mapdata = np.array([int(x.split('\t')[2]) for x in lines])
            mapdata[mapdata > 0] = 1
        
        g_pos = np.zeros(int(len(mapdata)/bppercell))
    
        for i in range(len(g_pos)-1):
            if (i+1)*bppercell < len(mapdata):
                g_pos[i] = np.mean(mapdata[int(i*bppercell):int((i+1)*bppercell)])
            else:
                g_pos[i] = np.mean(mapdata[int(i*bppercell):-1])
    
        g_pos = pd.Series(g_pos,index=[x for x in range(len(g_pos))],name = id)    
        
        if "G_pos" not in locals():
            G_pos = g_pos
        else:
            G_pos = pd.concat([G_pos,g_pos],axis = 1)
        gc.collect()
        #os.remove(bowtie2_data)
    G_pos = G_pos.T.sort_index()
    
    return G_pos


def extract_enriched_region(G_pos,ref_fna,savedir,threshold = 0.5, exclude_pos=[],
                            pvalue=0.01,bppercell=1000,method="top",ratio=0.25,
                            bp_per_line = 50):

    # Create summary file for extracting enriched regions
    summary = pd.DataFrame(0,index = G_pos.columns.values,columns=["MV_count","CDF"])
    
    total_MV = G_pos.shape[0]    # the number of total droplets analyzing MVs
    total_cells = total_MV*G_pos.shape[1]    # the number of total cells

    # Sorting mapped cells (i.e., more than 80% (default) of genomic region in a single cell
    # is mapped by the sequence reads).
    raveled = np.ravel(G_pos) 
    total_mapped = len(raveled[raveled > threshold])
    
    # Counting positive droplets (MV) in each genomic region
    for i in G_pos.columns.values:
        MVcount = G_pos.loc[G_pos[i] > threshold,:].shape[0]
        summary.loc[i,"MV_count"] = MVcount
        summary.loc[i,"frequency"] = MVcount/total_MV
        
        # Statistical screening of significantly enriched regions
        if method == "hygm":
            if MVcount == 0:
                MVcount = 1
            cdf = 1 - hypergeom.cdf(MVcount-1, total_cells, total_MV, total_mapped)
            summary.loc[i,"CDF"] = cdf
        elif method == "binom":
            if MVcount == 0:
                MVcount = 1
            cdf = binom_test(x=MVcount, n=total_MV, p=total_mapped/total_cells, alternative='greater')
            summary.loc[i,"CDF"] = cdf
    
    # These methods are also screening but without statistical test. 
    # "top" : just picking out the cells from the top by the frequency of positive droplets
    if method == "top":
        summary = summary.sort_values("MV_count",ascending = False)
        topgroup= [x for x in summary.loc[:int(len(summary)*ratio),:].sort_index().index]
    # "frequency" : Using the frequency of positive droplets as a threshold for siginficant enrichment
    elif method == "frequency":
        summary = summary.sort_values("MV_count",ascending = False)
        topgroup= [x for x in summary.sort_index().index if summary.loc[x,"frequency"] > ratio]
    # "negative" : Return the non-detected regions (i.e., the frequency is below the threshold.)
    elif method == "non-detected":
        summary = summary.sort_values("MV_count",ascending = False)
        topgroup= [x for x in summary.sort_index().index if summary.loc[x,"frequency"] < ratio]
    elif method =="hygm" or method=="binom":
        summary = summary.sort_values("CDF")
        topgroup= [x for x in summary.sort_index().index if summary.loc[x,"CDF"] < pvalue]
    else:
        print("[ERROR] Invalid method for the statistical screening of enriched regions.")
        sys.exit()
        
    enriched = [int(x) for x in G_pos.loc[:,topgroup].columns.values]

    # (optional) If there are genomic positions where the user don't want to include the analysis,
    # This position is excluded from the analysis.
    if len(exclude_pos) > 0:
        enriched = [x for x in enriched if x not in exclude_pos]
    
    print("%f of genomic region is overrepresented."%(len(enriched)/G_pos.shape[1]))
    #enriched_bp =[]
    
    #for i in enriched:
    #    enriched_bp = enriched_bp + [x for x in range(bppercell*i,bppercell*(i+1))]
    
        
    print("Found %s as a save directory..."%savedir)
    
    bowtie2_data = glob.glob(join(savedir,"*_depth.txt"))[0]
    
    print("%s used for the bowtie2_data..."%bowtie2_data)
    
    with open(bowtie2_data, mode='r') as f:
        lines = f.readlines()
    mapdata = [x for x in lines]
    
    
    # enriched_bp (list)
    # (position in dataframe, contig_name, position in reference genome)
    enriched_bp = []
    for i in enriched:
        bp_info = mapdata[min(bppercell*i,len(mapdata)-1)]
        contig_idx = bp_info.split('\t')[0]
        pos_idx = bp_info.split('\t')[1]
        enriched_bp.append((i,contig_idx,pos_idx))

    # refseq: Reference sequences (list)
    with open(ref_fna, mode='r') as f:
        refseq = [x for x in f.readlines()]
    
    # contig_pos: Start of the each contig position (dict)
    contigs_pos = {x:refseq.index(x) for x in refseq if ">" in x}
    target_contig = [refseq.index(x) for x in refseq if mapdata[0].split('\t')[0] in x][0]
    
    # Create fna file (enrichment sequences)
    outputfile = join(savedir, 
                      os.path.basename(ref_fna).replace("ref","").replace("fna","") 
                      + "_enriched_sequence_genome.fna")
    with open(outputfile, mode='wt') as nf:
        nf.write("")
    
    
    for n,i in enumerate(enriched_bp):
        # contig_s: contig start position(line) including the sequence of interest
        # Here, i[1] is the contig name, and searching for contigs including the sequences of interest
        contig_s = [refseq.index(x) for x in refseq if i[1].split('\t')[0] in x][0] + 1
        
        # next_contig: the next contig start position (line) 
        # contig_e: contig end position (line) including the sequence of interest
        next_contig = [refseq.index(x) for x in refseq[contig_s:] if ">" in x]
        if next_contig:
            contig_e = next_contig[0] - 1
        else:
            contig_e = len(refseq)
        
        # Here, i[2] is the contig position
        contig_interest = ''.join([x.strip('\n') for x in refseq[contig_s:contig_e]])
        enriched_seq = contig_interest[max(0,int(i[2])):min(len(contig_interest),int(i[2])+bppercell)]
        
        if enriched_bp[n][0] - enriched_bp[n-1][0] > 1 or n == 0:
            with open(outputfile, mode='a') as nf:
                nf.write('>NODE'+ str(i[0]) + '_' + str(i[1]) + '_' + str(i[2]) + '\n')
        
        
        if bppercell%bp_per_line != 0:
            print("[ERROR] Please specify a value that is divisible by the number of characters.")
            sys.exit()
        
        for p in range(int(len(enriched_seq)/bp_per_line)):
            with open(outputfile, mode='a') as nf:
                try:
                    nf.write(enriched_seq[p*bp_per_line:(p+1)*bp_per_line] + '\n')
                except:
                    nf.write(enriched_seq[p*bp_per_line:] + '\n')
    
    summary.to_csv(join(savedir, os.path.basename(ref_fna).replace("ref","").replace("fna",""))
                   + "_enriched_summary.csv")
    return summary



# Visualizing merged summary results
def draw_merged_summary(G_pos,summary,maindir,threshold=0.5,bppercell=20,pvalue=1e-06,xnticks="auto",xdtick="auto",
                        xbps="auto",ymax_h="auto",method="CDF",ratio=0.25,ylim_f = [],
                        SaveFig=True,log=False,vmax = 1,clustering = True):
    
    print("Creating a figure summary...")
    a = int(G_pos.shape[1]/1000)
    tick_width = 0.5
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['axes.linewidth'] = tick_width
    
    # creating subplot
    fig, ax = plt.subplots(3,1,dpi=300,figsize= (3,3),
                            gridspec_kw={'height_ratios': [4, 10, 1]})
    a_factor = 0.58
    fontsize_axis = 6
    fontsize_label = 8

    # setting xdticks
    if xdtick == "auto":
        xdtick = np.power(10,int(round(np.log10(G_pos.shape[1]*bppercell))/3)*3)
    
    xdtickpower = int(np.log10(xdtick))
    
    if xdtickpower == 3:
        scale = "k"
    elif xdtickpower == 6:
        scale = "M"
    elif xdtickpower == 0:
        scale = ""
    
    if xbps == "auto":
        xbps = math.floor((bppercell*G_pos.shape[1])/xdtick*5)/5

    if xnticks == "auto":
        xbps_m = int(xbps*10)
        xnticks = xbps_m
        for n in range(1,int(xbps_m/3)):
            if xbps_m%n == 0:
                xnticks = int(xbps_m/n + 1)

    # plot1: frequency
    droplet_frequency = []
    for i in G_pos.columns.values:
        pos_i = G_pos.loc[G_pos[i] > threshold,:]
        droplet_frequency.append(pos_i.shape[0]/G_pos.shape[0])
    ax[0].plot(np.linspace(0,G_pos.shape[1]-1,G_pos.shape[1]),np.array(droplet_frequency),
               color='black',lw=0.5)
    ax[0].set_xticks([x*xdtick/bppercell for x in np.linspace(0,xbps,xnticks)])
    ax[0].set_xlim([0,len(droplet_frequency)])
    ax[0].set_xticklabels([])
   
    ymax = math.ceil(max(droplet_frequency)*5)/5
    
    nyticks = 2
    for n in range(2,math.ceil(ymax*10/3)+1):
        if (ymax*10)%n == 0:
            nyticks = int(ymax*10/n + 1)

    ax[0].set_yticks([round(x*ymax/(nyticks-1),1) for x in range(nyticks)])
    ax[0].set_yticklabels([round(x*ymax/(nyticks-1),1) for x in range(nyticks)],fontsize=fontsize_axis)
    if len(ylim_f) > 1:
        ax[0].set_ylim(ylim_f)
    ax[0].set_ylabel("Frequency",fontsize = fontsize_label)
    ax[0].set_ylim([-0.01,0.25])

    # plot2 Heatmap of All data
    aspect = (G_pos.shape[1]/G_pos.shape[0])*a_factor
    
    # Clustering
    G_pos_bin = copy.copy(G_pos)
    G_pos_bin[G_pos_bin < max(np.array(np.mean(G_pos)))/10] = 0
    G_pos_bin[G_pos_bin >= max(np.array(np.mean(G_pos)))/10] = 1
    if G_pos.shape[0] > 2 and clustering:
        dist = pdist(G_pos_bin,metric="euclidean")
        ylinkage = linkage(dist,method = 'average',metric='euclidean')
        ydendro = dendrogram(ylinkage, orientation = "left", no_labels = True,
                            distance_sort = 'ascending', no_plot = True)

        G_pos = G_pos.loc[[G_pos.index[i] for i in ydendro['leaves']]]
        #G_pos = G_pos.iloc[::-1]
    
    if ymax_h == "auto":
        ymax_h = round(max(np.array(np.mean(G_pos)))/2,-2)

    if log:
        G_pos[G_pos < 1e-06] = 0.1
        G_pos = np.log10(G_pos)
        ymax = np.log10(ymax)
    
    
    G_pos[G_pos < threshold] = 0
    G_pos[G_pos >= threshold] = 1
    
    print(G_pos.shape)
    ax[1].imshow(G_pos, cmap = 'binary', vmin = 0, vmax = vmax, 
                 aspect = aspect, interpolation="nearest")

    ax[1].set_xticks([x*xdtick/bppercell for x in np.linspace(0,xbps,xnticks)])
    ax[1].set_xticklabels([])
    
    if G_pos.shape[0] < 20:
        yscale = 1
    else:
        yscale = 10

    ymax = math.floor(G_pos.shape[0]/yscale)*yscale
    

    nyticks = 4
   
    for n in range(2,math.ceil(ymax/yscale/3)+1):
        if (ymax/yscale)%n == 0:
            nyticks = int(ymax/yscale/n + 1)
    print([ymax,nyticks])
    ax[1].set_yticks([int(ymax/(nyticks-1))*x for x in range(nyticks)])
    ax[1].set_yticklabels([int(ymax/(nyticks-1))*x for x in range(nyticks)],fontsize=fontsize_axis)
    ax[1].set_ylabel("Droplets",fontsize = fontsize_label)

    # plot3: enrichment
    enrichment = pd.DataFrame(0,index=["enrichment"],columns=G_pos.columns.values)
    sort_summary = summary.sort_values("frequency",ascending=False)
    topgroup= [x for x in sort_summary.iloc[:int(len(sort_summary)*ratio),:].sort_index().index] 

    for i in summary.index:
        if method == "CDF":
            if summary.loc[i,"CDF"] < pvalue:
                enrichment.loc["enrichment",i] = 1
                for j in range(int(i)-a,int(i)+a):
                    if str(j) in enrichment.columns.values:
                        enrichment.loc["enrichment",str(j)] = 1
                    else:
                        print(j)
        elif method == "frequency":
            if summary.loc[i,"frequency"] > ratio:
                enrichment.loc["enrichment",i] = 1
                for j in range(int(i)-a,int(i)+a):
                    enrichment.loc["enrichment",j] = 1
        elif method == "non-detected":
            if summary.loc[i,"frequency"] < ratio:
                enrichment.loc["enrichment",i] = 1
                for j in range(int(i)-a,int(i)+a):
                    enrichment.loc["enrichment",j] = 1
        elif method == "top":
            if  i in topgroup:
                enrichment.loc["enrichment",i] = 1
                for j in range(int(i)-a,int(i)+a):
                    enrichment.loc["enrichment",j] = 1

    
    aspect = (enrichment.shape[1]/enrichment.shape[0])*(a_factor/10)

    ax[2].imshow(enrichment, cmap = "Reds", vmin = 0, vmax = 1, 
                 aspect = aspect, interpolation="nearest")

    ax[2].set_xlabel("Position in draft genome ({power}bp)".format(power=scale),fontsize=fontsize_label)
    ax[2].set_xticks([x*xdtick/bppercell for x in np.linspace(0,xbps,xnticks)])
    ax[2].set_xticklabels([round(x,1) for x in np.linspace(0,xbps,xnticks)], fontsize = fontsize_axis)
    ax[2].set_yticklabels([])
    ax[2].set_yticks([])
   
    if SaveFig:
        [x.tick_params(width=tick_width) for x in ax]
        plt.rcParams['font.family'] = 'Arial'
        fig.savefig(join(maindir, "summary_heatmap.pdf"),dpi = 300, 
                    bbox_inches="tight",pad_inches = 0.0)
    return enrichment


#################################
#### Main body of the script ####
#################################
if __name__ == "__main__":
    
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-d","--data_dir", type=str, default="None", help="data directory folder")
    p.add_argument("-r","--reffile", type=str, default="None", help="reference fna file")
    p.add_argument("-e","--enrichment", type=str, default="binom",
                   help="method for screening enriched regions on the genome. The default is 'binom'.")
    p.add_argument("-b","--bppercell", type=int, default=1000, help="base pair per cell (default: 1000)")
    p.add_argument("-i","--idlist",type=str,default="None",
                   help="droplet id list used for analysis (.txt filename should be specified)")
    p.add_argument("-p","--pvalue",type=float,default=0.01,
                   help="p-value for screening enriched regions on the genome (default:1000)")
    args = p.parse_args()
        
    if args.data_dir != "None":
        maindir = args.data_dir
        savedir = args.data_dir
    
    if args.reffile != "None":
        reffile = args.reffile
    else:
        print('Error: no reference file set', file=sys.stderr)
        sys.exit(1)

    if args.idlist != "None":
        if args.idlist.split(".")[-1] == "txt":
            with open(args.idlist, mode="r") as f:
                idlist = [x.strip("\n") for x in f.readlines()]
        elif args.idlist == "annotated":
            idlist = search_positive_droplets(maindir,tophit=True,GTDB=True,params="length",min_params=1000)
        else:
            print("[ERROR] idlist should be specified by txt format.")
            sys.exit()
    else:
        idlist = "all"

    bppercell = args.bppercell

    # Load config info
    config_ini = configparser.ConfigParser()
    configfile=join(str(pathlib.Path(__file__).parent.parent),'python_config.ini')
    config_ini.read(configfile, encoding='utf-8')

    if not os.path.exists(configfile):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), configfile)
        
    # Bowtie2 scripts and datadir
    bowtie2_script = config_ini['DEFAULT']['BOWTIE2_SCRIPT']
    bt_build_script = config_ini['DEFAULT']['BOWTIE2_BUILD_SCRIPT']
    bowtie2_dir = config_ini['DEFAULT']['BOWTIE2_DIR']
    Refseq_dir = config_ini['DEFAULT']['REFSEQ_DIR']
    Nthreads = int(config_ini['DEFAULT']['THREADS'])
    
    # Check bowtie2ref of target reference genome
    bowtie2ref = create_bowtie2ref(reffile)

    # create savedirectory
    savedir = join(maindir,os.path.basename(bowtie2ref).replace("ref",""))
    if not os.path.exists(savedir):
        os.mkdir(savedir)

    print(savedir)
    pvalue=args.pvalue

    aligndata_path = join(savedir,re.search("GC[A-Z]_[0-9]+.[0-9]",reffile).group(0)+"_gpos.txt")
    # Load align data if exists.
    if os.path.exists(aligndata_path):
        aligndata = pd.read_csv(aligndata_path, index_col=0, header=0, sep="\t")
    else:
        # if not exists, create from first
        aligndata = heatmap_bowtie2(maindir,bowtie2ref,bppercell=bppercell,droplets=idlist)
        aligndata.to_csv(aligndata_path, sep="\t")
    
    summary = extract_enriched_region(aligndata, reffile, savedir,
                                      threshold = 0.8, bppercell = bppercell, method = args.enrichment, pvalue=pvalue,
                                      exclude_pos = [],bp_per_line=80)

    if args.enrichment == "hygm" or args.enrichment == "binom":
        method_summary = "CDF"
    else:
        method_summary = args.enrichment

    enrichment = draw_merged_summary(aligndata,summary,maindir,pvalue=pvalue,threshold=0.8,
                                     bppercell=bppercell,method=method_summary,
                                     xnticks=4,xdtick="auto",xbps="auto",ylim_f=[-0.05,0.65],
                                     SaveFig=True,log=False,vmax = 1)
