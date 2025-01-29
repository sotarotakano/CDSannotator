#!/bin/bash

usage_exit() {
    echo "Usage: $0 [-d targetdir] [-i sample_id] [-r reference_file] item ..." 1>&2
    exit 1
}

# get script (parent) directory
script_dir=$(cd $(dirname $0); cd .. ; pwd)

# load contig file
source "$script_dir"/CDSannotator.config

# initiallizing conda before using virtual environment
source "$conda_profile"


while getopts "d:i:r:o:b:" opt
do    
    case "$opt" in
        d) target_dir=$OPTARG ;;
        i) sample_id=$OPTARG ;;
        r) reference_file=$OPTARG ;;
        o) output_dir=$OPTARG ;; 
        b) bowtie_option=$OPTARG ;;
        \?) usage_exit ;;
    esac
done


cd "$target_dir" || exit

bowtie2_mapping(){
    id=$1
    reference_file=$2
    savedir=$3
    option=$4

    fq1=./clean_fastq/${id}_1.fq.gz
    fq2=./clean_fastq/${id}_2.fq.gz

    if [ -f "${fq1}" ] && [ -f "${fq2}" ] ; then
        echo "Loading... ${fq1} and ${fq2}..."
        # There are several options: --sensitive-local, --very-sensitive-local ,--sensitive, --very-sensitive
        bowtie2 -x "$reference_file" -1 "${fq1}" -2 "${fq2}" -p $thread "${option}" | samtools view -@ $thread -bh > "$id".bam || return
        samtools sort -@ $thread -m 2G "$id".bam > sorted.bam
        samtools depth sorted.bam -aa > "$id"_depth.txt
        rm -f sorted.bam "$id".bam
        
        mv "$id"_depth.txt "$savedir"
      
    else
        echo "Oops! either of ${fq1} and ${fq2} does not exist..."
    fi
}


if [ -z "$output_dir" ] ;  then
    output_dir=$(basename "$reference_file" | grep -o '[a-zA-Z0-9]*.[0-9]*.*.bowtie2')
fi

if [ -z "$bowtie_option" ] ;  then
    bowtie_option="--sensitive-local"
fi

# main body
if [ ! -d "$output_dir" ];then
    mkdir "$output_dir"
fi

# mapping reads for bowtie2 ref files
bowtie2_mapping "$sample_id" "$reference_file" "$output_dir" "$bowtie_option"

