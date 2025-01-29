#!/bin/bash

usage_exit() {
    echo "Usage: $0 [-a assembly_file] [--filename get the path info of the generated bowtie2ref file] item ..." 1>&2
    exit 1
}

# get script (parent) directory
script_dir=$(cd $(dirname $0); cd .. ; pwd)

# load contig file
source "$script_dir"/CDSannotator.config

# initiallizing conda before using virtual environment
source "$conda_profile"



while getopts "a:-:" opt
do    
    case "$opt" in
        a) assembly_file=$OPTARG ;;

        -) 
            case $OPTARG in
            filename) stdout=1 ;;
             \?) usage exit ;;
        esac
        ;;
        \?) usage_exit ;;
    esac
done

bowtie2_build(){
    assembly_file="$1"
    prefix="$2"
    bowtie2-build -f "$assembly_file" "$bowtie2_refdir"/"$prefix"_bowtie2ref --quiet --threads 24
}


# main body
if [ ! -d "$bowtie2_refdir" ];then
    mkdir "$bowtie2_refdir"
fi


# creating bowtie2 ref file for ref genomes
prefix=$(echo "$assembly_file" | grep -o 'GC[A-Z]_[0-9]*.[0-9]*')
if [ -z "$prefix" ];then
    prefix=$(echo "$assembly_file" | awk -F"/" '{print $NF}' | awk -F"." '{print$(NF-1)}')
fi
bowtie2_build "$assembly_file" "$prefix"

if [ "$stdout" == 1 ]; then
    echo "$bowtie2_refdir"/"$prefix"_bowtie2ref >&1
fi