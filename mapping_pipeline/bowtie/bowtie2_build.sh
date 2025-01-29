#!/bin/bash

usage_exit() {
    echo "Usage: $0 [-a assembly_file] [--filename get the path info of the generated bowtie2ref file] item ..." 1>&2
    exit 1
}

# initiallizing conda before using virtual environment
source /home/soutacano/anaconda3/etc/profile.d/conda.sh


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
    bowtie2-build -f "$assembly_file" /home/soutacano/bowtie2ref/"$prefix"_bowtie2ref --quiet --threads 24
}


# main body
if [ ! -d "/home/soutacano/bowtie2ref" ];then
    mkdir /home/soutacano/bowtie2ref
fi


# creating bowtie2 ref file for ref genomes
prefix=$(echo "$assembly_file" | grep -o 'GC[A-Z]_[0-9]*.[0-9]*')
if [ -z "$prefix" ];then
    prefix=$(echo "$assembly_file" | grep -o '[a-zA-Z0-9]*.[0-9]*.fna')
fi
bowtie2_build "$assembly_file" "$prefix"

if [ "$stdout" == 1 ]; then
    echo /home/soutacano/bowtie2ref/"$prefix"_bowtie2ref >&1
fi