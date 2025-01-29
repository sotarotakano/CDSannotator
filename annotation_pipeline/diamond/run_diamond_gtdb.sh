#!/bin/bash
usage_exit() {
    echo "Usage: $0 [-i targetfile] item ..." 1>&2
    exit 1
}

# get script (parent) directory
script_dir=$(cd $(dirname $0); cd .. ; pwd)

# load contig file
source "$script_dir"/CDSannotator.config

# initiallizing conda before using virtual environment
source "$conda_profile"

while getopts "i:" opt
do    
    case "$opt" in
        i) targetfile=$OPTARG ;;
        \?) usage_exit ;;
    esac
done


db_dir=$GTDB_DATABASE


run_diamond_gtdb_ffn(){
    inputfile="$1"
    prefix=$(echo "$inputfile" | awk -F'[./]' '{print $(NF-1)}')
    curdir=$(dirname "$inputfile")
    outfile="${curdir}/diamond_gtdb_results/${prefix}_output.txt"

    if [ ! -d "${curdir}/diamond_gtdb_results" ];then
        mkdir "${curdir}/diamond_gtdb_results"
    fi
    
    if [ ! -d "${curdir}/tmp_files_gtdb" ];then
        mkdir "${curdir}/tmp_files_gtdb"
    fi
    
    echo "Checking the preexisting output file..."
    if [ -f "$outfile" ];then
        echo "$outfile exists."
        echo "The diamond result (GTDB) of ${inputfile} already exist. we can skip this process..."
        return 0
    fi
    
    seq_n=$(\seqkit stats -T "$inputfile" | tail -n 1 | awk -F "\t" '{print $4}')
    if [ "$seq_n" -gt 10000 ]; then
        split_n=$((seq_n/10000))
    else
        split_n=1
    fi


    if [ ! -d "$inputfile".split ];then
        seqkit split -p $split_n "$inputfile" || exit
    fi

    cd "$inputfile".split || exit

    for splitfile in $(\find . -maxdepth 1 -name "*.ffn"); do
        splitprefix=$(echo "$splitfile" | awk -F'[./]' '{print $(NF-1)}')
        tempout="${splitprefix}_temp_gtdb.txt"
        echo "$splitfile"
        if [ ! -f "$tempout"  ] || [ ! -s "$tempout" ] ;then
            diamond blastx -d "$db_dir" -q "$splitfile" --evalue 1e-10 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \
            --sensitive > "$tempout" || return 0
        fi
    done


    for splitfile in $(\find . -maxdepth 1 -name "*.faa"); do
        splitprefix=$(echo "$splitfile" | awk -F'[./]' '{print $(NF-1)}')
        tempout="${splitprefix}_temp_gtdb.txt"
        echo "$splitfile"
        if [ ! -f "$tempout"  ] || [ ! -s "$tempout" ] ;then
            diamond blastp -d "$db_dir" -q "$splitfile" --evalue 1e-10 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \
            --sensitive > "$tempout" || return 0
        fi
    done

    find . -maxdepth 1 -name "*_temp_gtdb.txt" | xargs -I% cat % > "$outfile"
    cd .. || exit

    echo "[OUTPUT] $outfile generated..."
    mv "$inputfile".split ./tmp_files_gtdb -f
}

run_diamond_gtdb_ffn "${targetfile}"