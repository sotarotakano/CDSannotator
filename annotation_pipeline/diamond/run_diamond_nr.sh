#!/bin/bash
usage_exit() {
    echo "Usage: $0 [-i inputdir] item ..." 1>&2
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
        i) input_file=$OPTARG ;;
        \?) usage_exit ;;
    esac
done


db_dir=$NR_DATABASE

run_diamond_nr_ffn(){
    inputfile="$1"
    echo "[DIAMOND] Loading $inputfile"
    prefix=$(echo "$inputfile" | awk -F'[./]' '{print $(NF-1)}')
    curdir=$(dirname "$inputfile")
    outfile="${curdir}/diamond_nr_results/${prefix}_output.txt"

    if [ ! -d "${curdir}/diamond_nr_results" ];then
        mkdir "${curdir}/diamond_nr_results"
    fi
    
    if [ ! -d "${curdir}/tmp_files_nr" ];then
        mkdir "${curdir}/tmp_files_nr"
    fi

    echo "[DIAMOND] Checking the preexisting output file..."
    if [ -f "$outfile" ];then
        echo "[DIAMOND] $outfile exists."
        echo "[DIAMOND] The diamond result of ${inputfile} already exist. we can skip this process..."
        return 0
    fi
    
    seq_n=$(\seqkit stats -T "$inputfile" | tail -n 1 | awk -F "\t" '{print $4}')

    if [ "$seq_n" -gt 6000 ]; then
        split_n=$((seq_n/6000))
    else
        split_n=1
    fi

    if [ ! -d "$inputfile".split ];then
        seqkit split -p $split_n "$inputfile" || exit
    fi

    cd "$inputfile".split || exit

    for splitfile in $(\find . -maxdepth 1 -name "*.ffn"); do
        splitprefix=$(echo "$splitfile" | awk -F'[./]' '{print $(NF-1)}')
        tempout="${splitprefix}_temp.txt"
        echo "[DIAMOND] run diamond for $splitfile"
        if [ ! -f "$tempout"  ] || [ ! -s "$tempout" ] ;then
            diamond blastx -d "$db_dir" -q "$splitfile" --evalue 1e-10 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \
            --sensitive > "$tempout" || return 0
        fi
    done
    
    tempout_all="${curdir}/diamond_nr_results/${prefix}_temp.txt"
    find . -maxdepth 1 -name "*_temp.txt" | xargs -I% cat % > "$outfile"
    cd .. || exit
    
    echo "$outfile"
    mv "$inputfile".split/ ./tmp_files_nr -f
}


# main body
if [ -z "$input_file" ]; then
    echo "Please input input file with -i option."
    exit
fi

target_dir=$(dirname "$input_file")
cd "$target_dir" || exit

echo "Running diamond_nr..."
run_diamond_nr_ffn "${input_file}"