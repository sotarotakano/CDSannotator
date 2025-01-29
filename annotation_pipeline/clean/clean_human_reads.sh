#!/bin/bash

usage_exit() {
    echo "Usage: $0 [-1 input fastqfile 1] [-2 input fastqfile 2] [-f reference fna file for cleaning] 
    [-o output_dir] [-r reference_file] [-b option of bowtie2] [-m method either of 'bowtie'(default) or 'samtools'] item ..." 1>&2
    exit 1
}

# get script (parent) directory
script_dir=$(cd $(dirname $0); cd .. ; pwd)

# load contig file
source "$script_dir"/CDSannotator.config

# initiallizing conda before using virtual environment
source "$conda_profile"


while getopts "1:2:r:o:b:m:-:" opt
do    
    case "$opt" in
        1) inputfile1=$OPTARG ;;
        2) inputfile2=$OPTARG ;;
        r) reference_file=$OPTARG ;;
        o) output_dir=$OPTARG ;; 
        b) bowtie_option=$OPTARG ;;
        m) method=$OPTARG ;;
        -) 
            case $OPTARG in
            get_output) output=1 ;;
            overwrite) ow=1 ;;
             \?) usage exit ;;
        esac
        ;;
        \?) usage_exit ;;
    esac
done

#cd "$target_dir" || exit


bowtie2_mapping_human(){
    fq1=$1
    fq2=$2
    reference_file=$3
    savedir=$4
    option=$5
    method=$6

    id=$(echo "$fq1" | awk -F"/" '{print $NF}' | awk -F"_" '{print $1}') 
    outfq1="./${savedir}/${id}-umhuman_1.fq.gz" 
    outfq2="./${savedir}/${id}-umhuman_2.fq.gz" 

    if [ -f "${outfq1}" ];then
        echo "Cleaned ${id} reads already exist. we can skip this process..."
        return 0
    fi

    if [ -f "${fq1}" ] && [ -f "${fq2}" ] ; then
        echo "Loading... ${fq1} and ${fq2}..."
        if [ "$method" == "bowtie" ] ; then
            bowtie2 -x "$reference_file" -1 "${fq1}" -2 "${fq2}" -p "$thread" "${option}" --un-conc-gz "$id"-umhuman_%.fq.gz --fr -S /dev/null || return
            find . -name "${id}-umhuman_*.fq.gz" | xargs -I% mv % "$savedir"
        elif [ "$method" == "samtools" ] ; then
            bowtie2 -x "$reference_file" -1 "${fq1}" -2 "${fq2}" -p "$thread" "${option}" | samtools view -@ "$thread" -bh -f 4 | samtools sort -@ "$thread" -m 2G > "$id"_umhuman.bam || return
            samtools fastq -@ "$thread" "$id"_umhuman.bam -1 "$outfq1" -2 "$outfq2" -0 /dev/null -s /dev/null -n
        else
            echo "Please specify bowtie or samtools by --method option"
        fi
    else
        echo "Oops! either of ${fq1} and ${fq2} does not exist..."
    fi

    #create summary
    echo id reads1 reads2 unmapped1 unmapped2 | sed s/' '/$'\t'/g > "$savedir"/"$id"_summary.txt
    reads_1=$(seqkit stats "$fq1" | awk -F" " '{print $4}' | tail -n 1 | sed 's/,//g')
    reads_2=$(seqkit stats "$fq2" | awk -F" " '{print $4}' | tail -n 1 | sed 's/,//g')
    
    unmapped_1=$(seqkit stats "$outfq1" | awk -F" " '{print $4}' | tail -n 1 | sed 's/,//g')
    unmapped_2=$(seqkit stats "$outfq2" | awk -F" " '{print $4}' | tail -n 1 | sed 's/,//g')

    echo -n -e "$id\t" "$reads_1\t" "$reads_2\t" "$unmapped_1\t" "$unmapped_2\n" >> "$savedir"/"$id"_summary.txt

}


if [ -z "$output_dir" ] ;  then
    output_dir="unmapped_human"
fi

if [ -z "$bowtie_option" ] ;  then
    bowtie_option="--sensitive"
fi

if [ -z "$reference_file" ] ;  then
    reference_file="$bowtie2_refdir"/GCF_000001405.40_bowtie2ref
fi

if [ -z "$method" ] ;  then
    method="bowtie"
fi

# main body
if [ ! -d "$output_dir" ];then
    mkdir "$output_dir"
fi

# mapping reads for bowtie2 ref files
#for file in $(\find . -maxdepth 1 -name "*.fastq.gz" | awk -F'[._/]' '{print $3}' | sort -u); do
#    bowtie2_mapping_human "$file" "$reference_file" "$output" "$output_dir" "$bowtie_option" "$method"
#done

bowtie2_mapping_human "$inputfile1" "$inputfile2" "$reference_file" "$output_dir" "$bowtie_option" "$method"

# Third step: check the number of mapped and unmapped reads

#check_mapped_reads "$output_dir"