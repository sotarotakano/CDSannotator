#!/bin/bash
# CDSannotator coded by Sotaro Takano
# Updated in 2023/11/1


usage_exit() {
    echo "Usage: $0 [-i inputdir] [-a adapter type]  [-d decoy database (.fasta)] \
    [--nonhuman_only clean human reads] [--skip_fastp skip fastp process] [--skip_nr skip diamond against nr database] ..." 1>&2
    exit 1
}

set -o pipefail

# get script directory
script_dir=$(cd $(dirname $0); pwd)
source "$script_dir"/CDSannotator.config

# initiallizing conda before using virtual environment
source "$conda_profile"

non_human=0
skip_fastp=0
adapter=0

while getopts "i:a:d:-:" opt
do    
    case "$opt" in
        i) target_dir=$OPTARG ;;
        a) adapter=$OPTARG ;;
        d) decoy_data=$OPTARG ;;
        -) 
            case $OPTARG in
            nonhuman_only) non_human=1 ;;
            skip_fastp) skip_fastp=1;;
            skip_nr) skip_nr=1;;
            \?) usage exit ;;
        esac
        ;;
        \?) usage_exit ;;
    esac
done


# Function to create contigs from fastq
fastq_to_spades () {
    prefix="$1"
    non_human="$2"
    skip_fastp="$3"
   
    if [ -n "$(find . -maxdepth 1 -name "${prefix}*.fastq.gz")" ]; then
        suffix="fastq.gz"
    elif [ -n "$(find . -maxdepth 1 -name "${prefix}*.fq.gz")" ]; then
        suffix="fq.gz"
    else
        echo "invalid extension of ${prefix}"
        return 0
    fi
    
    input1=${prefix}_1.${suffix}
    input2=${prefix}_2.${suffix}

    # Create a directory
    if [ ! -d "done_fastp" ];then
        mkdir done_fastp
    fi

    output1=./done_fastp/${prefix}_1_out.fq.gz
    output2=./done_fastp/${prefix}_2_out.fq.gz
   

    # Quality control by fastp
    if [ "$skip_fastp" == 1 ]; then
        echo "[fastp] Skip trimming and quality control by fastp."
        mv "$input1" "$output1"
        mv "$input2" "$output2"

    else
        if [ -f "${output1}" ];then
            echo "The fastp-processed fastq of ${prefix} already exist. we can skip this process..."
        else
            echo "[fastp] Running fastp..."
            if [ -f "${input1}" ] && [ -f "${input2}" ] ; then
                echo "Loading ${input1} and ${input2}..."
                fastp -i  "${input1}" -I "${input2}" -3 -o "${output1}" -O "${output2}" \
                -h report.html -j report.json -n 10 -t 1 -T 1 -l 20 -w 16 \
                --adapter_sequence "${adapter1}" --adapter_sequence_r2 "${adapter2}"
            elif [ -f "${input1}" ] ; then
                echo "Found ${input1} but not ${input2}..."
                fastp -i  "${input1}" -3 -o "${output1}" \
                -h report.html -j report.json -n 10 -t 1 -T 1 -l 20 -w 16 --adapter_sequence "${adapter1}"
            elif [ -f "${input2}" ] ; then
                echo "Found ${input2} but not ${input1}..."
                fastp -i  "${input2}" -3 -o "${output2}" \
                -h report.html -j report.json -n 10 -t 1 -T 1 -l 20 -w 16 --adapter_sequence_r2 "${adapter2}"
            else
                echo "Oops! both ${input1} and ${input2} does not exist..."
                return 0
            fi
        fi

    fi

    # Optional: Clean human reads 
    if [ "$non_human" == 1 ];then
        if [ ! -d "non_human_fastq" ];then
            mkdir non_human_fastq
        fi
        
        #clean human reads
        if [ ! -f ./non_human_fastq/"${prefix}"-umhuman_1.fq.gz ] || [ ! -f ./non_human_fastq/"${prefix}"-umhuman_1.fq.gz ] ; then
            "$script_dir"/clean/clean_human_reads.sh -1 "${output1}" -2 "${output2}" -o "non_human_fastq" -r $bowtie2_humanref || exit
        else
            echo "[CLEAN_HUMAN] The cleaned fastq of ${prefix} already exist. we can skip this process..."
        fi
        output1=./non_human_fastq/${prefix}-umhuman_1.fq.gz
        output2=./non_human_fastq/${prefix}-umhuman_2.fq.gz
    fi

    # Optional: Cleaned by decoy data (should be specified as fasta (contigs))
    if [ -n "$decoy_data" ]; then
        if [ ! -d "clean_fastq" ];then
            mkdir clean_fastq
        fi
        # eliminate contamination reads      
        if [ ! -f ./clean_fastq/"${prefix}"-cleaned_1.fq.gz ] || [ ! -f ./clean_fastq/"${prefix}"-cleaned_2.fq.gz ] ; then
            "$script_dir"/clean/clean_contam_reads.sh -1 "${output1}" -2 "${output2}" -o "clean_fastq" -r "$decoy_data" || exit
        else
            echo "[CLEAN_DECOY] The cleaned fastq of ${prefix} already exist. we can skip this process..."
        fi
    fi

    # Finally specify which fastq files were used for the downstream analysis
    if [ ! -d "clean_fastq" ];then
        mkdir clean_fastq
    fi

    # Copy the fastq files to where the programs will seach for
    if [ -n "$decoy_data" ]; then
        cp ./clean_fastq/${prefix}-cleaned_1.fq.gz ./clean_fastq/${prefix}_1.fq.gz -f
        cp ./clean_fastq/${prefix}-cleaned_2.fq.gz ./clean_fastq/${prefix}_2.fq.gz -f
    elif [ "$non_human" == 1 ]; then
        cp ./non_human_fastq/${prefix}-umhuman_1.fq.gz ./clean_fastq/${prefix}_1.fq.gz -f
        cp ./non_human_fastq/${prefix}-umhuman_2.fq.gz ./clean_fastq/${prefix}_2.fq.gz -f
    else
        cp ./done_fastp/${prefix}_1_out.fq.gz ./clean_fastq/${prefix}_1.fq.gz -f
        cp ./done_fastp/${prefix}_2_out.fq.gz ./clean_fastq/${prefix}_2.fq.gz -f
    fi

     # Specify file names for SPAdes 
    input1=./clean_fastq/${prefix}_1.fq.gz
    input2=./clean_fastq/${prefix}_2.fq.gz
    output=./${prefix}_spades   
    
    
    # running fastqc
    echo "[fastqc] running fastqc..."
    if [ -f "${input1}" ] || [ -f "${input2}" ] ; then
        if [ -z "$(find ./output_fastqc -maxdepth 2 -name "${prefix}*fastqc.html")" ] && [ "$skip_fastp" == 0 ] ; then
            echo running fastqc qcqcqcq
            echo $skip_fastp
            fastqc "${input1}"
            fastqc "${input2}"
            # moving fastqc file to other folder
            find ./clean_fastq -maxdepth 1 -type f -name "${prefix}*_fastqc*" | xargs -I% mv -f % output_fastqc/
            #mv "report.html" "./output_fastqc/${prefix}_report.html" 
            #rm "report.json"
        else
            echo "The output fastqc of ${prefix} already exist. we can skip this process..."
        fi
    fi


    if [ -f "./contigs/${prefix}_contigs.fasta" ];then
        echo "The contig of ${prefix} already exist. we can skip this process..."
        return 0
    fi
    
    # Assembly of reads by SPAdes
    if [ -f "${input1}" ] && [ -f "${input2}" ] ; then
        echo "Loading ${input1} and ${input2}..."
        spades.py -1  "${input1}" -2 "${input2}" -k auto -t $thread --disable-rr -o "${output}" --meta || return 0
    elif [ -f "${input1}" ] ; then
        echo "Found ${input1} but not ${input2} so assembly will be performed using single-end..."
        spades.py -s "${input1}"--sc --careful -t $thread --disable-rr -o "${output}"
    elif [ -f "${input2}" ] ; then
        echo "Found ${input2} but not ${input1} so assembly will be performed using single-end..."
        spades.py -s "${input2}" --sc --careful -t $thread --disable-rr -o "${output}"
    else
        echo "Oops! both ${input1} and ${input2} does not exist..."
        return 0
    fi

    if [ -d "${output}" ] ; then
        find ./"${prefix}"* -type f -name "contigs.fasta" | xargs -I% cp -f % ./contigs/"${prefix}_contigs.fasta"
        find . -type d -name "*${prefix}*_spades" | xargs -I% mv -f % output_spades/
    else
        echo "[ERROR] output files of SPAdes does not exist."
    fi

}


## Function to run prokka based on contigs
contigs_to_prokka(){
    inputfile="$1"
    prefix=$(echo "$inputfile" | awk -F"/" '{print $NF}' | awk -F"_" '{print $1}' | sort -u)
    
    if [ -f ./prokka_faa/"${prefix}_prokka.faa" ];then
        echo "The prokka result of ${inputfile} already exist. we can skip this process..."
        return 0
    fi
    
    if [ ! -d "prokka_faa" ];then
        mkdir prokka_faa
    fi

    # optional
    cat "$inputfile" | seqkit seq -m $min_contig_length > tmpfile && mv tmpfile "$inputfile"

    prokka "$inputfile" --outdir ./"$prefix"_prokka --locustag "$prefix" --prefix "$prefix" --force --compliant --cpus $thread || return 0

    find ./"${prefix}"* -type f -name "*.faa" | xargs -I% cp % ./prokka_faa/"${prefix}_prokka.faa"

    if [ ! -d "prokka_results" ];then
        mkdir prokka_results
    fi
    
    find . -type d -name "${prefix}*_prokka" -maxdepth 1 | xargs -I% mv % prokka_results/
}




# main body

cd "$target_dir" || exit

# Setting adapter sequences
if [ "$adapter" == "nextera" ];then
    adapter1=CTGTCTCTTATACACATCT
    adapter2=CTGTCTCTTATACACATCT
elif [ "$adapter" == "truseq" ]; then
    adapter1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    adapter2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
elif [ "$adapter" == "bgi" ]; then
    adapter1=AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA
    adapter2=AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG
elif [ -z "$adapter" ]; then
    adapter1=CTGTCTCTTATACACATCT
    adapter2=CTGTCTCTTATACACATCT
else
    echo "[ERROR] Invalid adapter info."
    exit 1
fi


if [ -z "$target_dir" ]; then
    echo "Please input target directory with -i option."
    exit
fi

if [ ! -d "contigs" ];then
    mkdir contigs
fi

if [ ! -d "output_spades" ];then
    mkdir output_spades
fi


if [ ! -d "output_fastqc" ];then
    mkdir output_fastqc
fi

conda activate $ASSEMBLYenv

for file in $(\find . -maxdepth 1 -name "*.fastq.gz" -or -name "*.fq.gz" | awk -F"/" '{print $NF}' | awk -F"_" '{print $1}' | sort -u); do
    fastq_to_spades "${file}" "${non_human}" "${skip_fastp}"
done

conda deactivate


conda activate $PROKKAenv
for file in $(\find ./contigs -maxdepth 2 -name "*contigs.fasta"); do
    contigs_to_prokka "${file}"
done
conda deactivate

#>>>>>> Clustering by cd-hit
expid=$(echo "$target_dir" | awk -F'/' '{print $NF}' | sed -e 's/_/-/g')
cdhitoutput="${expid}_prokkaffn_98qc.ffn"
if [ ! -f "$cdhitoutput" ] || [ ! -s "$cdhitoutput" ]; then
    find ./prokka_results -maxdepth 2 -name "*.ffn" | xargs -I% cat % | sed s/'*'/''/g > prokka_ffn_concat.fasta
    if [ -s "prokka_ffn_concat.fasta" ]; then
        echo "[CD-HIT] Running CD-HIT..."
        cd-hit -i prokka_ffn_concat.fasta -c 0.98 -T 10 -M 16000 -o "$cdhitoutput" -s 0.5 -aS 0.9 -d 30 -T 16 || exit
    else
        echo "[CD-HIT] inputfile ./prokka_ffn_concat.fasta is empty."
        echo "Terminating this process..."
        exit
    fi
fi

if [ ! -f "$cdhitoutput" ]; then 
    echo "${cdhitoutput} is empty. Something wrong... terminating the process."
    rm "$cdhitoutput" prokka_ffn_concat.fasta "$cdhitoutput".clstr
    exit 
fi
#<<<<<<


#>>>>>> Set up diamond input and output files
if [ ! -d  "$diamond_outputdir" ]; then
    mkdir "$diamond_outputdir" || exit
fi

cdhitfile="$diamond_outputdir/${expid}_prokkaffn_98qc.ffn"
if [ ! -f "$cdhitfile" ] || [ ! -s "$cdhitfile" ]; then
    cp "${expid}_prokkaffn_98qc.ffn" $diamond_outputdir
fi
#<<<<<<


diamondfile="$diamond_outputdir/diamond_nr_results/${expid}_prokkaffn_98qc_output.txt"
basta_output="$diamond_outputdir/basta_output/basta_${expid}_output.txt"
if [ "$skip_nr" == 1 ]; then
    echo "[DIAMOND] DIAMOND aginst nr database is skipped."
    echo "" > "$diamondfile"
    echo "" > "$basta_output"
else
    #>>>>>> Run DIAMOND (nr)
    if [ -f "$diamondfile" ]; then
        echo "[DIAMOND] DIAMOND output file already exists. ${diamondfile}"
        if [ ! -s "$diamondfile" ]; then
            echo "[DIAMOND] output file is empty...running DIAMOND again..."
            #rm "$diamondfile"
            echo "[DIAMOND] running DIAMOND." 
            "$script_dir"/diamond/run_diamond_nr.sh -i "$cdhitfile" || exit
        else
            echo "[DIAMOND] Skip this process."
        fi
    else
        echo "[DIAMOND] running DIAMOND." 
        "$script_dir"/diamond/run_diamond_nr.sh -i "$cdhitfile" || exit
        # To check the generated output by diamond
        diamondfile=$(\find "$diamond_outputdir"/diamond_nr_results -name "${expid}_*")
        echo "[DIAMOND] DIAMOND output file generated! ${diamondfile}"
    fi


    if [ ! -f "$diamondfile" ]; then 
        echo "[DIAMOND] DIAMOND output file is empty..."
        echo "Terminating this process..."
        exit
    fi
    #<<<<<<


    #>>>>> run BASTA
    conda activate $BASTAenv

    if [ -f "$basta_output" ]; then
        echo "[BASTA] BASTA output file already exists. ${basta_output}"
        if [ ! -s "$basta_output" ]; then
            echo "[BASTA] output file is empty...running BASTA again..."
            rm "$basta_output"
            echo "running BASTA." 
            basta sequence "$diamondfile" "$basta_output" prot -l 100 -p 80 -i 80 -b True -z True || exit
        else
            echo "[BASTA] Skip this process."
        fi
    else
        echo "running basta." 
        basta sequence "$diamondfile" "$basta_output" prot -l 100 -p 80 -i 80 -b True -z True || exit
        echo "BASTA output file generated! $basta_output"
    fi

    if [ ! -s "$basta_output" ]; then 
        echo "[BASTA] BASTA output file is empty..."
        echo "Terminating this process..."
        exit
    fi
    conda deactivate
fi

#<<<<<


#>>>> Run bowtie2mapping to CDS
clstrfile="${expid}_prokkaffn_98qc.ffn.clstr"
if [ -f "$clstrfile" ] && [ -s "$clstrfile" ]; then
    "$script_dir"/mapping/read_mapping_prokka.sh -i "$target_dir" -b "$basta_output" -c "$clstrfile" -d "$diamondfile"
else
    echo "clstr output from CD-HIT for ${cdhitfile} is missing or empty..."
    echo "terminate this process..."
    exit
fi
#<<<<

#>>>> run diamond for GTDB
cdhitfile="${diamond_outputdir}/${expid}_prokkaffn_98qc.ffn"
if [ ! -f "$cdhitfile" ] || [ ! -s "$cdhitfile" ]; then
    cp "${expid}_prokkaffn_98qc.ffn" $diamond_outputdir
fi

# run diamond gtdb
diamond_gtdb_output="${diamond_outputdir}/diamond_gtdb_results/${expid}_prokkaffn_98qc_output.txt"
if [ ! -f "$diamond_gtdb_output" ] || [ ! -s "$diamond_gtdb_output" ]; then
    "$script_dir"/diamond/run_diamond_gtdb.sh -i "$cdhitfile" || exit
else
    echo "$diamond_gtdb_output exists. Skip this process..."
fi
mv "$cdhitfile" "${diamond_outputdir}/done"
python "$script_dir"/annotation/GTDB_getTaxonomy.py -i "$target_dir" -d "$diamond_gtdb_output" -mp 80 -ml 100
#<<<<
