# CDSannotator
Python scripts for generating and extracting bacterial CDSs from fastq files,
and performing taxonomic annotations based on the homology search against NCBI nr and GTDB databases.

The codes are specifically designed for single-microbial vesicle(MV) sequencing analysis study,
which deposited in biorxiv doi: https://doi.org/10.1101/2024.09.18.613607 (Takano et al., 2024).

The main program is "./annotation_pipeline/CDSannotator.sh", which automatically cleans, generates contigs, extracts CDSs,
and then annotates them under used defined configurations.
Before running, all the parameters and environment paths should be set in "python_config.ini" and "CDSannotator.config".
