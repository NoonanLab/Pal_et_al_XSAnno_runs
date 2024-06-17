#!/bin/bash

# Load necessary modules
module load HISAT2
module load SAMtools
module load StringTie

# Change to the appropriate directory
cd ~/palmer_scratch/XSAnno/Sim_testing/

# Define the arrays for the different species
a=("hg38" "panTro6")
b=("hs" "pt")

# Loop over the species and the iterations
for ((s=0; s<${#a[@]}; s++)); do
    for i in {1..10}; do
        # Run HISAT2
        hisat2 --dta-cufflinks -p 10 -x /vast/palmer/scratch/noonan/ap2549/XSAnno/Hisat2_indexes/${a[$s]}_Hisat2_index/${a[$s]} -q simReads.cDNA.hg38TopanTro6.${a[$s]}.${i}.fq -S hisat_alignment/${b[$s]}${i}.sam

        # Sort SAM file using SAMtools
        samtools sort hisat_alignment/${b[$s]}${i}.sam -o hisat_alignment/${b[$s]}${i}_sort.sam

        # Run StringTie
        stringtie hisat_alignment/${b[$s]}${i}_sort.sam -o stringtie/${b[$s]}${i}_stringTie.gtf -e -G blatFiltered.hg38TopanTro6.exon.${a[$s]}.gtf -A stringtie/${b[$s]}${i}_stringTie.tsv
    done
done
