#!/bin/bash

# Usage: ./convert_gtf_to_bed.sh input.gtf.gz output.bed
# input can be input.gtf.gz or input.gtf

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.gtf.gz output.bed"
    exit 1
fi

input_gtf=$1
output_bed=$2

# Check if the input file is gzipped
if [[ $input_gtf == *.gz ]]; then
    zcat $input_gtf | awk -F '\t' '
    BEGIN { OFS = "\t" }
    {
        if ($0 ~ /^#/ || $3 != "exon") next
        split($9, attributes, "; ")
        gene_id = transcript_id = gene_name = transcript_name = exon_id = ""
        for (i in attributes) {
            split(attributes[i], kv, " ")
            gsub(/"/, "", kv[2])
            if (kv[1] == "gene_id") gene_id = kv[2]
            else if (kv[1] == "transcript_id") transcript_id = kv[2]
            else if (kv[1] == "gene_name") gene_name = kv[2]
            else if (kv[1] == "transcript_name") transcript_name = kv[2]
            else if (kv[1] == "exon_id") exon_id = kv[2]
        }
        annotation = gene_id "|" transcript_id "|" gene_name "|" transcript_name "|" exon_id
        print $1, $4-1, $5, annotation,"1", $7
    }' > $output_bed
else
    awk -F '\t' '
    BEGIN { OFS = "\t" }
    {
        if ($0 ~ /^#/ || $3 != "exon") next
        split($9, attributes, "; ")
        gene_id = transcript_id = gene_name = transcript_name = exon_id = ""
        for (i in attributes) {
            split(attributes[i], kv, " ")
            gsub(/"/, "", kv[2])
            if (kv[1] == "gene_id") gene_id = kv[2]
            else if (kv[1] == "transcript_id") transcript_id = kv[2]
            else if (kv[1] == "gene_name") gene_name = kv[2]
            else if (kv[1] == "transcript_name") transcript_name = kv[2]
            else if (kv[1] == "exon_id") exon_id = kv[2]
        }
        annotation = gene_id "|" transcript_id "|" gene_name "|" transcript_name "|" exon_id
        print $1, $4-1, $5, annotation,"1", $7
    }' $input_gtf > $output_bed
fi

echo "Conversion complete. BED file saved to $output_bed"
