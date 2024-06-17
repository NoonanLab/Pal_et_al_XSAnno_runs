#!/bin/bash

# Function to sort a CSV file while keeping the header intact
sort_file() {
    local input_file=$1
    local output_file=$2
    
    # Clean up Windows-style line endings
    dos2unix "$input_file"
    
    # Extract the header
    head -n 1 "$input_file" > "$output_file"
    
    # Sort the rest of the file and append to the output
    tail -n +2 "$input_file" | sort >> "$output_file"
}

# Define input and temporary sorted file names
file1="exon_count_matrix_hs.csv"
file2="exon_count_matrix_pt.csv"
sorted_file1="sorted_file1.csv"
sorted_file2="sorted_file2.csv"
output_file="exon_count_matrix_hs_pt.txt"

# Sort the files
sort_file "$file1" "$sorted_file1"
sort_file "$file2" "$sorted_file2"

# Merge the sorted files by the first column
# Convert comma-separated to tab-separated and use join
join -t, -1 1 -2 1 <(awk 'BEGIN {FS=OFS=","} {print $0}' "$sorted_file1") <(awk 'BEGIN {FS=OFS=","} {print $0}' "$sorted_file2") | tr ',' '\t' > "$output_file"

# Clean up
rm "$sorted_file1" "$sorted_file2"
