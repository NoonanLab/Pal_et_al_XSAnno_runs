#!/bin/bash

cd /vast/palmer/scratch/noonan/ap2549/XSAnno/Sim_testing/
# Define the input and output file prefixes
inputs=("liftover.cDNA.hg38TopanTro6.hg38.fa" "liftover.cDNA.hg38TopanTro6.panTro6.fa")
outputs=("simReads.cDNA.hg38TopanTro6.hg38" "simReads.cDNA.hg38TopanTro6.panTro6")

# Loop through each input and output file prefix
for ((j=0; j<${#inputs[@]}; j++)); do
    input_file=${inputs[j]}
    output_prefix=${outputs[j]}
    
    # Loop from 1 to 10
    for i in {1..10}; do
        # Construct the full command with the current value of i and the current input/output prefix
        full_command="/vast/palmer/scratch/noonan/ap2549/simNGS/bin/simLibrary -r 101 -i 100 -x 10 -p ${input_file} |/vast/palmer/scratch/noonan/ap2549/simNGS/bin/simNGS -I -o \"fastq\" s_3_4x.runfile >${output_prefix}.${i}.fq"
        
        # Run the command
        eval $full_command
    done
done
