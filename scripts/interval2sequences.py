#!/usr/bin/env python3

import sys
import os
import subprocess
import tempfile
from intervaltree import Interval, IntervalTree

def usage(prog_name):
    print(f"Usage: {prog_name} <file.2bit> <file.annotation> <exonic|genomic>")
    sys.exit(1)

def read_interval_file(interval_file):
    intervals = {}
    with open(interval_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            #if len(fields) < 8:
                #print(f"Error parsing line: {line.strip()}")
                #continue  # Skip lines that do not have enough columns
            transcript_id = fields[0]
            chromosome = fields[1]
            strand = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            num_exons = int(fields[5])
            exon_starts = list(map(int, fields[6].split(',')))
            exon_ends = list(map(int, fields[7].split(',')))
            
            if transcript_id not in intervals:
                intervals[transcript_id] = []

            for exon_start, exon_end in zip(exon_starts, exon_ends):
                intervals[transcript_id].append((chromosome, exon_start, exon_end, strand))
    return intervals

def write_targets_file(intervals, interval_type):
    targets_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt')
    if interval_type == 'genomic':
        for transcript_id, exons in intervals.items():
            for exon in exons:
                chromosome, start, end, strand = exon
                targets_file.write(f"{chromosome}:{start}-{end}\n")
    elif interval_type == 'exonic':
        for transcript_id, exons in intervals.items():
            for exon in exons:
                chromosome, start, end, strand = exon
                targets_file.write(f"{chromosome}:{start}-{end}\n")
    targets_file.close()
    return targets_file.name

def run_twoBitToFa(two_bit_file, targets_file):
    cmd = f"twoBitToFa {two_bit_file} stdout -noMask -seqList={targets_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    sequences = result.stdout.split('>')[1:]  # Skip the first empty split
    parsed_sequences = {}
    for seq in sequences:
        header, sequence = seq.split('\n', 1)
        parsed_sequences[header.strip()] = sequence.replace('\n', '')
    return parsed_sequences

def main():
    if len(sys.argv) != 4:
        usage(sys.argv[0])
    
    two_bit_file = sys.argv[1]
    interval_file = sys.argv[2]
    interval_type = sys.argv[3]

    if interval_type not in ["genomic", "exonic"]:
        usage(sys.argv[0])
    
    intervals = read_interval_file(interval_file)
    targets_file = write_targets_file(intervals, interval_type)
    sequences = run_twoBitToFa(two_bit_file, targets_file)

    output_file = f"./blat/{interval_file.split('/')[-1].replace('.Interval', '.fa')}"
    with open(output_file, 'w') as out_f:
        if interval_type == "genomic":
            for transcript_id, exons in intervals.items():
                for exon in exons:
                    chromosome, start, end, strand = exon
                    seq_key = f"{chromosome}:{start}-{end}"
                    sequence = sequences.get(seq_key, '')
                    out_f.write(f">{transcript_id}|{chromosome}|{strand}|{start}|{end}\n{sequence}\n")
        elif interval_type == "exonic":
            for transcript_id, exons in intervals.items():
                full_sequence = ''
                for exon in exons:
                    chromosome, start, end, strand = exon
                    seq_key = f"{chromosome}:{start}-{end}"
                    full_sequence += sequences.get(seq_key, '')
                # Find the overall start and end for the combined sequence
                combined_start = min(exon[1] for exon in exons)
                combined_end = max(exon[2] for exon in exons)
                out_f.write(f">{transcript_id}|{chromosome}|{strand}|{combined_start}|{combined_end}\n{full_sequence}\n")

    os.remove(targets_file)

if __name__ == "__main__":
    main()
