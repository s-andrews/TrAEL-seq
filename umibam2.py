#!/usr/bin/env python

# Usage: module load python 
# python ./umibam2.py --verbose bamfile
import getopt
import sys
import os
import argparse
import gzip
import pysam
from argparse import RawTextHelpFormatter
import hashlib
import random

input_file = ""
verbose = False

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description = '''For deduplicating a bam file. \nIf the mapq value for a read is > 20, it is deduplicated by position and umi. \nIf the mapping quality is <= 20, only the UMI is used for deduplication. For 8bp UMIs, there would be a maximum of 65,536 (4^8) unique sequences. To increase this number we take the first 10 bases of the actual sequence and add that to the umi.''')
parser.add_argument('input_file', type=str, default="", help='sorted and indexed bam file containing reads to be deduplicated')
parser.add_argument('--verbose', default=False, action='store_true', help='verbose processing')
args=parser.parse_args()

input_file = args.input_file
output_file = input_file.replace("_sorted", "")
output_file = output_file.replace(".bam", "")
output_file = f'{output_file}_dedup_umibam2.1.bam'
trimming_report = input_file.replace(".bam", "")
trimming_report = trimming_report.replace("_sorted", "")
trimming_report = f'{trimming_report}_umibam2.1_dedup_report.txt'

samfile = pysam.AlignmentFile(input_file, "rb")
outfile = pysam.AlignmentFile(output_file, "wb", template=samfile)

print('\n-----------------------------------------------------------------------------------------------')
print(f'This only works for single ended data. For paired end, use the original UmiBam script for now.')
print('-----------------------------------------------------------------------------------------------\n')
print(f'Now processing {input_file}...' "\n")
print(f'Deduplication summary will be written to {trimming_report}.'"\n")

pos_umi_dup_count = 0
pos_not_umi_dup_count = 0
reads_processed = 0
first_occurrence = 0
umi_dup = 0
low_mapq = 0
low_mapq_keep = 0
lowqual_reads_processed = 0

pos_umi_dict = {}
umi_counts = {}
umis_written = set()
umis_written_lowqual = set()
dup_pos_dict = {}
all_lowqual = {}

all_chr = samfile.references

# The original UmiBam (this is umibam2) only uses start position for single ended data as it may have been trimmed to different lengths
# but then it does parse the cigar string. We're just going for the start position, and not dealing with paired end reads

# A read here is an AlignedSegment.
# https://pysam.readthedocs.io/en/latest/usage.html
for chr in all_chr:
    chr_dict = {}
    dup_pos_dict[chr] = set()
    all_reads_in_chr = samfile.fetch(reference=chr)
    for read in all_reads_in_chr:
        reads_processed+=1
        if reads_processed % 1000000 == 0:
            if args.verbose:
                print(f'{reads_processed:,} reads processed...')
        mapq = read.mapping_quality
        if mapq > 20: 
            ref_start = read.reference_start
            #print(f'reference_start =  {ref_start}')
            #print(f'query_sequence =  {read.query_sequence}')
            #print(f'query_alignment_sequence =  {read.query_alignment_sequence}')            
            if read.flag == 0: # read is on forward strand
                ref_start = ref_start
            elif read.flag == 16: # read is on reverse strand
                ref_start = read.reference_end # get the end of the mapped read
            else:
                print(f'Unknown flag value found {read.flag}')

            umi = read.query_name.split(":").pop()
            if ref_start in chr_dict: # duplicate read position found
                if umi in chr_dict[ref_start]: # duplicate position and umi found, adding umi
                    pos_umi_dup_count+=1 
                    chr_dict[ref_start][umi].append(read)
                    dup_pos_dict[chr].add(ref_start) # to keep track of the number of duplicated positions we find - this is a python set so won't allow duplicates                                     
                else:  # duplicate position, but new umi
                    chr_dict[ref_start][umi] = [read]                  
                    pos_not_umi_dup_count+=1

            else: # new read position, add the umi and read
                chr_dict[ref_start] = {}
                chr_dict[ref_start][umi] = [read]
                first_occurrence+=1

        else: # low quality read
            lowqual_reads_processed+=1
            umi = read.query_name.split(":").pop() # ignoring position for mapq <= 20, instead using UMI + first 10 bases of sequence
            query_seq =  read.query_sequence
            longer_umi = umi+query_seq[:10]
            
            if longer_umi in all_lowqual:
                all_lowqual[longer_umi].append(read) # this is appending to a list

            else: 
                all_lowqual[longer_umi] = [read] # adding new long umi

    # Go through all the high qual reads for each position. 
    # If the position is unique, (only 1 umi) write out the read. 
    # If the position is duplicated but umi is unique, write out the read.
    # If the position and umi are duplicated select a random read from the list.                 
    # The seed has been set dependent on the umi sequence, so the same read should be selected if the script is run multiple times.

    chosen_read = None
    for ref_start in chr_dict:
        for umi in chr_dict[ref_start]:
            no_of_reads = len(chr_dict[ref_start][umi])
            if no_of_reads > 1:
                umi_seed = int(hashlib.sha256(umi.encode('utf-8')).hexdigest(), 16) % 10**3
                random.seed(umi_seed)
                random_int = random.randint(0, no_of_reads-1)
                #print(f'random int = {random_int}')
                chosen_read = chr_dict[ref_start][umi][random_int]
            else: 
                chosen_read = chr_dict[ref_start][umi][0] # Only 1 read with that umi and position so we write it out
                
            outfile.write(chosen_read)
            umis_written.add(umi)


    # Go through all the low qual reads.
    # If the extended umi is unique, write out the read.
    # If the extended umi is duplicated, select a random read to write out.                 
    # The seed has been set dependent on the umi sequence, so the same read should be selected if the script is run multiple times.
chosen_read = None
for longer_umi in all_lowqual:
    no_of_reads = len(all_lowqual[longer_umi])
    if no_of_reads >= 1: 
        if no_of_reads > 1: 
            umi_seed = int(hashlib.sha256(longer_umi.encode('utf-8')).hexdigest(), 16) % 10**3
            random.seed(umi_seed)
            random_int = random.randint(0, no_of_reads-1)
            chosen_read = all_lowqual[longer_umi][random_int]
            removed_reads=no_of_reads-1
            umi_dup+=removed_reads
            
        else:
            chosen_read = all_lowqual[longer_umi][0] # There's only 1 read so we write it out 

        low_mapq_keep+=1
        outfile.write(chosen_read)
                

outfile.close()

if args.verbose:  
    print(f'Processed file {input_file}')
    print("\n"f'{first_occurrence+pos_not_umi_dup_count:,} deduplicated reads with mapq values > 20 have been written out.')
    print(f'Number of unique 8bp umis in high quality reads = {len(umis_written):,}')

# Count up how many duplicated positions we find - don't necessarily need to keep this in but it allows another comparison with the original umibam script.
dup_position_count = 0
for chr in dup_pos_dict:
    dup_position_count += len(dup_pos_dict[chr])

total_removed = umi_dup+pos_umi_dup_count
total_retained = first_occurrence+pos_not_umi_dup_count+low_mapq_keep

report_out = open(trimming_report, "w")

report_out.write("\n"f'Total reads processed: {reads_processed:,}')
report_out.write("\n"f'Total reads removed: {total_removed:,}, {(total_removed/reads_processed)*100:.1f}%')
report_out.write("\n"f'Total reads retained: {total_retained:,}, {(total_retained/reads_processed)*100:.1f}%'"\n")
report_out.write("\n"f'Number of unique 8bp umis in high quality reads = {len(umis_written):,}')
report_out.write("\n"f'Duplicated alignments of high quality reads were found at: {dup_position_count} different positions'"\n")
report_out.write("\n-------------more details-------------------------\n")
report_out.write("\n"f'Unique reads or first occurrences found: {first_occurrence:,}')
report_out.write("\n"f'Duplicated sequences with different UMIs and retained: {pos_not_umi_dup_count:,}')
report_out.write("\n"f'Duplicated sequences with identical UMIs and removed: {pos_umi_dup_count:,}')
report_out.write("\n"f'Low mapq reads : {lowqual_reads_processed:,}')
report_out.write("\n"f'Duplicated UMIs removed (that may have mapped to multiple regions): {umi_dup:,}')
report_out.write("\n"f'Low mapq reads where one UMI has been kept (that may have mapped to multiple regions): {low_mapq_keep:,}')

report_out.close()

if args.verbose:
    print("\n================================================\n")
    print(f'Total reads processed: {reads_processed:,}')
    print(f'Total reads removed: {total_removed:,}, {(total_removed/reads_processed)*100:.1f}%')
    print(f'Total reads retained: {total_retained:,}, {(total_retained/reads_processed)*100:.1f}%')
    print("\n"f'Duplicated alignments of high quality reads were found at: {dup_position_count:,} different positions'"\n")
    print("\n-------------more details-------------------------\n")
    print(f'Unique reads or first occurrences found: {first_occurrence:,}')
    print(f'Duplicated sequences with different UMIs and retained: {pos_not_umi_dup_count:,}')
    print(f'Duplicated sequences with identical UMIs and removed: {pos_umi_dup_count:,}')
    print(f'Low mapq reads : {lowqual_reads_processed:,}')
    print(f'Duplicated UMIs removed (that may have mapped to multiple regions): {umi_dup:,}')
    print(f'Low mapq reads where one UMI has been kept (that may have mapped to multiple regions): {low_mapq_keep:,}')
    print("\n===============================================\n")

