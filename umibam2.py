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

input_file = ""
verbose = False

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description = '''For deduplicating a bam file. \nIf the mapq value for a read is > 20, it is deduplicated by position and umi. \nIf the mapping quality is <= 20, only the UMI is used for deduplication. For 8bp UMIs, there would be a maximum of 65,536 (4^8) unique sequences. To increase this number we take the first 10 bases of the actual sequence and add that to the umi.''')
parser.add_argument('input_file', type=str, default="", help='sorted and indexed bam file containing reads to be deduplicated')
parser.add_argument('--verbose', default=False, action='store_true', help='verbose processing')
args=parser.parse_args()

input_file = args.input_file
output_file = input_file.replace(".bam", "")
output_file = f'{output_file}_dedup_umibam2.bam'
trimming_report = input_file.replace(".bam", "")
trimming_report = f'{trimming_report}_umibam2_dedup_report.txt'

samfile = pysam.AlignmentFile(input_file, "rb")
outfile = pysam.AlignmentFile(output_file, "wb", template=samfile)

print('\n-----------------------------------------------------------------------------------------------')
print(f'This only works for single ended data. For paired end, use the original UmiBam script for now.')
print('-----------------------------------------------------------------------------------------------\n')
print(f'Now processing {input_file}...' "\n")
print(f'Deduplication summary will to written to {trimming_report}.'"\n")

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
chr_dict = {}
low_mapq_reads = []
umis_written = set()
umis_written_lowqual = set()
dup_pos_dict = {}

all_chr = samfile.references

# UmiBam only uses start position for single ended data as it may have been trimmed to different lengths
# but then it does parse the cigar string. We're just going for the start position, then using the ref_end method to get the start of a read on the reverse strand.

# A read here is an AlignedSegment.
# https://pysam.readthedocs.io/en/latest/usage.html
for chr in all_chr:
    chr_dict[chr] = {}
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
            if ref_start in chr_dict[chr]: # duplicate read position found
                if umi in chr_dict[chr][ref_start]: # duplicate position and umi found, not writing out
                    pos_umi_dup_count+=1                  
                    dup_pos_dict[chr].add(ref_start) # to keep track of the number of duplicated positions we find - this is a python set so won't allow duplicates                                     
                else:
                    chr_dict[chr][ref_start].append(umi) # umi is different, adding to set and writing out
                    pos_not_umi_dup_count+=1
                    outfile.write(read)
                    umis_written.add(umi)
            else:
                chr_dict[chr][ref_start] = [umi] # adding new read position
                first_occurrence+=1
                outfile.write(read)
                umis_written.add(umi)
        else:
            umi = read.query_name.split(":").pop() # ignoring position for mapq <= 20, instead using UMI + first 10 bases of sequence
            query_seq =  read.query_sequence
            longer_umi = umi+query_seq[:10]
            
            lowqual_reads_processed+=1
            if lowqual_reads_processed % 1000000 == 0:
                if args.verbose:
                    print(f'{lowqual_reads_processed:,} lowqual_reads processed...')
            if longer_umi in umis_written_lowqual:
                umi_dup+=1 # discard 
            else: 
                umis_written_lowqual.add(longer_umi)
                low_mapq_keep+=1
                outfile.write(read)

outfile.close()

if args.verbose:  
    print(f'Processed file {input_file}')
    print("\n"f'{first_occurrence+pos_not_umi_dup_count:,} deduplicated reads with mapq values > 20 have been written out.')
    print(f'Number of unique 8bp umis in high quality reads = {len(umis_written):,}')
    print(f'Number of unique 18bp umis in low quality reads = {len(umis_written_lowqual):,}')

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
report_out.write("\n"f'Number of unique 18bp umis in low quality reads = {len(umis_written_lowqual):,}')
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

