#!/bin/bash

INPUT_DIR="data"

#READS="ecoli-reference.fasta-reads-len1000-cov10.fasta"
#READS="ecoli-reference.fasta-reads-len10000-cov40.fasta"
#READS="SRR10971019-ecoli.fastq_20x_sampled.fasta"
READS="SRR10971019-ecoli.fastq_40x_sampled.fasta"
#READS="SRR10971019-ecoli.fastq_60x_sampled.fasta"

READS_PATH=$INPUT_DIR"/"$READS

K=101
A=2
#INITIAL_SLOTS=41000000 # 0.8
#INITIAL_SLOTS=86000000 # 0.8
#INITIAL_SLOTS=55000000 # 0.6
#INITIAL_SLOTS=115000000 # 0.6
SIZE_FLOOR=55000000
KMERS_PATH="output/myht-"$READS"-"$K"kmers-a"$A".txt"
KMERS_SORTED_PATH="output/myht-"$READS"-"$K"kmers-a"$A"-sorted.txt"
OUTPUT_1=$KMERS_PATH"-STATS-1.txt"
OUTPUT_2=$KMERS_PATH"-STATS-2.txt"

#perf record --output "PERF-RECORD" --call-graph dwarf -e cpu-clock,faults ./build/KMerHashtable -k $K -a $A -s $SIZE_FLOOR -p $READS_PATH -v
/usr/bin/time -v ./build/KMerHashtable -k $K -a $A -s $SIZE_FLOOR -p $READS_PATH -v -o $KMERS_PATH #> $OUTPUT_1
#/usr/bin/time -v sort -o $KMERS_SORTED_PATH $KMERS_PATH #> $OUTPUT_2

