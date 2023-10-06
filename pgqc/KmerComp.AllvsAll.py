#!/usr/bin/env python3

### Authors: Max Marin (maximillian_marin@hms.harvard.edu)
# Purpose: Perform an all vs all comparison of the k-mer profiles of input sequences.

# Input Paths: 
### pan_genome_reference.fa

# Output Files:
### AvA.MaxJC.tsv

# Example Usage: KmerComp.AllvsAll.py -i pan_genome_reference.fa -o pan_genome_reference.AvA.MaxJC.tsv -k 31

__doc__ = """Perform an all vs all comparison of k-mer profiles of a `pan_genome_reference.fasta` output file of the Panaroo pipeline."""

import sys
import argparse
import time
from kmerlib import all_vs_all_kmer_MaxJC, read_kmers_from_file_ToHashesDict

import logging
# Set the logging level to INFO
logging.basicConfig(level=logging.INFO)



def main():
    parser = argparse.ArgumentParser(description="Calculate the maximum Jaccard Containment of k-mer content between all pairs of input nucleotide sequences.")
    parser.add_argument('-i', '--in_pg_ref', type=str, required=True, help="Input Panaroo pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo (FASTA)")
    parser.add_argument('-o', '--out_ava_tsv',type=str, required=True, help="All vs all comparison of sequence k-mer profiles using the maximum Jaccard Containment metric. (TSV)")
    parser.add_argument('-k', '--kmer_size',type=int, default=31, help="k-mer size (bp) to use for generating profile of unique k-mers for each sequence")

    args = parser.parse_args()

    ## 1) Set input parameters and PATHs ####

    input_PG_Ref_FA = args.in_pg_ref

    output_AvA_TSV = args.out_ava_tsv

    kmer_size = args.kmer_size
    
    ## 2) Parse and hash all k-mers for each representative nucleotide sequence
    logging.info(" Beginning parsing of input FASTA")

    start = time.time()

    Ref_DictOf_Hashes, Ref_DictOf_SeqLen = read_kmers_from_file_ToHashesDict(input_PG_Ref_FA, kmer_size)             

    All_SeqIDs = list(Ref_DictOf_Hashes.keys())

    end = time.time()
    time_diff = end - start
    logging.info(f" Time to parse and hash all k-mers: {round(time_diff, 2)} seconds")


    ## 3) Calculate the maximum Jaccard Containment (JC) between all pairs of sequences.
    ### NOTE: The maximum JC between sets a and b will always be symetrical, while JC is not

    start = time.time()

    PG_AvA_DF = all_vs_all_kmer_MaxJC(All_SeqIDs, Ref_DictOf_Hashes, Ref_DictOf_SeqLen) 

    end = time.time()
    time_diff = end - start
    logging.info(f" Time for all vs all comparison of k-mer profiles: {round(time_diff, 2)} seconds")

    ## 4) Output the AvA comparison table
    logging.info(f" Saving the All vs All comparison table (TSV) to: {output_AvA_TSV}")
    PG_AvA_DF.to_csv(output_AvA_TSV, sep = "\t", index = False)


if __name__ == "__main__":
    sys.exit(main())

