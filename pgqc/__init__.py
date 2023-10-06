#!/usr/bin/env python3

### Authors: Max Marin (maximillian_marin@hms.harvard.edu)
# Pan-genome QC toolkit (PQGC)


import sys
import argparse

from ._version import __version__

import logging
# Set the logging level to INFO
logging.basicConfig(level=logging.INFO)

from .ava import ava


def _ava_cli(args):
    ## 1) Set input parameters and PATHs ####
    input_PG_Ref_FA = args.in_pg_ref

    output_AvA_TSV = args.out_ava_tsv

    kmer_size = args.kmer_size

    ## Run the All vs All comparison function 
    PG_AvA_DF = ava(input_PG_Ref_FA, kmer_size)

    ## Output the AvA comparison table
    logging.info(f" Saving the All vs All comparison table (TSV) to: {output_AvA_TSV}")
    PG_AvA_DF.to_csv(output_AvA_TSV, sep = "\t", index = False)


def _pgcluster_cli(args):
    print("not yet implemented")


def main():
    parser = argparse.ArgumentParser(description="Toolkit for focused on augmenting Panaroo's pan-genome analysis with nucleotide sequence comparison.")
    sub_parser = parser.add_subparsers(required=True, help='Please select one of the pipelines of the PGQC toolkit.')

    ava_parser = sub_parser.add_parser("ava")
    ava_parser.add_argument('-i', '--in_pg_ref', type=str, required=True, help="Input Panaroo pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo (FASTA)")
    ava_parser.add_argument('-o', '--out_ava_tsv',type=str, required=True, help="All vs all comparison of sequence k-mer profiles. (TSV)")
    ava_parser.add_argument('-k', '--kmer_size',type=int, default=31, help="k-mer size (bp) to use for generating profile of unique k-mers for each sequence")
    ava_parser.set_defaults(func=_ava_cli)

    cluster_parser = sub_parser.add_parser("nscluster")
    cluster_parser.set_defaults(func=_pgcluster_cli)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())

