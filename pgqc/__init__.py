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
from .nscluster import run_nscluster
from .asm_gene_search import asmseqcheck_frompaths, get_AsmSeqCheck_QCStats


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


def _nscluster_cli(args):

    #run_nscluster()
    print("not yet implemented")


def _asmseqcheck_cli(args):

    ## 1) Set input parameters and PATHs ####
    
    # Define input paths
    input_PG_Ref_FA = args.in_pg_ref

    input_AsmFA_TSV = args.in_assemblies

    input_PresAbs_CSV = args.in_gene_matrix

    # Define output path
    output_PresAbs_WiDNASeqCheck = args.out_gene_matrix_wi_geneseqcheck
    
    # Set the minimum query coverage and sequence identity based on user input 
    min_query_cov = args.min_query_cov
    min_seq_id = args.min_seq_id

    ## 2) Run the assembly sequence check function ####

    Gene_PresAbs_WiAsmSeqCheck_DF = asmseqcheck_frompaths(input_PresAbs_CSV,
                                                          input_PG_Ref_FA,
                                                          input_AsmFA_TSV,
                                                          min_query_cov,
                                                          min_seq_id)

    # 3) Print the general QC Stats

    _ = get_AsmSeqCheck_QCStats(Gene_PresAbs_WiAsmSeqCheck_DF)


    # 4) Output TSV 
    Gene_PresAbs_WiAsmSeqCheck_DF.to_csv(output_PresAbs_WiDNASeqCheck,
                                         sep = "\t",
                                         index = False)



def main():
    parser = argparse.ArgumentParser(description="Toolkit for focused on augmenting Panaroo's pan-genome analysis with nucleotide sequence comparison.")
    sub_parser = parser.add_subparsers(required=True, help='Please select one of the pipelines of the PGQC toolkit.')

    asmseqcheck_parser = sub_parser.add_parser("asmseqcheck")
    asmseqcheck_parser.add_argument('-a', '--in_assemblies', type=str, required=True,
                                    help="Paths to input assemblies. (TSV)")

    asmseqcheck_parser.add_argument('-r', '--in_pg_ref', type=str, required=True,
                                    help="Input pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo/Roary (FASTA)")

    asmseqcheck_parser.add_argument('-m', '--in_gene_matrix', type=str, required=True,
                                    help="Input pan-genome gene presence/absence matrix. Typically output as `gene_presence_absence.csv` by Panaroo/Roary (CSV)")

    asmseqcheck_parser.add_argument('-o', '--out_gene_matrix_wi_geneseqcheck',type=str, required=True,
                                    help="Output pan-genome gene presence/absence matrix with updated gene presence/absence calls. (CSV). \n NOTE: 2 reflects that similar gene sequence is present at the nucleotide level (CSV)")

    asmseqcheck_parser.add_argument('-c', '--min_query_cov', type=float, default=0.9,
                            help="Minimum query coverage to classify a gene as present within an assembly (0-1)")

    asmseqcheck_parser.add_argument('-i', '--min_seq_id', type=float, default=0.9,
                            help="Minimum sequence identity to classify a gene as present within an assembly (0-1)")

    asmseqcheck_parser.set_defaults(func=_asmseqcheck_cli)


    ava_parser = sub_parser.add_parser("ava")
    ava_parser.add_argument('-i', '--in_pg_ref', type=str, required=True,
                            help="Input Panaroo pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo (FASTA)")
    ava_parser.add_argument('-o', '--out_ava_tsv',type=str, required=True,
                            help="All vs all comparison of sequence k-mer profiles. (TSV)")
    ava_parser.add_argument('-k', '--kmer_size',type=int, default=31,
                            help="k-mer size (bp) to use for generating profile of unique k-mers for each sequence (Default: 31))")
    ava_parser.set_defaults(func=_ava_cli)

    cluster_parser = sub_parser.add_parser("nscluster")
    cluster_parser.set_defaults(func=_nscluster_cli)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())

