#!/usr/bin/env python3

### Authors: Max Marin (maximillian_marin@hms.harvard.edu)
# Purpose: Perform an all vs all comparison of the k-mer profiles of input sequences.

# Input Paths: 
### pan_genome_reference.fa

# Output Files:
### AvA.MaxJC.tsv

# Example Usage: KmerComp.AllvsAll.V3.py -i pan_genome_reference.fa -o pan_genome_reference.AvA.MaxJC.tsv -k 31


__doc__ = """Perform an all vs all comparison of k-mer content of a `pan_genome_reference.fasta` output file of the Panaroo pipeline."""

import sys
import argparse
import pandas as pd
import numpy as np
import time
import screed
import mmh3


###### Define Functions ######

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
        
    return kmers


def read_kmers_from_file(filename, ksize):
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence
        
        kmers = build_kmers(sequence, ksize)
        all_kmers += kmers

    return all_kmers


def hash_kmer(kmer):
    # calculate the reverse complement
    rc_kmer = screed.rc(kmer)
    
    # determine whether original k-mer or reverse complement is lesser
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
        
    # calculate murmurhash using a hash seed of 42
    hash = mmh3.hash64(canonical_kmer, 42)[0]
    if hash < 0: hash += 2**64
        
    # done
    return hash


def hash_kmers_ToSet(kmers):
    hashes = set()
    for kmer in kmers:
        hashes.add(hash_kmer(kmer))
    return hashes


def read_kmers_from_file_CreateDicts_Set(filename, ksize):
    all_kmers_Dict = {}
    all_hashes_Set_Dict = {}
    seqLen_Dict = {}
    
    NumParsedRecords = 0
    
    for record in screed.open(filename):

        NumParsedRecords += 1
        sequence = record.sequence

        kmers = build_kmers(sequence, ksize)
        hashes_Set = hash_kmers_ToSet(kmers)
        
        #all_kmers_Dict[record.name] = kmers
        all_hashes_Set_Dict[record.name] = hashes_Set
        seqLen_Dict[record.name] = len(sequence)

    print(NumParsedRecords, " total records were parsed")
    
    return all_kmers_Dict, all_hashes_Set_Dict, seqLen_Dict


def read_kmers_from_file_ToHashesDict(filename, ksize):

    all_hashes_Set_Dict = {}
    seqLen_Dict = {}
    
    NumParsedRecords = 0
    
    for record in screed.open(filename):

        NumParsedRecords += 1
        sequence = record.sequence

        kmers = build_kmers(sequence, ksize)
        hashes_Set = hash_kmers_ToSet(kmers)
        
        all_hashes_Set_Dict[record.name] = hashes_Set
        seqLen_Dict[record.name] = len(sequence)

    print(NumParsedRecords, " total records were parsed")
    
    return all_hashes_Set_Dict, seqLen_Dict





def jaccard_containment_FromSets(a, b):
    '''
    This function returns the Jaccard Containment between sets a and b.
    '''
    
    intersection = len(a.intersection(b))
    
    return intersection / len(a)

def jaccard_similarity_FromSets(a, b):
    '''
    This function returns the Jaccard Similarity between sets a and b.
    '''
    intersection = len(a.intersection(b))
    union = len(a.union(b))
    
    return intersection / union

def jaccard_containment_MaxVal_FromSets(a, b):
    '''
    This function returns the maximum possible Jaccard Containment between sets a and b.
    '''
    
    intersection = len(a.intersection(b))

    min_Len = min(len(a), len(b) )

    return intersection / min_Len


def all_vs_all_kmer_JC(all_SeqIDs, dictOf_Hashes_Set, dictOf_SeqLen):

    listOfTuples = []
    
    for record_Name_1 in all_SeqIDs:
        for record_Name_2 in all_SeqIDs: 
            if record_Name_1 != record_Name_2:
        
                record_1and2_JC = jaccard_containment_FromSets(dictOf_Hashes_Set[record_Name_1],
                                                                  dictOf_Hashes_Set[record_Name_2] )


                if record_1and2_JC != 0:
                    seqlen_1 = dictOf_SeqLen[record_Name_1]
                    seqlen_2 = dictOf_SeqLen[record_Name_2]

                    listOfTuples.append((record_Name_1, record_Name_2, seqlen_1, seqlen_2, record_1and2_JC) )
    if listOfTuples != []:
        PG_AvA_DF = pd.DataFrame(listOfTuples)
        PG_AvA_DF.columns = ["RecordID_1", "RecordID_2", "Record1_Len", "Record2_Len", "JaccardContain"]
        PG_AvA_DF = PG_AvA_DF.sort_values(["JaccardContain", "RecordID_1", "RecordID_2"], ascending=False)
    else:
        PG_AvA_DF = pd.DataFrame(columns=["RecordID_1", "RecordID_2", "Record1_Len", "Record2_Len", "JaccardContain"])

    return PG_AvA_DF



def all_vs_all_kmer_MaxJC(all_SeqIDs, dictOf_Hashes_Set, dictOf_SeqLen):

    listOfTuples = []
    
    for i, record_Name_1 in enumerate(all_SeqIDs) :
        for j, record_Name_2 in enumerate(all_SeqIDs) : 
            if i < j: # Check the seqID index so that the same pair of sequences is not compared twice

                record_1and2_JC = jaccard_containment_MaxVal_FromSets(dictOf_Hashes_Set[record_Name_1],
                                                                      dictOf_Hashes_Set[record_Name_2] )

                if record_1and2_JC != 0:
                    seqlen_1 = dictOf_SeqLen[record_Name_1]
                    seqlen_2 = dictOf_SeqLen[record_Name_2]    

                    listOfTuples.append((record_Name_1, record_Name_2, seqlen_1, seqlen_2, record_1and2_JC) )
    
    if listOfTuples != []:
        PG_AvA_DF = pd.DataFrame(listOfTuples)
        PG_AvA_DF.columns = ["RecordID_1", "RecordID_2", "Record1_Len", "Record2_Len", "MaxJaccardContain"]
        PG_AvA_DF = PG_AvA_DF.sort_values(["MaxJaccardContain", "RecordID_1", "RecordID_2"], ascending=False)
    else:
        PG_AvA_DF = pd.DataFrame(columns=["RecordID_1", "RecordID_2", "Record1_Len", "Record2_Len", "MaxJaccardContain"])

    return PG_AvA_DF





###### End of Functions ######

def main():
    parser = argparse.ArgumentParser(description="Calculate the maximum Jaccard Containment of k-mer content between all pairs of input nucleotide sequences.")
    parser.add_argument('-i', '--in_pg_ref', type=str, required=True, help="Input Panaroo pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo (FASTA)")
    parser.add_argument('-o', '--out_ava_tsv',type=str, required=True, help="All vs all comparison of sequence k-mer profiles using the maximum Jaccard Containment metric. (TSV)")
    parser.add_argument('-k', '--kmer_size',type=int, default=31, help="k-mer size (bp) to use for generating profile of unique k-mers for each sequence")
    #parser.add_argument('-t', '--threads',type=int, default=1, help="Number of threads to use")


    args = parser.parse_args()

    ## 1) Set input parameters and PATHs ####

    input_PG_Ref_FA = args.in_pg_ref

    output_AvA_TSV = args.out_ava_tsv

    kmer_size = args.kmer_size

    #n_threads = args.threads
    
    ## 2) Parse and hash all k-mers for each representative nucleotide sequence
    start = time.time()

    Ref_DictOf_Hashes, Ref_DictOf_SeqLen = read_kmers_from_file_ToHashesDict(input_PG_Ref_FA, kmer_size)             

    All_SeqIDs = list(Ref_DictOf_Hashes.keys())

    end = time.time()
    time_diff = end - start
    print("Time to parse and hash all k-mers:", round(time_diff, 2), "seconds")


    ## 3) Calculate the maximum Jaccard Containment (JC) between all pairs of sequences.

    ### NOTE: The maximum JC between sets a and b will always be symetrical, while JC is not

    # print(f"Performing all vs all comparison with {n_threads} threads.")

    start = time.time()

    PG_AvA_DF = all_vs_all_kmer_MaxJC(All_SeqIDs, Ref_DictOf_Hashes, Ref_DictOf_SeqLen) 

    end = time.time()
    time_diff = end - start
    print("Time to perform all vs all comparison of Max Jaccard Containment of k-mer profiles:", round(time_diff, 2), "seconds")
    print("")

    ## 4) Output the AvA comparison table
    print("Saving the All vs All comparison table (TSV) to:", output_AvA_TSV)
    PG_AvA_DF.to_csv(output_AvA_TSV, sep = "\t", index = False)



if __name__ == "__main__":
    sys.exit(main())



