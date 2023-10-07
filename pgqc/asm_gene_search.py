# pgqc/asm_gene_search.py

import pandas as pd
import tqdm as tqdm

import mappy as mp


#### Define function for parsing Mappy alignment hits

def parse_AlnHits_To_DF(i_AsmAlner, QuerySeq):

    hits = list(i_AsmAlner.map(QuerySeq, cs = True))
    
    i_QueryLen = len(QuerySeq)
    
    listOfAlnRows = []

    if hits != []:
        
        for hit in hits:
        
            SeqID = hit.mlen / hit.blen
            LenQueryAligned = hit.q_en - hit.q_st
            
            QueryCov = LenQueryAligned / i_QueryLen
            
            i_AlnRow = [hit.ctg, hit.r_st, hit.r_en, i_QueryLen, LenQueryAligned,
                        QueryCov,  hit.mlen, hit.blen, SeqID, hit.is_primary]
            listOfAlnRows.append(i_AlnRow)
            
    Aln_DF = pd.DataFrame(listOfAlnRows, columns = ["RefSeq", "Ref_Start", "Ref_End", "QueryLen", "LenQueryAligned", "QueryCoverage", "AlnMatches", "AlnLen", "AlnSeqID", "IsPriAln"])

    return Aln_DF

######################################################



# General function for searching for gene DNA seqs in each assembly

def PresAbsQC_CheckAsmForGeneSeq(i_Gene_PresAbs_DF, i_PG_Ref_NucSeqs,
                                 i_AsmFA_Dict, i_SampleIDs):

                                     
    i_Gene_PresAbs_DF_Updated = i_Gene_PresAbs_DF.copy().set_index("Gene", drop=False)
    
    for i_SampleID in tqdm(i_SampleIDs) :

        i_Alner_Asm = mp.Aligner(i_AsmFA_Dict[i_SampleID], preset="asm10")  # load or build index
        if not i_Alner_Asm: raise Exception(f"ERROR: failed to load/build index for SR Asm - {i_SampleID}")
        
        i_Sample_GenePres = i_Gene_PresAbs_DF.set_index("Gene")[i_SampleID]
        
        i_AbsentGenes = list( i_Sample_GenePres[i_Sample_GenePres == 0].index)
        
        for gene in i_AbsentGenes:
        
            gene_seq = i_PG_Ref_NucSeqs[gene]
        
            i_AlnToAsm_DF = parse_AlnHits_To_DF(i_Alner_Asm, gene_seq)

            i_AlnToAsm_Filt = i_AlnToAsm_DF[(i_AlnToAsm_DF["QueryCoverage"] > 0.9) & (i_AlnToAsm_DF["AlnSeqID"] > 0.9)]
            
            Num_Hits_Pass = i_AlnToAsm_Filt.shape[0]
    
            if Num_Hits_Pass > 0: # Update Pres/Abs matrix to have 2, meaning that the Protein-level annotation is Not present, BUT the gene sequence can be found with high confidence.
                i_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 2

    i_Gene_PresAbs_DF_Updated["NumAsm_WiGene_DNASeq"] = i_Gene_PresAbs_DF_Updated[ListOf_SampleID_Cols].applymap(lambda x: 1 if x > 0 else 0).sum(axis = 1)

    return i_Gene_PresAbs_DF_Updated




def SRAsm_PresAbsQC_CheckInLRAsm(i_SR_Gene_PresAbs_DF, i_SR_PG_Ref_NucSeqs,
                                 i_SR_AsmFA_Dict, i_LR_AsmFA_Dict, i_SampleIDs,
                                 MinQueryCov = 0.9, MinQuerySeqID = 0.9):


    TotalNum_All_SR_MissingGenes = 0
    TotalNum_All_Abs_InSRandLR = 0
    TotalNum_All_Abs_InSR_NotInLR = 0
    
    i_SR_Gene_PresAbs_DF_Updated = i_SR_Gene_PresAbs_DF.copy().set_index("Gene", drop=False)
    
    PresAbs_Updates = []
    for i_SampleID in tqdm(i_SampleIDs) :
    
        i_SRAsm_GenePres = i_SR_Gene_PresAbs_DF.set_index("Gene")[i_SampleID]
        
        i_Alner_SRAsm = mp.Aligner(i_SR_AsmFA_Dict[i_SampleID], preset="asm10")  # load or build index
        if not i_Alner_SRAsm: raise Exception(f"ERROR: failed to load/build index for SR Asm - {i_SampleID}")
        
        i_Alner_LRAsm = mp.Aligner(i_LR_AsmFA_Dict[i_SampleID], preset="asm10")  # load or build index
        if not i_Alner_LRAsm: raise Exception(f"ERROR: failed to load/build index for LR Asm - {i_SampleID}")
        
        i_SRAsm_AbsentGenes = list( i_SRAsm_GenePres[i_SRAsm_GenePres == 0].index)
        
        NumSRAsm_NotInSR_InLR = 0
        NumAbsent_ButInSRandLR = 0
        NumAbsent_InSR_NotInLR = 0

        
        for gene in i_SRAsm_AbsentGenes:
        
            gene_seq = i_SR_PG_Ref_NucSeqs[gene]
        
            i_AlnToSRAsm_DF = parse_AlnHits_To_DF(i_Alner_SRAsm, gene_seq)
            i_AlnToLRAsm_DF = parse_AlnHits_To_DF(i_Alner_LRAsm, gene_seq)
        
            #i_AlnToSRAsm_Filt = i_AlnToSRAsm_DF[(i_AlnToSRAsm_DF["QueryCoverage"] > 0.9) & (i_AlnToSRAsm_DF["AlnSeqID"] > 0.9)]
            #i_AlnToLRAsm_Filt = i_AlnToLRAsm_DF[(i_AlnToLRAsm_DF["QueryCoverage"] > 0.9) & (i_AlnToLRAsm_DF["AlnSeqID"] > 0.9)]

                
            Num_SR_Hits = i_AlnToSRAsm_DF[(i_AlnToSRAsm_DF["QueryCoverage"] > 0.9) & (i_AlnToSRAsm_DF["AlnSeqID"] > 0.9)].shape[0]
            Num_LR_Hits = i_AlnToLRAsm_DF[(i_AlnToLRAsm_DF["QueryCoverage"] > 0.9) & (i_AlnToLRAsm_DF["AlnSeqID"] > 0.9)].shape[0]
    
            #if Num_SR_Hits > 0: # Update Pres/Abs matrix to have 2, meaning that the Protein-level annotation is Not present, BUT the gene sequence can be found with high confidence.
            #    i_SR_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 2
    
            
            if (Num_SR_Hits == 0) & (Num_LR_Hits > 0):
                #i_SR_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 3 # Means Not in SR, In LR Asm

                PresAbs_Updates.append((gene, i_SampleID, 3))
                
                NumSRAsm_NotInSR_InLR += 1
                
            if (Num_SR_Hits > 0) & (Num_LR_Hits == 0):  
                #i_SR_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 4 # Means In SR Asm, NOT In LR Asm
                PresAbs_Updates.append((gene, i_SampleID, 4))
                
                NumAbsent_InSR_NotInLR += 1
    
            if (Num_SR_Hits > 0) & (Num_LR_Hits > 0):
                #i_SR_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 5 # Means In SR Asm, In LR Asm
                PresAbs_Updates.append((gene, i_SampleID, 5))
                
                NumAbsent_ButInSRandLR += 1
            
        TotalNum_All_SR_MissingGenes += NumSRAsm_NotInSR_InLR
    
        TotalNum_All_Abs_InSR_NotInLR += NumAbsent_InSR_NotInLR
        
        TotalNum_All_Abs_InSRandLR += NumAbsent_ButInSRandLR
        
        #print("------------")

    # Apply all updates
    for gene, sample, value in PresAbs_Updates:
        i_SR_Gene_PresAbs_DF_Updated.loc[gene, sample] = value
            
    print("Across all samples, total missing genes - Not in SR Asm, but in LR Asm:", TotalNum_All_SR_MissingGenes)
    print("Across all samples, total missing genes - In SR Asm, Not in LR Asm:", TotalNum_All_Abs_InSR_NotInLR)
    print("Across all samples, total missing genes - In BOTH SR Asm and LR Asm:", TotalNum_All_Abs_InSRandLR)
    
    i_SR_Gene_PresAbs_DF_Updated["NumAsm_WiGene_DNASeq"] = i_SR_Gene_PresAbs_DF_Updated[ListOf_SampleID_Cols].applymap(lambda x: 1 if x in [1, 4, 5] else 0).sum(axis = 1)
    
    return i_SR_Gene_PresAbs_DF_Updated






def get_SRAsm_Vs_LRAsm_QCStats(i_Gene_PresAbs_DF, i_SampleIDs, print_stats = True):

    N_AbsentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] == 0).sum().sum()  
    N_PresentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] == 1).sum().sum()  
    
    N_AbsentCDS_DNASeq_NotInSR_InLR = (i_Gene_PresAbs_DF[i_SampleIDs] == 3).sum().sum()  
    
    N_AbsentCDS_DNASeq_InSR_NotInLR = (i_Gene_PresAbs_DF[i_SampleIDs] == 4).sum().sum()  
    
    N_AbsentCDS_DNASeq_InSR_InLR = (i_Gene_PresAbs_DF[i_SampleIDs] == 5).sum().sum()  

    N_AnyType_AbsentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] != 1).sum().sum()  

    
    if print_stats:
        print("# of absent genes (CDS and DNA level):", N_AbsentCDS, round(N_AbsentCDS/ N_AnyType_AbsentCDS, 4))
        print("# of present genes (CDS):", N_PresentCDS)
    
        print("# of absent genes (DNA only found in LR Asm):", 
              N_AbsentCDS_DNASeq_NotInSR_InLR,
              round(N_AbsentCDS_DNASeq_NotInSR_InLR/ N_AnyType_AbsentCDS, 4))
    
        print("# of absent genes (DNA only found in SR Asm):", N_AbsentCDS_DNASeq_InSR_NotInLR, round(N_AbsentCDS_DNASeq_InSR_NotInLR/ N_AnyType_AbsentCDS, 4))
    
        print("# of absent genes (DNA only found in BOTH SR & LR Asm):",
              N_AbsentCDS_DNASeq_InSR_InLR,
              round(N_AbsentCDS_DNASeq_InSR_InLR / N_AnyType_AbsentCDS, 4))

        print("# of TOTAL absent genes (At Any Level):", N_AnyType_AbsentCDS)

    
    return N_AbsentCDS, N_PresentCDS, N_AbsentCDS_DNASeq_NotInSR_InLR, N_AbsentCDS_DNASeq_InSR_NotInLR, N_AbsentCDS_DNASeq_InSR_InLR















