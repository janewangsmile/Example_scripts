#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:01:15 2019

@author: jun
"""
import os
os.chdir("/home/jun/Desktop/lncRNAs_ASD/")
# paths
import pandas as pd
transcript = pd.read_csv("Dataset/feature_full/1_exp_seq_with_category3_genes_filtered_low_genes.csv",header=0)
# lncRNAs:
transcript = pd.read_csv("raw_dataset_refinement/refined_input_data/lncRNAs_sequence_data.csv",header=0)
transcript.columns = ['seq_name','longest_ENSG']
# DE lncRNAs:
transcript = pd.read_csv("Dataset/DE_ASD_lncRNAs/lncRNAs_cerebellum_expression_sequence.csv",header=0)

#transcript['Gene.Type'] = [1 if x =='ASD' else 0 for x in transcript['Gene.Type']]
# 
gene_list = ['A', 'C', 'G', 'T']

def backtracking(k, comb, res):
    """
    Use backtracking algorithm to find all possible combinations of length k.
    
    Parameters
    ---------------
    k: int
        Length of combination.
    comb: string
        Current combination.
    res: list
        The list that contains all possible combinations.

    Returns
    ---------------
    None
    """
    if k == 0:
        res.append(comb)
        return

    for gene in gene_list:
        comb += gene
        backtracking(k-1, comb, res)
        comb = comb[:-1]


def get_k_mer_list(k):
    """
    Get the k-mer list with a given k.

    Parameters
    ---------------
    k: int
        Length of combination.

    Returns
    ---------------
    res: list
        The list that contains all possible combinations of length k.
    """
    res = []
    backtracking(k, "", res)
    return res

#################################################
# Extract k-mer features
#################################################

# get 4/3/2/1-mer list
k_mer_list =[]
for k in [1,2,3,4]:
    k_mer_list += get_k_mer_list(k)
    
import regex as re
# translate the gene transcript to k-mer features
k_mer_features = []
for line in transcript['longest_ENSG']:
    feature = [len(re.findall(x,line,overlapped=True))/len(line) for x in k_mer_list]
    k_mer_features.append(feature)

k_mer_features = pd.DataFrame(k_mer_features)
k_mer_features.columns=k_mer_list
# combine the kmers with the transcript
transcript = pd.concat([transcript,k_mer_features],axis=1,sort=False)
transcript = transcript.drop(['longest_ENSG'],axis=1)
transcript.to_csv("Dataset/feature_full/ASD_with_cate3_transcript_seqs_with_kmers.csv",index=False)
transcript.to_csv("raw_dataset_refinement/refined_input_data/lncRNAs_transcript_kmers.csv",index=False)
transcript.to_csv("Dataset/DE_ASD_lncRNAs/lncRNAs_cortex_expression_kmers.csv",index=False)
transcript.to_csv("Dataset/DE_ASD_lncRNAs/lncRNAs_cerebellum_expression_kmers.csv",index=False)




#################################################
# Extract k-mer features for hypothetic loci sets
#################################################
N = [101,201,401]
for n in N:
    files = os.listdir("hypothetic_loci/"+str(n)+"/testSet/")
    for f in files:
        myData = pd.read_csv("hypothetic_loci/"+str(n)+"/testSet/"+f, header = 0)
        k_mer_features = []
        for line in myData['longest_ENSG']:
            feature = [len(re.findall(x,line,overlapped=True))/len(line) for x in k_mer_list]
            k_mer_features.append(feature)
            
        k_mer_features = pd.DataFrame(k_mer_features)
        k_mer_features.columns=k_mer_list
        # combine the kmers with the transcript
        myData = pd.concat([myData,k_mer_features],axis=1,sort=False)
        myData.to_csv("hypothetic_loci/"+str(n)+"/testSet/"+f, index=False)







