"""
helper function
"""
import numpy as np
import pandas as pd
import seaborn as sns
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
import os

def raw2abun(df):
    # Convert a pandas df frame with raw counts into one with relative abundances of each OTU in each sample
    # Assumes that the OTUs are in the columns and samples are in rows
    df1 = df.copy()
    df2 = df1.sum(axis=1)
    abundances = df1.copy()
    for i in df1.index:
        abundances.loc[i] = df1.loc[i]/df2.loc[i]
    return abundances

def otu2tax(taxanomy=''):
	## 1. Make a dictionary, that maps the otu number to taxanomy.
	otu2tax = {}
	p = {}
	## 2. Read as input the taxanomy file:
	with open(taxanomy) as g:
		lines = g.readlines()
		for line in lines:
			sep_line = line.split('\t')
			otu2tax[sep_line[0]] = sep_line[1][:-1]
			phyloinfo = sep_line[1].split(';')
			phylum = phyloinfo[1]
			if p.has_key(phylum)==False:
				p[phylum]=0
	return otu2tax
	
def otu2g(otu2tax=''):
	otu2g = {}
	for otu in otu2tax:
		A = otu2tax[otu].split(';')
		s = ';'
		B = s.join(A[0:6])
		otu2g[otu] = B
	return otu2g

def otu2f(otu2tax=''):
	otu2f = {}
	for otu in otu2tax:
		A = otu2tax[otu].split(';')
		s = ';'
		B = s.join(A[0:5])
		otu2f[otu] = B
	return otu2f
	
def otu2o(otu2tax=''):
	otu2o = {}
	for otu in otu2tax:
		A = otu2tax[otu].split(';')
		s = ';'
		B = s.join(A[0:4])
		otu2o[otu] = B
	return otu2o
	
def otu2c(otu2tax=''):
	otu2c = {}
	for otu in otu2tax:
		A = otu2tax[otu].split(';')
		s = ';'
		B = s.join(A[0:3])
		otu2c[otu] = B
	return otu2c
	
def otu2p(otu2tax=''):
	otu2p = {}
	for otu in otu2tax:
		A = otu2tax[otu].split(';')
		s = ';'
		B = s.join(A[0:2])
		otu2p[otu] = B
	return otu2p	
	
def otu2species(otu_df='',otu2tax=''):
	# To change OTU table into matrix that represent the abundances of species
	## 1. Collect the information of all species from the OTU table
	Species_list={}   
	## 1.1. Go through the OTU numbers in the OTU table, collect all the species names
	OTU_list = otu_df.T.index
	for otu in OTU_list:
		t=otu
		if otu2tax.has_key(str(otu)):
			speciesinfo = otu2tax[str(otu)]
			if Species_list.has_key(speciesinfo)==False:
				Species_list[speciesinfo] = otu_df[otu]
			else:
				Species_list[speciesinfo] = otu_df[otu] + Species_list[speciesinfo]
	## 2. Create dataframe
	Species_df = pd.DataFrame(Species_list)
	## Give a name to the index column, so that later on, when write the df to a file, the matrix is complete
	Species_df.columns.name = 'Species'
	
	## 3. Return dataframe, that have the phylum information
	return Species_df	

## Calculate Shannon Index, and create a new dataframe for that
def SDI(abun_df):
    def ShannonIndex(numList):   ## Calculate Shannon Entropy
        SU = sum(numList)
        SDI = 0.0
        for num in numList:
            freq = float(num)/SU
            if freq>0:
        	    SDI = SDI - freq * np.log(freq)    
        return SDI
    # Calculate shannon entropy for each sample
    SDIs = pd.DataFrame(index=abun_df.index, columns=['SDI'])
    for sample in abun_df.index:
        SDIs.loc[sample, 'SDI'] = ShannonIndex(abun_df.loc[sample])
    return SDIs
    
## For any dataframe, subset based on the meta_dict
def SubSetDataFrame(df='', label='', value='',meta_dict=''):
	# Create a new dataframe, that contains the information of the value we specified
	frame = []
	int_df = df.copy()  ## Intermediate dataframe, not sure if that's necessary... 
	for i in range(len(value)):
		df1 = int_df.loc[int_df.index.isin(meta_dict[label][value[i]])]	 ## find rows, that the value is in values[i]
		frame.append(df1) 
	output_df = pd.concat(frame)        ## Convert a matrix into dataframe...
	return output_df
	


## Read OTU table and metadata files, give a combined file
def read_raw_files(otu_table, metadata):
	df1 = pd.read_csv(otu_table, sep='\t', index_col=0)
	df2 = pd.read_csv(metadata, sep='\t', comment='#', index_col=0)
# 	result = pd.concat([df1,df2],axis=1)
	return df1,df2

## 3. Extract metadata, create a dictionary, with contains many many dictionaries #**********************************************************
def ExtractMetadata(meta_df):
    meta_dict = defaultdict(dict)
    # Description: reads metadata file and returns a dict with {label: {value1: list_of_samples, value2: list_ofsamples}, label2: {value1: list_of_samples...}}
    # Indices of meta_df should be SampleID        
    for sample in meta_df.index:
    # All subsequent columns are metadata labels
    	for label in meta_df.columns:
    		val = meta_df.loc[sample, label]
    		# Translation of following line: meta_dict['DiseaseState]['H'].append('S1')
    		meta_dict[label].setdefault(val, []).append(sample)
    return meta_dict
    
