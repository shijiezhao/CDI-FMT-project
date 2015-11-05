"""
Name: 1.BasicAnalysis.py
@Author: Shijie Zhao, Sep/26/15, Cambridge

Goals:
(1). Preprocessing OTU table into a great format...
(2). Just... can't stand other people's code...
(3). Not to re-invent the wheels ever again...
"""
from Helperfunctions import *


def parse_args():
    ## 1. Parse command line arguments                                    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--otu', help = 'Input OTU Table', default = '', required = True)
    parser.add_argument('-M','--metadata', help = 'Input formated metadata file', default = '', required = True)
    parser.add_argument('-T','--taxanomy', help = 'Input otu-taxanomy file', default = '', required = False)
    ## 2. Assign args
    args = parser.parse_args()
    return args




"""
Single species abundance: underlying labels
"""

#######################################################################################################
## 1. Read the arguments #**********************************************************
args = parse_args()
## 2. raw_df: DF with metadata and otu raw abundances
[otu_df, meta_df] = read_raw_files(otu_table=args.otu, metadata=args.metadata)
meta_dict = ExtractMetadata(meta_df)
## 3. Relative:
abun_df = raw2abun(otu_df)
## 4. SDI
SDIs = SDI(abun_df)
## 5. Species abundance
otu2tax = otu2tax(taxanomy=args.taxanomy)
Species_df = otu2species(otu_df=abun_df,otu2tax=otu2tax)
## 6. Choose samples from a single donor
# Post = SubSetDataFrame(df=Species_df, label='type', value=['post'],meta_dict=meta_dict)
# CD_SD = SubSetDataFrame(df=Post, label='donor_subject', value=['MGH03D'],meta_dict=meta_dict)

"""
All phylogeny level information
"""
otu2g = otu2g(otu2tax=otu2tax)
otu2f = otu2f(otu2tax=otu2tax)
otu2o = otu2o(otu2tax=otu2tax)
otu2c = otu2c(otu2tax=otu2tax)
otu2p = otu2p(otu2tax=otu2tax)
G_df = otu2species(otu_df=abun_df,otu2tax=otu2g)
F_df = otu2species(otu_df=abun_df,otu2tax=otu2f)
O_df = otu2species(otu_df=abun_df,otu2tax=otu2o)
C_df = otu2species(otu_df=abun_df,otu2tax=otu2c)
P_df = otu2species(otu_df=abun_df,otu2tax=otu2p)

Complete_phylo = pd.concat([abun_df,Species_df, G_df, F_df, O_df,C_df,P_df], axis=1, join='inner')

"""
(1). otu_df
(2). meta_df
(3). abun_df
(4). SDIs
(5). Species_df
"""




