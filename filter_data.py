import sys
print(f'Working with Python {sys.version}')
import numpy as np
import pandas as pd
import importlib
import collections

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import rdkit.Chem.Descriptors as Descriptors
from rdkit.Chem import PandasTools
import torch
seed = 42
torch.manual_seed(seed)
torch.set_printoptions(precision=2, sci_mode=False)

# Deep Learning
import sklearn
import scipy

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Load Modules

from functions.preprocess_functions import *
    
import argparse
parser = argparse.ArgumentParser()

#-db DATABASE -u USERNAME -p PASSWORD -size 20
parser.add_argument("-i", "--input", help="Target values with metadata")
parser.add_argument("-o", "--output", help="Output file name/path")
parser.add_argument("-m", "--mona", help="Output file name/path", default=False)
args = parser.parse_args()


if __name__ == '__main__':
    
    
    
# All data

    df=pd.read_csv(args.input)#, index_col=0)
   
    #if args.mona:
      #  df=df.drop(index=85072)
    print("Step 0 of 8")
    df=restore_dict(df)
    df=count_peaks(df)
    df=format_df(df)
    df=filter_peaks(df, 2)
    print("Step 1 of 8")
    if args.mona:
        df=filter_low_intensity(df,0.025,99)# Remove peaks in MS/MS if intensity is bellow 5% of max for that spectrum
    else: df=filter_low_intensity(df,0.025,5000)
    df=filter_type(df)
    print("Step 2 of 8")
    #df=count_peaks(df)
    df=format_df(df)
    print("Step 3 of 8")
    if args.mona:
        df=convert_info(df,args.mona)
    else: 
        df=convert_info(df)
    
   
    #f=convert_info(df)
    print("Step 4 of 8")
    df=check_InChIKey(df)
    #df=filter_peaks(df, 2,["[M+H]+", "[M-H]-"])
    
    print("Step 5 of 8")
    df=CE_filtering(df,["[M+H]+", "[M-H]-"],args.mona)
    df=compute_graph(df) # takes some time
    
    print("Step 6 of 8")
    df=match_fragments_and_peaks(df)
    print("Step 7 of 8")
    df=add_metadata(df,args.mona)
    df=add_identifiers(df)
    print("Step 8 of 8")
    df=precursor_processing(df)
    if 'FEATURE_ID' in df.columns:
        df=df.drop(columns=['SCANS', 'FEATURE_ID','RETENTIONTIME','RTINMINUTES'])
    print("---------- Number of matched MS/MS spectrum: ", len(df), "--------------")

    #output_name=args.well+"_"+str(args.decimals)+"_preprocessed.csv"#output#../dummy_E01_test.csv"
    output_name=args.output
    print("Saving as", output_name)
    df.to_csv(output_name)














    