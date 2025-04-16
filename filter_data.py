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
args = parser.parse_args()


if __name__ == '__main__':
    
    
    
# All data

    df=pd.read_csv(args.input)#, index_col=0)
    
    df=restore_dict(df)
    
    df=filter_low_intensity(df,0.05) # Remove peaks in MS/MS if intensity is bellow 5% of max for that spectrum
    
    df=filter_type(df)
    df=count_peaks(df)
    df=format_df(df)
    df=convert_info(df)
    df=check_InChIKey(df)
    df=filter_peaks(df, 2,["[M+H]+", "[M-H]-"])
    df=compute_graph(df)
    df=CE_filtering(df)
    df=match_fragments_and_peaks(df)
    df=add_metadata(df)
    df=add_identifiers(df)
    df=precursor_processing(df)
    if 'FEATURE_ID' in df.columns:
        df=df.drop(columns=['SCANS', 'FEATURE_ID','RETENTIONTIME','RTINMINUTES'])
    print("---------- Number of matched MS/MS spectrum: ", len(df), "--------------")

    #output_name=args.well+"_"+str(args.decimals)+"_preprocessed.csv"#output#../dummy_E01_test.csv"
    output_name=args.output
    print("Saving as", output_name)
    df.to_csv(output_name)














    