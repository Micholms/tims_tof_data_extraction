import pandas as pd
import matplotlib.pyplot as plt
from functions.match_functions import *
from functions.preprocess_functions import *
import seaborn as sns



import argparse
parser = argparse.ArgumentParser()

#-db DATABASE -u USERNAME -p PASSWORD -size 20
parser.add_argument("-i", "--input", help="Raw data")
parser.add_argument("-t", "--targets", help="Target values with metadata")
parser.add_argument("-o", "--output", help="Output file name/path")

args = parser.parse_args()

df_restapi=pd.read_csv(args.input,index_col=0)
values_df=pd.read_csv(args.targets, sep=";")
df_restapi.columns=['CCS', 'Precursor_type', 'PrecursorMZ', 'POLARITY',
      'CHARGE', 'RT', 'peaks', 'Well']



match_all=[]
ppm=5
not_found_all=[]
i_match_all=[]
i_not_found_all=[]
data_combined=pd.DataFrame(columns=df_restapi.columns)
t=0
for well in df_restapi["Well"].unique():
    match,not_found,i_match,i_not_found,data_out=match_one_well(df_restapi,well,ppm,values_df,int_s=True)
    for m in match:
        match_all.append(m)
    for n in not_found:
        not_found_all.append(n)
    for im in i_match:
        i_match_all.append(im)
    for i in i_not_found:
        i_not_found_all.append(i)
        
    t=t+len(values_df[values_df["Position"]==well])
    data_combined=pd.concat([data_combined,data_out])
    
data_filtered=data_combined.dropna(subset="Target")#.iloc[:,[-1,0,1,2,3,7,-4,-3,-2]]
print("Found: ", round(len(data_filtered)/t*100,3), "%")
data_filtered.columns=['CCS', 'Precursor_type', 'PrecursorMZ', 'POLARITY', 'CHARGE', 'RT',
       'peaks', 'Well', 'M', 'Target', 'SMILES', 'InChIKey', 'Name']

data_filtered.to_csv(args.output+str(ppm)+"ppm.csv")