import sys
print(f'Working with Python {sys.version}')
import numpy as np
import pandas as pd
from functions.match_functions import find_closest, find_match, count_matches, match_one_well
import fiora.IO.mgfReader as mgfReader
import argparse

def read_rawfile(input_df):
    df_mgf = mgfReader.read(input_df)
    df = pd.DataFrame(df_mgf)
    return df


parser = argparse.ArgumentParser()

#-db DATABASE -u USERNAME -p PASSWORD -size 20
parser.add_argument("-t", "--target", help="Target values with metadata")
parser.add_argument("-d", "--dir", help="Path to folder w spectrum files (csv, mgf)")
parser.add_argument("-w", "--well", help="Well number, ex A01")
#parser.add_argument("-o", "--output", help="Output file name/path")
parser.add_argument("-p", "--ppm", help="PPM for tolerance level", type=float)
args = parser.parse_args()


if __name__ == '__main__':
    ppm=args.ppm
    msms_csv=pd.read_csv(args.dir+args.well+"/"+args.well+".msmsonly.csv")
    values_df=  pd.read_csv(args.target ,sep=";")
    
    spectrum_mgf=args.dir+args.well+"/"+args.well+".gnps.mgf"
# All data
   
    target_data=args.well
    

    print("Targeting", target_data)

    #Sort out target
    target_list=values_df[values_df["Position"]==args.well]

    #Get list of found values in MetaboScape
    data_in=msms_csv
    

    match,not_found,data_out=count_matches(data_in,target_list,ppm)
    matched_msms=data_out.dropna(subset="Target").sort_values("Target")
    matched_msms=matched_msms.dropna(axis=1, how="all")


    
    test=read_rawfile(spectrum_mgf)  #"../../Tims_Tof_data/Tims_Dummy_test/E01.gnps.mgf")
    test["FEATURE_ID"]=test["FEATURE_ID"].astype(int)
    
    df=test.merge(matched_msms, right_on="FEATURE_ID", left_on="FEATURE_ID", how="right")
    df.columns=['SCANS', 'FEATURE_ID', 'PrecursorMZ', 'MSLEVEL', 'CHARGE', 'POLARITY',
       'RTINMINUTES', 'Precursor_type', 'peaks', 'SHARED_NAME', 'RT', 'PrecursorM', 'CCS','ADDUCT',
       'MaxIntensity', '5_P2-E-1_1_9', 'Target', 'SMILES','InChIKey',"Name"]

    output_name=args.well+"_"+str(args.ppm)+"_matched.csv"#output#../dummy_E01_test.csv"
    print("Saving as", output_name)
    df.to_csv(output_name)












    