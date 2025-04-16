import sys
print(f'Working with Python {sys.version}')
import numpy as np
import pandas as pd
from functions.match_functions import *
import functions.mgfReader as mgfReader
import argparse



    
def read_rawfile(input_df):
    df_mgf = mgfReader.read(input_df)
    df = pd.DataFrame(df_mgf)
    return df


parser = argparse.ArgumentParser()

#-db DATABASE -u USERNAME -p PASSWORD -size 20
#parser.add_argument("-d", "--dir", help="Parent directory to spectrum files (csv, mgf)")
#parser.add_argument("-t", "--target", help="Name of target file")
#parser.add_argument("-i", "--input", help="Name of target well(s)")
#parser.add_argument("-w", "--well", help="Well number, ex A01")
#parser.add_argument("-o", "--output", help="Output file name/path")
#parser.add_argument("-p", "--ppm", help="PPM for tolerance level", type=float)
#args = parser.parse_args()

parser.add_argument("-m", "--msmsonly", help="Path to msmsonly.csv file")
parser.add_argument("-t", "--target", help="Path to file with targets")
parser.add_argument("-s", "--spectrum", help="Path to msms spectrum file")
#parser.add_argument("-d", "--dir", help="Parent directory to spectrum files (csv, mgf)")
parser.add_argument("-o", "--output", help="Output file name/path")
parser.add_argument("-p", "--ppm", help="PPM for tolerance level", type=float)
args = parser.parse_args()


if __name__ == '__main__':
    
    ppm=args.ppm
    #msms=pd.read_csv(args.dir+"/Exported_data/"+args.input+"/"+args.input+".msmsonly.csv")
    #values_df= pd.read_csv(args.dir+"/"+args.target ,sep=";")
    #spectrum_mgf=args.dir+"/Exported_data/"+args.input+"/"+args.input+".gnps.mgf"

    msms=pd.read_csv(args.msmsonly)
    values_df= pd.read_csv(args.target,sep=";")
    spectrum_mgf=args.spectrum
    #spectrum_mgf
    #import glob, os
    #os.chdir("/mydir")
    #for file in glob.glob("*.txt"):
      #  print(file)
    #spectrum_mgf=args.dir+args.dir+".gnps.mgf"
    msms=process_col_names(msms,5000)
    if "A01" in args.msmsonly:
        cols=list(range(1,74))
        cols.remove(12)
        cols.remove(13)
        msms=msms.iloc[:,cols]
    msms=add_id(msms)
    data_combined=pd.DataFrame(columns=msms.columns)
    
    for well in msms["Well"].unique():
        
        data_out=match_one_well(msms,well,ppm,values_df)[-1]
        
        data_combined=pd.concat([data_combined,data_out])
    
    data_filtered=data_combined.dropna(subset="Target").iloc[:,[-1,0,1,2,3,7,-4,-3,-2]]
    print("Found ", len(data_filtered))
    #spectrum=read_rawfile("first_plate/Exported_data/A01_B23_no_recursive_fitting/A01-B23.gnps.mgf")
    spectrum=read_rawfile(spectrum_mgf)
    spectrum["FEATURE_ID"]=spectrum["FEATURE_ID"].astype(int)
    
    df=spectrum.merge(data_filtered, right_on="FEATURE_ID", left_on="FEATURE_ID", how="right")
    df.columns=['SCANS', 'FEATURE_ID', 'PrecursorMZ', 'MSLEVEL', 'CHARGE', 'POLARITY',
       'RTINMINUTES', 'Precursor_type', 'peaks', "Name", 'RT', 'PrecursorM', 'CCS','ADDUCT',
        'Target', 'SMILES','InChIKey']




    



    #output_name=args.well+"_"+str(args.decimals)+"_matched.csv"#output#../dummy_E01_test.csv"
    output_name=args.output
    print("Saving as", output_name)
    df.to_csv(output_name)












    