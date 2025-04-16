
import numpy as np
import pandas as pd
import sys
sys.path.append("../")

def find_closest(numbers, target):
    return min(numbers, key=lambda x: abs(x- target))

def find_match(numbers,target, ppm):
    closest=find_closest(numbers, target)
  
    if abs(closest-target)<target*ppm*10**-6:
        return closest
        
def mz_to_m(adduct):
    mass_dict={"[M+H]+":-1.0078+0.0005,"[M+H+H]2+":2*(-1.0078+0.0005),"[M+Na]+":-22.990+0.0005,"[M+K]+":-39.098+0.0005,
              "[M+H-H2O]+":2*(-1.0078-15.999)+0.0005,"[M+H+H2]3+":3*(-1.0078+0.0005)}
    charge_dict={"[M+H]+":1,"[M+H+H]2+":2,"[M+Na]+":1,"[M+K]+":1,
              "[M+H-H2O]+":1,"[M+H+H2]3+":3}
    
    return mass_dict[adduct],charge_dict[adduct]

def count_matches(data,target_list,ppm):#, m_type="m"):
    data=data.copy()
    match=[]
    not_found=[]
    if "PrecursorMZ" in data.columns:
         print("Converting m/z to m for matching")
         for m in (data.index):
            m_change,c_change=mz_to_m(data.loc[m,'Precursor_type'])
            data.loc[m,'M']=(data.loc[m,'PrecursorMZ']+m_change)*c_change #Turn into M from m/z
    else:
        data['M']=data["PEPMASS"]
    all_masses=data['M']
    target_mass=target_list["Exact Mass"]    
    target_smiles=target_list["SMILES"]
    target_INCHIKEY=target_list["INCHIKEY"]
    target_name=target_list["Sample Name"]
    for i in target_list.index:
        m=find_match(all_masses,target_mass[i],ppm)
        if m:
            match.append(target_mass[i])
            j=data.index[data['M']==m].tolist()
            data.loc[j,"Target"]=float(target_mass[i])
            data.loc[j,"SMILES"]=target_smiles[i]
            data.loc[j,"INCHIKEY"]=target_INCHIKEY[i]
            data.loc[j,"Name"]=target_name[i]
            
           
        else:
            not_found.append(target_mass[i])
    print("Found", len(match), "/", len(target_mass), "matches at", ppm,"ppm")
    return match,not_found,data

def match_one_well(df,well,ppm, values_df,int_s=False):

    target_data=well
    print("----- Targeting", target_data)
    
    #Sort out target
    target_list=values_df[values_df["Position"]==target_data]
    #Get list of found values in MetaboScape
    data_in=df[df["Well"]==target_data]
  
    match,not_found,data_out=count_matches(data_in,target_list,ppm)
    if int_s:
        internal_s=pd.read_csv("first_plate/internal_standards.csv", sep=";")
        print("Internal Standard")
        i_match,i_not_found=count_matches(data_in,internal_s,ppm)[0:-1]
        return match,not_found,i_match,i_not_found,data_out
    else:
        return match,not_found, data_out

def process_col_names(df, cut_off):
    msms=df.copy()
    col=msms.columns
    new=[]
    for name in col:
        if "5_P" in name:
            new_name=name[5:9].replace("-","0").strip("_")
            msms[name] = msms[name].apply(lambda x: x if x > cut_off else 0 )
            if len(new_name)==4:
                if new_name[-1]=="0": 
                    new_name=new_name.replace("0","")
                    new_name=new_name+"0"
                    new.append(new_name)
                else:
                    new.append(new_name.replace("0",""))
            else:
                new.append(new_name)
        else:
            new.append(name)
    
    msms.columns=new
    return msms

def add_id(df):
    c=0
    msms=df.copy()
    msms["Well"]=None
    for i, name in enumerate (msms.index):
        value=list(msms.iloc[i,12:-1])
        if value.count(0) ==len(msms.iloc[i,12:-1])-1:
            df=msms.iloc[i,12:-1]
            msms.loc[name,"Well"]=(str(df.loc[~(df==0)].index.values).strip("[").strip("]")[1:-1])
            c=c+1

        else:
            msms.loc[i,"Well"]="Nan"
    print("Number of mz found in just one sample: ",c )
    msms_filtered=msms[msms["Well"]!="Nan"]
    return msms_filtered