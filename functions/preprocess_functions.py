import ast
from rdkit import Chem
import sys
sys.path.append("..")
from os.path import expanduser
home = expanduser("~")
import fiora.IO.mgfReader as mgfReader
import fiora.visualization.spectrum_visualizer as sv
import fiora.IO.molReader as molReader
from fiora.MOL.constants import ADDUCT_WEIGHTS, DEFAULT_MODES
from fiora.MOL.collision_energy import align_CE
from fiora.MOL.Metabolite import Metabolite
from fiora.IO.LibraryLoader import LibraryLoader
from fiora.MOL.FragmentationTree import FragmentationTree 
from fiora.GNN.AtomFeatureEncoder import AtomFeatureEncoder
from fiora.GNN.BondFeatureEncoder import BondFeatureEncoder
from fiora.GNN.SetupFeatureEncoder import SetupFeatureEncoder
import numpy as np

def read_rawfile(input_df):
    df_mgf = mgfReader.read(input_df)
    df = pd.DataFrame(df_mgf)
    return df
    
def filter_type(df):
    target_precursor_type = ["[M+H]+", "[M-H]-", "[M+H-H2O]+", "[M+Na]+"]
    a=len(df)
    df = df[df["Precursor_type"].apply(lambda ptype: ptype in target_precursor_type)]
    b=len(df)
    print("Removing ", a-b, "/", a, "due to precursor out of spec")
    return df

def format_df(df):
   # print("> Formating dataset")
    df_mona=df.copy()
    #df_mona['PrecursorMZ'] = df_mona["PrecursorMZ"].str.replace(',', '.')
    df_mona['PrecursorMZ'] = df_mona["PrecursorMZ"].astype('float')
    df_mona['Num Peaks'] = df_mona["Num Peaks"].astype('int')
    return df_mona

def get_mol(df,ID):
    x = df.loc[ID]
    smiles=x["SMILES"]
    if smiles:        
        mol=Chem.MolFromSmiles(smiles)
    else: 
        mol=np.nan
        smiles=np.nan
    return mol,smiles

def convert_info(df):
   # print("> Fetch SMILES and interpret InChIKeys")
    df_mona=df.copy()
    
    df_mona= df_mona[~df_mona["InChIKey"].isnull()] 
    
    MOL=[]
    SMILES=[]
    for i in df_mona.index:
        mol,smiles=get_mol(df_mona,i)
           
        MOL=MOL+[mol]
        SMILES=SMILES+[smiles]
    df_mona["MOL"] = MOL
    df_mona["SMILES"]=SMILES
    print(f"Successfully interpreted {sum(df_mona['MOL'].notna())} from {df_mona.shape[0]} entries. Dropping the rest.")
    df_mona = df_mona[df_mona['MOL'].notna()]
    df_mona["InChI"] = df_mona["MOL"].apply(Chem.MolToInchi)
    df_mona["K"] = df_mona["MOL"].apply(Chem.MolToInchiKey)
    df_mona["ExactMolWeight"] = df_mona["MOL"].apply(Chem.Descriptors.ExactMolWt) # get exact mol weight from "MOL"
    return df_mona

def check_InChIKey(df_mona):
  #  print("> Checking if interpreted InChIKeys are correct")
    correct_keys = df_mona.apply(lambda x: x["InChIKey"] == x["K"], axis=1) #Check so that computed keys are same as original
    s = "confirmed!" if correct_keys.all() else "not confirmed !! Attention!"
    #print(f"Confirming whether computed and provided InChI-Keys are correct. Result: {s} ({correct_keys.sum()/len(correct_keys):0.2f} correct)")
    half_keys = df_mona.apply(lambda x: x["InChIKey"].split('-')[0] == x["K"].split('-')[0], axis=1)
    s = "confirmed!" if half_keys.all() else "not confirmed !! Attention!"
    print(f"Checking if main layer InChI-Keys are correct. Result: {s} ({half_keys.sum()/len(half_keys):0.3f} correct)")
    
    print("Dropping all other.")
    df_mona["matching_key"] = df_mona.apply(lambda x: x["InChIKey"] == x["K"], axis=1)
    df_mona = df_mona[df_mona["matching_key"]]
    return df_mona
    
def filter_peaks(df,MIN_PEAKS,PRECURSOR_TYPES):
    #print("> Filtering based on number of peaks")
    df_mona=df.copy()
    df_mona = df_mona[df_mona["Num Peaks"] > MIN_PEAKS]
    df_mona["theoretical_precursor_mz"] = df_mona["ExactMolWeight"] + df_mona["Precursor_type"].map(ADDUCT_WEIGHTS) # Calculate theoretical mz
    df_mona = df_mona[df_mona["Precursor_type"].apply(lambda ptype: ptype in PRECURSOR_TYPES)]
    df_mona["precursor_offset"] = df_mona["PrecursorMZ"] - df_mona["theoretical_precursor_mz"] # Caluclate difference from true value
    #print(f"Shape {df_mona.shape}")
    return df_mona
    
def compute_graph(df):
    #print("> Computing molecular structure graph")
    df_mona=df.copy()
    TOLERANCE = 200 * PPM
    df_mona["Metabolite"] = df_mona["SMILES"].apply(Metabolite)
    df_mona["Metabolite"].apply(lambda x: x.create_molecular_structure_graph())
    df_mona["Metabolite"].apply(lambda x: x.compute_graph_attributes())
    df_mona["Metabolite"].apply(lambda x: x.fragment_MOL(depth=1))
    df_mona.apply(lambda x: x["Metabolite"].match_fragments_to_peaks(x["peaks"]["mz"], x["peaks"]["intensity"], tolerance=TOLERANCE), axis=1)
    return df_mona    
    
def CE_filtering(df):
    #print("> Filtering on Collision energy")
    df_mona=df.copy()
    df_mona["Collision_energy"] ="20eV" 
    df_mona["CE"] = df_mona.apply(lambda x: align_CE(x["Collision_energy"], x["theoretical_precursor_mz"]), axis=1) #modules.MOL.collision_energy.align_CE) 
    df_mona["CE_type"] = df_mona["CE"].apply(type)
    df_mona["CE_derived_from_NCE"] = df_mona["Collision_energy"].apply(lambda x: "%" in str(x))
    
    print("Distinguish CE absolute values (eV - float) and normalized CE (in % - str format)")
   # print("Removing all but absolute values")
    df_mona = df_mona[df_mona["CE_type"] == float]
    df_mona = df_mona[~df_mona["CE"].isnull()]

    print(f'Detected {len(df_mona["CE"].unique())} unique collision energies in range from {np.min(df_mona["CE"])} to {max(df_mona["CE"])} eV')
    
    df_mona=compute_graph(df_mona)
    return df_mona


def add_metadata(df):
    #print("> Adding various metadata")
    orbitrap_nametags = ["Orbitrap"]
    qtof_nametags = ["QTOF", "LC-ESI-QTOF", "ESI-QTOF"]
    #df["Instrument_type"]="TimsTof"
    df["Instrument_type"]="Q-TOF"
    df["Instrument_type"] = df["Instrument_type"].apply(lambda x: "HCD" if x in orbitrap_nametags else "Q-TOF" if x in qtof_nametags else x)
    df["RETENTIONTIME"] = np.nan
    df["CCS"] = np.nan
    df["PPM_num"] = 50
    df["ppm_peak_tolerance"] = df["PPM_num"] * PPM
    df["lib"] = "TimsTof"
    df["origin"] = "TimsTof"
    df["Ionization"] = "ESI" #??????
   
    return df

def add_identifiers(df):
    #print("> Assigning unique metabolite identifiers.")
    
    metabolite_id_map = {}
    
    for metabolite in df["Metabolite"]:
        is_new = True
        for id, other in metabolite_id_map.items():
            if metabolite == other:
                metabolite.set_id(id)
                is_new = False
                break
        if is_new:
            new_id = len(metabolite_id_map)
            metabolite.id = new_id
            metabolite_id_map[new_id] = metabolite
    
    df["group_id"] = df["Metabolite"].apply(lambda x: x.get_id())
    df["num_per_group"] = df["group_id"].map(df["group_id"].value_counts())
    
    for i, data in df.iterrows():
        data["Metabolite"].set_loss_weight(1.0 / data["num_per_group"])
    print(f"Found {len(metabolite_id_map)} unique molecular structures.")
    return df
    
def precursor_processing(df):
    
    print("> Processing precursors")
    df["loss_weight"] = df["Metabolite"].apply(lambda x: x.loss_weight)
    df["Precursor_offset"] = df["PrecursorMZ"] - df.apply(lambda x: x["Metabolite"].ExactMolWeight + ADDUCT_WEIGHTS[x["Precursor_type"]], axis=1)
    df["Precursor_abs_error"] = abs(df["Precursor_offset"])
    df["Precursor_rel_error"] = df["Precursor_abs_error"] / df["PrecursorMZ"]
    df["Precursor_ppm_error"] = df["Precursor_abs_error"] / (df["PrecursorMZ"] * PPM)
    print((df["Precursor_ppm_error"] > df["PPM_num"]).sum(), "found with misaligned precursor. Removing these.")
    
    df = df[df["Precursor_ppm_error"] <= df["PPM_num"]]
    return df
def count_peaks(df):
    for i in (df.index):
        df.loc[i,"Num Peaks"]=int(len(df.loc[i,"peaks"]["mz"]))
    return df

def restore_dict(df):
    dict_columns = ["peaks"]#, "summary"]
    for col in dict_columns:
        df[col] = df[col].apply(lambda x: ast.literal_eval(str(x).replace('nan', 'None')))
   
    
    return df
    
def filter_low_intensity(df,level):
    peak_list=[]
    saved=0
    tot=0
    for j,name in enumerate(df.index):
        
        mz=[]
        inte=[]
        tot=tot+len(df.loc[name,"peaks"]["mz"])
        for i,I in enumerate(df.loc[name,"peaks"]["intensity"]):
            max_i=max(df.loc[name,"peaks"]["intensity"])
            
            if I>=level*max_i:
                mz.append(df.loc[name,"peaks"]["mz"][i])
                inte.append(I)
                saved=saved+1
        peaks={"mz":mz,"intensity":inte}
        peak_list.append(peaks)   
    df["peaks"]=peak_list
    print("Dropped", str(round((tot-saved)/tot*100,3)), "% of the peaks due to low intensity")
    return df
    
def match_fragments_and_peaks(df):
    #print(">   Matching fragments to peaks")
    df_mona=df.copy()
    df_mona["peak_matches"] = df_mona["Metabolite"].apply(lambda x: getattr(x, "peak_matches"))
    df_mona["num_peaks_matched"] = df_mona["peak_matches"].apply(len)
    
    def get_match_stats(matches, mode_count={m: 0 for m in DEFAULT_MODES}):
        num_unique, num_conflicts = 0, 0
        for mz, match_data in matches.items():
            ion_modes = match_data["ion_modes"]
            if len(ion_modes) == 1:
                num_unique += 1
            elif len(ion_modes) > 1:
                num_conflicts += 1
            for c in ion_modes:
                mode_count[c[0]] += 1
        return num_unique, num_conflicts, mode_count
    
    df_mona["match_stats"] = df_mona["peak_matches"].apply(lambda x: get_match_stats(x))
    df_mona["num_unique_peaks_matched"] = df_mona.apply(lambda x: x["match_stats"][0], axis=1)
    df_mona["num_conflicts_in_peak_matching"] = df_mona.apply(lambda x: x["match_stats"][1], axis=1)
    df_mona["match_mode_counts"] = df_mona.apply(lambda x: x["match_stats"][2], axis=1)
    u= df_mona["num_unique_peaks_matched"].sum() 
    s= df_mona["num_conflicts_in_peak_matching"].sum() 
    print(f"Total number of uniquely matched peaks: {u} , conflicts found within {s} matches ({100 * s / (u+s):.02f} %))")
    print(f"Total number of conflicting peak to fragment matches: {s}")
    return df_mona 